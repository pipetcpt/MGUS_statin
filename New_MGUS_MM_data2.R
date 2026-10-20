# -- library -- #
library(MatchIt)
library(cobalt)
library(cmprsk)
library(dplyr)
library(survminer)
library(stringr)
library(ggplot2)
library(haven)
library(mstate)
library(nnet)
library(splines)
library(lubridate)

# -- setting -- #
setwd("/vol/userdata14/sta_room462/backup/")
FU_period = 10
MM_landmark_month = 6
landmark_month = 3
landmark_day = 30.42*landmark_month
source("New_MGUS_MM_data_base.R")

select = dplyr::select

level1_ths = 30
level2_ths = 365
p_ths = 365
exposure_year = 1
FU_year = 6
# -- data -- #
# db3$ageg4 <- ifelse(db3$first_d472_age<60,"< 60", ifelse(db3$first_d472_age<70, "60 - 69", ifelse(db3$first_d472_age<80, "70 - 79",">= 80")))
# db3$ageg4 <- factor(db3$ageg4, levels = c("< 60","60 - 69","70 - 79",">= 80"))

# db3$ageg4 <- case_when(db3$first_d472_age<50 ~ "< 50", db3$first_d472_age<60 ~ "50 - 59",
#                        db3$first_d472_age<70 ~ "60 - 69", db3$first_d472_age<80 ~ "70 - 79", T ~ ">= 80")
# db3$ageg4 <- factor(db3$ageg4, levels = c("< 50","50 - 59","60 - 69", "70 - 79", ">= 80"))

db3$ageg4 <- case_when(db3$first_d472_age<60 ~ "< 60", db3$first_d472_age<70 ~ "60 - 69", T ~ ">= 70")
db3$ageg4 <- factor(db3$ageg4, levels = c("< 60","60 - 69", ">= 70"))

db3$ld_yes = as.numeric(db3$mld_yes == 1 | db3$sld_yes == 1)
db3$db_yes = as.numeric(db3$dwoc_yes == 1 | db3$dwcc_yes == 1)

prep_file_list = paste0(list.files(".", pattern = "_prep_"))[c(4,7,9,12,13,15,16)]
prep_dat_list = lapply(prep_file_list, read_sas)
ATC_DDD_dat = read.csv("my_drug_ATC_DDD.csv", stringsAsFactors = F)
ATC_DDD_dat$DDD_rate = ATC_DDD_dat$amount / ATC_DDD_dat$DDD

prep_dat_list1 = 
  lapply(prep_dat_list, function(x){
    res = 
      x %>% 
      select(-MID) %>%
      unique %>%
      left_join(ATC_DDD_dat %>% select(DIV_CD, DDD_rate), by = "DIV_CD") %>%
      mutate(DDD = TOT_USE_QTY_OR_EXEC_FQ * DDD_rate) %>%
      select(JID, RECU_FR_DD, RECU_TO_DD, TOT_INJC_DDCNT_EXEC_FQ, FQ1_MDCT_QTY, DDD) %>%
      mutate(RECU_TO_DD = ifelse(RECU_TO_DD != "", RECU_TO_DD, RECU_FR_DD))
    
    return(res)
  })

prep_dat_list1[[1]]$drug = "ATORVASTATIN"
prep_dat_list1[[2]]$drug = "FLUVASTATIN"
prep_dat_list1[[3]]$drug = "LOVASTATIN"
prep_dat_list1[[4]]$drug = "PITAVASTATIN"
prep_dat_list1[[5]]$drug = "PRAVASTATIN"
prep_dat_list1[[6]]$drug = "ROSUVASTATIN"
prep_dat_list1[[7]]$drug = "SIMVASTATIN"

# prep_dat %>% 
#   mutate(over10 = as.numeric(FQ1_MDCT_QTY >= 10)) %>%
#   group_by(drug) %>%
#   summarize(over10_prop = mean(over10) * 100,
#             over10_DDD_mean = sum(DDD * over10) / sum(over10))

prep_dat = do.call(bind_rows, prep_dat_list1) %>% filter(FQ1_MDCT_QTY < 10)
prep_dat$RECU_FR_DD = as.Date(prep_dat$RECU_FR_DD, format = "%Y%m%d")
prep_dat$RECU_TO_DD = prep_dat$RECU_FR_DD + prep_dat$TOT_INJC_DDCNT_EXEC_FQ
prep_dat$JID = as.numeric(prep_dat$JID)

db3$start_date = db3$first_d472_date-365*exposure_year # start
db3$end_date = db3$first_d472_date + landmark_day # end
db3$JID = as.numeric(db3$JID)

prep_dat_STATIN_all = 
  prep_dat %>%
  group_by(JID, RECU_FR_DD) %>%
  summarize(period = max(TOT_INJC_DDCNT_EXEC_FQ)) %>%
  ungroup %>%
  mutate(RECU_TO_DD = RECU_FR_DD + period) %>%
  inner_join(db3 %>% select(JID, start_date, end_date), by = "JID") %>%
  group_by(JID) %>%
  mutate(next_drug = lead(RECU_FR_DD)) %>%
  ungroup %>%
  mutate(RECU_TO_DD = case_when(RECU_TO_DD >= next_drug ~ next_drug - 1, T ~ RECU_TO_DD)) %>%
  mutate(real_start = case_when(((RECU_FR_DD <= start_date) & (RECU_TO_DD >= start_date)) ~ start_date, T ~ RECU_FR_DD),
         real_end = case_when(((RECU_FR_DD <= end_date) & (RECU_TO_DD >= end_date)) ~ end_date, T ~ RECU_TO_DD),
         period = as.numeric(real_end - real_start),
         real_prop = case_when((real_start > end_date) ~ 0, (real_end < start_date) ~ 0, T ~ 1)) %>%
  group_by(JID) %>% 
  summarize(real_EXP = sum(period * real_prop))

DDD_db3 = db3 %>% select(JID)
# each_drug_set = unique(substr(colnames(db3)[46:66], 5, str_length(colnames(db3)[46:66])-3))
each_drug_set = sapply(prep_dat_list1, function(x) unique(x$drug))
i=1
pb = txtProgressBar(max = length(each_drug_set), style = 3)
for(i in 1:length(each_drug_set)){
  
  temp_drug = each_drug_set[i]
  
  temp_db3 = db3[,c("JID", "start_date", "end_date")]
  
  prep_dat1 = 
    prep_dat %>%
    filter(drug == temp_drug) %>%
    inner_join(temp_db3, by = "JID") %>%
    mutate(real_start = RECU_FR_DD, real_end = RECU_TO_DD, real_prop = 1)
  
  prep_dat1$real_prop[prep_dat1$RECU_FR_DD > prep_dat1$end_date] = 0 # drug after end date
  prep_dat1$real_prop[prep_dat1$RECU_TO_DD < prep_dat1$start_date] = 0 # drug before start date
  
  # drug start before index & drug end after index
  prep_dat1$real_start[prep_dat1$RECU_FR_DD <= prep_dat1$start_date & prep_dat1$RECU_TO_DD >= prep_dat1$start_date] = 
    prep_dat1$start_date[prep_dat1$RECU_FR_DD <= prep_dat1$start_date & prep_dat1$RECU_TO_DD >= prep_dat1$start_date]
  
  # drug start before failure & drug end after failure
  prep_dat1$real_end[prep_dat1$RECU_FR_DD <= prep_dat1$end_date & prep_dat1$RECU_TO_DD >= prep_dat1$end_date] = 
    prep_dat1$end_date[prep_dat1$RECU_FR_DD <= prep_dat1$end_date & prep_dat1$RECU_TO_DD >= prep_dat1$end_date]
  
  prep_dat1$origin_dur = as.numeric(prep_dat1$RECU_TO_DD - prep_dat1$RECU_FR_DD)
  prep_dat1$real_dur = as.numeric(prep_dat1$real_end - prep_dat1$real_start)
  
  prep_dat1$real_prop[prep_dat1$real_prop != 0] = (prep_dat1$real_dur / prep_dat1$origin_dur)[prep_dat1$real_prop != 0]
  
  prep_dat1$real_DDD = prep_dat1$DDD * prep_dat1$real_prop
  prep_dat1$real_PRD = prep_dat1$TOT_INJC_DDCNT_EXEC_FQ * prep_dat1$real_prop
  
  prep_dat2 = prep_dat1 %>% group_by(JID) %>% summarize(real_cDDD = sum(real_DDD), real_EXP = sum(real_PRD))
  colnames(prep_dat2)[2:3] = paste0(c("DDD_", "EXP_"), temp_drug)
  
  DDD_db3 = DDD_db3 %>% left_join(prep_dat2, by = "JID")
  DDD_db3[[colnames(prep_dat2)[2]]][is.na(DDD_db3[[colnames(prep_dat2)[2]]])] = 0
  DDD_db3[[colnames(prep_dat2)[3]]][is.na(DDD_db3[[colnames(prep_dat2)[3]]])] = 0
  
  setTxtProgressBar(pb, value = i)
}

prep_file_list2 = paste0(list.files(".", pattern = "_prep_"))[c(2,3,6,10)]
prep_oth_list = lapply(prep_file_list2, read_sas)

prep_oth_list1 = 
  lapply(prep_oth_list, function(x){
    res = 
      x %>% 
      select(-MID) %>%
      unique %>%
      mutate(DDD = TOT_USE_QTY_OR_EXEC_FQ) %>%
      select(JID, RECU_FR_DD, RECU_TO_DD, TOT_INJC_DDCNT_EXEC_FQ, FQ1_MDCT_QTY, DDD) %>%
      mutate(RECU_TO_DD = ifelse(RECU_TO_DD != "", RECU_TO_DD, RECU_FR_DD))
    
    return(res)
  })

prep_oth_list1[[1]]$drug = "ARB"
prep_oth_list1[[2]]$drug = "ASPIRIN"
prep_oth_list1[[3]]$drug = "EZETIMIBE"
prep_oth_list1[[4]]$drug = "METFORMIN"

prep_oth = do.call(bind_rows, prep_oth_list1) %>% filter(FQ1_MDCT_QTY < 10)

prep_oth$RECU_FR_DD = as.Date(prep_oth$RECU_FR_DD, format = "%Y%m%d")
prep_oth$RECU_TO_DD = prep_oth$RECU_FR_DD + prep_oth$TOT_INJC_DDCNT_EXEC_FQ
prep_oth$JID = as.numeric(prep_oth$JID)

OTH_db3 = db3 %>% select(JID)
each_drug_set = c("ARB", "ASPIRIN", "EZETIMIBE", "METFORMIN")
i=1
pb = txtProgressBar(max = length(each_drug_set), style = 3)
for(i in 1:length(each_drug_set)){
  
  temp_drug = each_drug_set[i]
  
  prep_oth_temp = 
    prep_oth %>%
    filter(drug == temp_drug) %>%
    group_by(JID, RECU_FR_DD) %>%
    summarize(period = max(TOT_INJC_DDCNT_EXEC_FQ)) %>%
    ungroup %>%
    mutate(RECU_TO_DD = RECU_FR_DD + period) %>%
    inner_join(db3 %>% select(JID, start_date, end_date), by = "JID") %>%
    group_by(JID) %>%
    mutate(next_drug = lead(RECU_FR_DD)) %>%
    ungroup %>%
    mutate(RECU_TO_DD = case_when(RECU_TO_DD >= next_drug ~ next_drug - 1, T ~ RECU_TO_DD)) %>%
    mutate(real_start = case_when(((RECU_FR_DD <= start_date) & (RECU_TO_DD >= start_date)) ~ start_date, T ~ RECU_FR_DD),
           real_end = case_when(((RECU_FR_DD <= end_date) & (RECU_TO_DD >= end_date)) ~ end_date, T ~ RECU_TO_DD),
           period = as.numeric(real_end - real_start),
           real_prop = case_when((real_start > end_date) ~ 0, (real_end < start_date) ~ 0, T ~ 1)) %>%
    group_by(JID) %>% 
    summarize(real_EXP = sum(period * real_prop))
  
  colnames(prep_oth_temp)[2] = paste0("EXP_", temp_drug)
  
  OTH_db3 = OTH_db3 %>% left_join(prep_oth_temp, by = "JID")
  OTH_db3[[colnames(prep_oth_temp)[2]]][is.na(OTH_db3[[colnames(prep_oth_temp)[2]]])] = 0
  
  
  setTxtProgressBar(pb, value = i)
}

DDD_db3$DDD_STATIN = rowSums(DDD_db3 %>% select(contains("DDD")))

DDD_db3 = DDD_db3 %>% left_join(prep_dat_STATIN_all %>% rename(EXP_STATIN = real_EXP), by = "JID") %>% 
  mutate(EXP_STATIN = case_when(is.na(EXP_STATIN) ~ 0, T ~ EXP_STATIN))

# period-based
# DDD_db3$STATIN_P_yes = case_when(DDD_db3$EXP_STATIN < level1_ths ~ 0, DDD_db3$EXP_STATIN <= p_ths ~ NA_real_, T ~ 1) %>% factor
# OTH_db3$ARB_P_yes = case_when(OTH_db3$EXP_ARB < level1_ths ~ 0, OTH_db3$EXP_ARB <= p_ths ~ NA_real_, T ~ 2) %>% factor
# OTH_db3$ASPIRIN_P_yes = case_when(OTH_db3$EXP_ASPIRIN < level1_ths ~ 0, OTH_db3$EXP_ASPIRIN <= p_ths ~ NA_real_, T ~ 1) %>% factor
# OTH_db3$EZETIMIBE_P_yes = case_when(OTH_db3$EXP_EZETIMIBE < level1_ths ~ 0, OTH_db3$EXP_EZETIMIBE <= p_ths ~ NA_real_, T ~ 1) %>% factor
# OTH_db3$METFORMIN_P_yes = case_when(OTH_db3$EXP_METFORMIN < level1_ths ~ 0, OTH_db3$EXP_METFORMIN <= p_ths ~ NA_real_, T ~ 1) %>% factor

DDD_db3$STATIN_P_yes = case_when(DDD_db3$EXP_STATIN <= p_ths ~ 0, T ~ 1) %>% factor
OTH_db3$ARB_P_yes = case_when(OTH_db3$EXP_ARB <= p_ths ~ 0, T ~ 1) %>% factor
OTH_db3$ASPIRIN_P_yes = case_when(OTH_db3$EXP_ASPIRIN <= p_ths ~ 0, T ~ 1) %>% factor
OTH_db3$EZETIMIBE_P_yes = case_when(OTH_db3$EXP_EZETIMIBE <= p_ths ~ 0, T ~ 1) %>% factor
OTH_db3$METFORMIN_P_yes = case_when(OTH_db3$EXP_METFORMIN <= p_ths ~ 0, T ~ 1) %>% factor


# DDD-based
DDD_db3$STATIN_D_yes = case_when(DDD_db3$DDD_STATIN <= level1_ths ~ 0, DDD_db3$DDD_STATIN <= level2_ths ~ 1, T ~ 2) %>% factor
DDD_db3$ATORVA_D_yes = case_when(DDD_db3$DDD_ATORVASTATIN <= level1_ths ~ 0, DDD_db3$DDD_ATORVASTATIN <= level2_ths ~ 1, T ~ 2) %>% factor
DDD_db3$ROSUVA_D_yes = case_when(DDD_db3$DDD_ROSUVASTATIN <= level1_ths ~ 0, DDD_db3$DDD_ROSUVASTATIN <= level2_ths ~ 1, T ~ 2) %>% factor

db4 = 
  db3 %>% 
  select(-contains("EXP")) %>% 
  left_join(DDD_db3, by = "JID") %>%
  left_join(OTH_db3, by = "JID") %>%
  mutate(C_EVENT1 = ifelse(C_TIME >= 12*FU_year, 0, C_EVENT) %>% factor,
         C_TIME1 = ifelse(C_TIME >= 12*FU_year, 12*FU_year, C_TIME)) %>%
  as.data.frame


rm(list = ls()[str_detect(ls(), "prep_dat|prep_oth")])
