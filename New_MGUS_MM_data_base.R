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

# -- setting -- #
setwd("/vol/userdata14/sta_room462/backup/")
select = dplyr::select

# -- data -- #
db_mgus <- read_sas("mgus_index_final.sas7bdat")
db_mgus$first_d472_date = as.Date(db_mgus$first_d472_date1, format = "%Y%m%d")

dat_C90 = read_sas("c90_diag_data.sas7bdat")
dat_C90$RECU_FR_DD = as.Date(dat_C90$RECU_FR_DD, format = "%Y%m%d")
double_C90_JID = dat_C90 %>% filter(substr(PRCL_SYM_TP_CD, 1, 4) == "V193") %>% group_by(JID) %>% summarize(nn = n()) %>% filter(nn >= 2) %>% .$JID
dat_C90_sum = dat_C90 %>% group_by(JID) %>% summarize(first_c90_date = min(RECU_FR_DD))

dat_mm_drug = read_sas("spj_t20_mm_drug.sas7bdat")
dat_mm_drug$RECU_FR_DD = as.Date(dat_mm_drug$RECU_FR_DD, format = "%Y%m%d")

dat_mm = 
  dat_mm_drug %>%
  inner_join(dat_C90_sum %>% filter(JID %in% double_C90_JID), by = "JID") %>%
  # filter(RECU_FR_DD >= first_c90_date) %>%
  group_by(JID) %>%
  summarize(first_mm_date = min(RECU_FR_DD))

db_mgus1 = 
  db_mgus %>% 
  filter(first_d472_age >= 19) %>%
  left_join(dat_C90_sum, by = "JID") %>%
  left_join(dat_mm, by = "JID") %>%
  filter(first_d472_date + MM_landmark_month*30.42 <= first_c90_date | is.na(first_c90_date)) # 5869

# make death data
dat_death1 = read_sas("yes_death1.sas7bdat")
dat_death1$death_date1 = as.Date(dat_death1$death_date1, format = "%Y%m%d")
dat_death1$death_date2 = as.Date(dat_death1$death_date2, format = "%Y%m%d")
dat_death1$death_date = as.Date(apply(dat_death1[,2:3], 1, max))

dat_death2 = 
  read_sas("yes_death2.sas7bdat") %>%
  mutate_at(c("last_date1", "last_date2"), function(x) as.Date(x, format = "%Y%m%d"))

dat_death2$last_date = as.Date(apply(dat_death2[,2:3], 1, max))
dat_death2$death_date2 = dat_death2$last_date + 365
dat_death2$death_date2[dat_death2$last_date + 365 > "2022-12-31"] = NA

dat_death = dat_death2 %>% left_join(dat_death1 %>% select(JID, death_date), by = "JID")
dat_death$death_date[is.na(dat_death$death_date)] = dat_death$death_date2[is.na(dat_death$death_date)]

db = db_mgus1 %>% left_join(dat_death %>% select(JID, death_date), by = "JID")

# c90 death 6month new MM
c90_6mo_death = which(is.na(db$first_mm_date) & (db$first_c90_date + 182 > db$death_date))
db$first_mm_date[c90_6mo_death] = db$death_date[c90_6mo_death]

db$censor_date = db$first_d472_date + 365.25*FU_period
db$censor_date[db$censor_date > "2022-12-31"] = "2022-12-31"
db$OBS_END_DD = apply(db %>% select(censor_date, first_mm_date, death_date), 1, min, na.rm = T) %>% as.Date
db$C_EVENT = 0
db$C_EVENT[db$OBS_END_DD == db$death_date] = 2
db$C_EVENT[db$OBS_END_DD == db$first_mm_date] = 1

db$C_TIME = as.numeric(db$OBS_END_DD - db$first_d472_date) / 30.42

db3 = db %>% filter(C_TIME >= MM_landmark_month)

rm(list = ls()[!(ls() %in% c("db3", "landmark_month", "landmark_day", "MM_landmark_month"))])