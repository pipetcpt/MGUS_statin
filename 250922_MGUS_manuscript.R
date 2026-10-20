# -- library -- #
library(MatchIt)
library(cobalt)
library(cmprsk)
library(dplyr)
library(tidyr)
library(survminer)
library(stringr)
library(ggplot2)
library(haven)
library(mstate)
library(nnet)
library(survey)
library(ggsci)
library(ggpubr)
library(tableone)


# -- setting -- #
setwd("/vol/userdata14/sta_room462/backup/")

# -- source -- #
source("New_MGUS_MM_data2.R")
source("MGUS_MM_function.R")
# write.csv(db4 %>% select(JID, first_d472_date), "spj_MGUS_data.csv", row.names = F)

# MGUS_V193 = read_sas("spj_mgus_v193.sas7bdat") %>% 
#   mutate(JID = as.numeric(JID)) %>%
#   inner_join(db4 %>% select(JID, first_d472_date), by = "JID")

# aa = MGUS_V193 %>%
#   mutate(diag_date = as.Date(RECU_FR_DD, format = "%Y%m%d")) %>%
#   filter(diag_date >= first_d472_date, diag_date <= first_d472_date + 30.42*6,
#          substr(PRCL_SYM_TP_CD, 1, 4) == "V193")

# db3_drug = read.csv("spj_MGUS_P1_var.csv")

# -- set variable -- #
x_var = c("SEX_TP_CD", "ageg4", "mi_yes", "chf_yes", "pvd_yes", "cvd_yes", "dem_yes", "cpd_yes",
          "rhe_yes", "pud_yes", "ld_yes", "db_yes", "hp_yes", "rd_yes", "cancer_yes")

x_var2 = c(x_var, c("STATIN_P_yes", "ARB_P_yes", "ASPIRIN_P_yes", "EZETIMIBE_P_yes", "METFORMIN_P_yes"))
x_var3 = c(x_var, c("ARB_P_yes", "ASPIRIN_P_yes", "EZETIMIBE_P_yes", "METFORMIN_P_yes"))
x_var4 = c(x_var, c("ARB_P_yes", "ASPIRIN_P_yes", "METFORMIN_P_yes"))

sum(db4$DDD_ATORVASTATIN > 0) ; mean(db4$DDD_ATORVASTATIN > 0) * 100
sum(db4$DDD_FLUVASTATIN > 0) ; mean(db4$DDD_FLUVASTATIN > 0) * 100
sum(db4$DDD_LOVASTATIN > 0) ; mean(db4$DDD_LOVASTATIN > 0) * 100
sum(db4$DDD_PITAVASTATIN > 0) ; mean(db4$DDD_PITAVASTATIN > 0) * 100
sum(db4$DDD_PRAVASTATIN > 0) ; mean(db4$DDD_PRAVASTATIN > 0) * 100
sum(db4$DDD_ROSUVASTATIN > 0) ; mean(db4$DDD_ROSUVASTATIN > 0) * 100
sum(db4$DDD_SIMVASTATIN > 0) ; mean(db4$DDD_SIMVASTATIN > 0) * 100

# --------------------- #
#### Part 1 analyis  ####
# --------------------- #

#### Table S4: Total baseline characteristics ####
x_var_base = c("C_EVENT1", x_var, "aids_yes", 
               "STATIN_P_yes", "ARB_P_yes", "ASPIRIN_P_yes", "EZETIMIBE_P_yes", "METFORMIN_P_yes")
db_t1 = db4 %>% mutate_at(x_var_base, as.factor)

base_res_overall = CreateTableOne(vars = c("first_d472_age", x_var_base), data = db_t1)

write.csv(print(base_res_overall), "spj_MGUS_TableS4_short.csv")
write.csv(print(base_res_overall, showAllLevels = T), "spj_MGUS_TableS4_long.csv")

survfit(Surv(C_TIME, C_EVENT != 1) ~ 1, data = db4) # median follow-up time
survfit(Surv(C_TIME, C_EVENT != 1) ~ STATIN_D_yes, data = db4) # median follow-up time
survfit(Surv(C_TIME, C_EVENT != 1) ~ ATORVA_D_yes, data = db4) # median follow-up time
survfit(Surv(C_TIME, C_EVENT != 1) ~ ROSUVA_D_yes, data = db4) # median follow-up time

survdiff(Surv(C_TIME1, C_EVENT1 != 1) ~ STATIN_D_yes, data = db4) %>% summary # log-rank
survdiff(Surv(C_TIME1, C_EVENT1 != 1) ~ ATORVA_D_yes, data = db4) %>% summary # log-rank
survdiff(Surv(C_TIME1, C_EVENT1 != 1) ~ ROSUVA_D_yes, data = db4) %>% summary # log-rank


#### Figure S2: Cumulative incidence MM in total / MM subcohort ####
fg_dat_S2 = crprep("C_TIME1", "C_EVENT1", data = db4, trans = 1, cens = 0, id = "JID")
fit_S2A = survfit(Surv(Tstart, Tstop, status == 1) ~ 1, data = fg_dat_S2, weights = weight.cens)
fit_S2B = survfit(Surv(Tstart, Tstop, status == 1) ~ 1, data = fg_dat_S2 %>% filter(status == 1), weights = weight.cens)

fit_S2A$n.risk[1] = nrow(db4)
fit_S2B$n.risk[1] = nrow(db4 %>% filter(C_EVENT1 == 1))

plot_S2A = 
  fit_S2A %>%
  ggsurvplot(risk.table = T, legend = "none", conf.int = T, conf.int.fill = pal_lancet()(2)[2], xlab = "Time (Month)", ylab = "Survival probability", 
             xlim = c(0, 72), break.x.by = 12, risk.table.height = 0.3, fontsize = 4.5,
             legend.title = "", palette = pal_lancet()(2)[2], fun = 'event', title = "(A)")

plot_S2A_add = data.frame(
  x = c(12, 24, 48, 72),
  y = 1 - c(plot_S2A$data.survplot %>% filter(time <= 12) %>% filter(time == max(time)) %>% .$surv,
            plot_S2A$data.survplot %>% filter(time <= 24) %>% filter(time == max(time)) %>% .$surv,
            plot_S2A$data.survplot %>% filter(time <= 48) %>% filter(time == max(time)) %>% .$surv,
            plot_S2A$data.survplot %>% filter(time <= 72) %>% filter(time == max(time)) %>% .$surv),
  y_lci = 1 - c(plot_S2A$data.survplot %>% filter(time <= 12) %>% filter(time == max(time)) %>% .$upper,
                plot_S2A$data.survplot %>% filter(time <= 24) %>% filter(time == max(time)) %>% .$upper,
                plot_S2A$data.survplot %>% filter(time <= 48) %>% filter(time == max(time)) %>% .$upper,
                plot_S2A$data.survplot %>% filter(time <= 72) %>% filter(time == max(time)) %>% .$upper),
  y_uci = 1 - c(plot_S2A$data.survplot %>% filter(time <= 12) %>% filter(time == max(time)) %>% .$lower,
                plot_S2A$data.survplot %>% filter(time <= 24) %>% filter(time == max(time)) %>% .$lower,
                plot_S2A$data.survplot %>% filter(time <= 48) %>% filter(time == max(time)) %>% .$lower,
                plot_S2A$data.survplot %>% filter(time <= 72) %>% filter(time == max(time)) %>% .$lower)
)
plot_S2A_add$l = sprintf("%.1f%% (%.1f%%, %.1f%%)", plot_S2A_add$y * 100, plot_S2A_add$y_lci * 100, plot_S2A_add$y_uci * 100)

plot_S2A$plot = 
  plot_S2A$plot + 
  geom_segment(aes(x = x, xend = x, y = -Inf, yend = y + 0.02), data = plot_S2A_add) +
  scale_y_continuous(breaks = seq(0, 0.10, 0.02), labels = seq(0, 0.10, 0.02)*100) +
  geom_text(aes(x = x-4.2, y = y + 0.025, label = l), data = plot_S2A_add) +
  labs(y = "Cumulative incidence (%)")

ggsave(filename = "spj_MGUS_FigureS2A.pdf", plot = plot_S2A, device = "pdf", width = 7, height = 5.5, dpi = 600)
ggsave(filename = "spj_MGUS_FigureS2A.png", plot = plot_S2A, device = "png", width = 7, height = 5.5, dpi = 600)

plot_S2B = 
  fit_S2B %>%
  ggsurvplot(risk.table = T, legend = "none", conf.int = T, conf.int.fill = pal_lancet()(2)[2],xlab = "Time(Month)", ylab = "Survival probability", 
             xlim = c(0, 72), break.x.by = 12, risk.table.height = 0.3, fontsize = 4.5,
             legend.title = "", palette = pal_lancet()(2)[2], fun = 'event', surv.median.line = "hv", title = "(B)")

med_res = sprintf("%.1f months (%.1f months, %.1f months)", surv_median(fit_S2B)[2], surv_median(fit_S2B)[3], surv_median(fit_S2B)[4])

plot_S2B$plot = 
  plot_S2B$plot + 
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 1, 0.25)*100) +
  geom_text(aes(x = 38, y = 0.35, label = med_res), data = plot_S2A_add) +
  labs(y = "Cumulative incidence (%)")

ggsave(filename = "spj_MGUS_FigureS2B.pdf", plot = plot_S2B, device = "pdf", width = 7, height = 5.5, dpi = 600)
ggsave(filename = "spj_MGUS_FigureS2B.png", plot = plot_S2B, device = "png", width = 7, height = 5.5, dpi = 600)

#### Table S5: Univariate, Multivariate HR ####
trt1_var = paste0(c("STATIN", "ARB", "ASPIRIN", "EZETIMIBE", "METFORMIN"), "_P_yes")
trt1_var = paste0(c("STATIN", "ARB", "ASPIRIN", "METFORMIN"), "_P_yes")
# table(floor(db4$EXP_STATIN/10), db4$C_EVENT1 == 1)
fit_S5_uni <- fit_S5_mul <- fit_S5_mul2 <- list()
i=1
for(i in 1:length(trt1_var)){
  target_trt = trt1_var[i]
  
  form_temp_uni = as.formula(paste0("Surv(C_TIME1, C_EVENT1 == 1) ~ ", target_trt))
  form_exp1 = as.formula(paste0(target_trt, "~ +SEX_TP_CD+ageg4"))
  form_exp2 = as.formula(paste0(target_trt, "~", paste0(c(trt1_var[-i], x_var), collapse = "+")))
  
  fit_S5_uni[[i]] = coxph(form_temp_uni, data = db4)
  
  ww = get_sw(form_exp1, db4)$weight
  fit_S5_mul[[i]] = coxph(form_temp_uni, data = db4, weights = ww)
  
  ww = get_sw(form_exp2, db4)$weight
  fit_S5_mul2[[i]] = coxph(form_temp_uni, data = db4, weights = ww)
}
form_temp_mul2 = as.formula(paste0("Surv(C_TIME1, C_EVENT1 == 1) ~ ", paste0(c(trt1_var[-4], x_var), collapse = "+")))
fit_S5_mul2 = coxph(form_temp_mul2, data = db4)

fit_S5_res = data.frame(var = trt1_var)
fit_S5_res$uni = sapply(fit_S5_uni, function(x) sprintf("%.2f (%.2f, %.2f)", summary(x)$conf.int[,1], summary(x)$conf.int[,3], summary(x)$conf.int[,4]))
fit_S5_res$mul = sapply(fit_S5_mul, function(x) sprintf("%.2f (%.2f, %.2f)", summary(x)$conf.int[,1], summary(x)$conf.int[,3], summary(x)$conf.int[,4]))
fit_S5_res$mul2 = sapply(fit_S5_mul2, function(x) sprintf("%.2f (%.2f, %.2f)", summary(x)$conf.int[,1], summary(x)$conf.int[,3], summary(x)$conf.int[,4]))
fit_S5_res$mul2 = apply(summary(fit_S5_mul2)$conf.int[c(1,2,3,NA,4),], 1, function(x) sprintf("%.2f (%.2f, %.2f)", x[1], x[3], x[4]))
fit_S5_res
write.csv(fit_S5_res, "spj_MGUS_TableS5.csv", row.names = F)

# Table S. Sensitivity analysis of Part I analysis
db4_sens = db4 %>% left_join(db3_drug, by = "JID")

form_temp_sens1 = as.formula(paste0("Surv(C_TIME1, C_EVENT1 == 1) ~ ", 
                                    paste0(c(colnames(db4_sens_temp)[65:69][-4], x_var), collapse = "+")))

form_temp_sens2 = as.formula(paste0("Surv(C_TIME1, C_EVENT1 == 1) ~ ", 
                                    paste0(c(colnames(db4_sens_temp)[70:74][-4], x_var), collapse = "+")))

form_temp_sens3 = as.formula(paste0("Surv(C_TIME1, C_EVENT1 == 1) ~ ", 
                                    paste0(c(colnames(db4_sens_temp)[75:79][-4], x_var), collapse = "+")))

db4_sens_temp04 <- db4_sens_temp06 <- db4_sens_temp08 <- db4_sens
# 40%
day04 = round(c((365+91), (730+91), (1095+91)) * c(0.7)) + 1
db4_sens_temp04[,65:69] = (db4_sens_temp04[,65:69] >= day04[1]) + 0
db4_sens_temp04[,70:74] = (db4_sens_temp04[,70:74] >= day04[2]) + 0
db4_sens_temp04[,75:79] = (db4_sens_temp04[,75:79] >= day04[3]) + 0

# round(mapply(table, db4_sens_temp04[65:79]) / 5131 * 100, 1)
# 60%
day06 = round(c((365+91), (730+91), (1095+91)) * c(0.8)) + 1
db4_sens_temp06[,65:69] = (db4_sens_temp06[,65:69] >= day06[1]) + 0
db4_sens_temp06[,70:74] = (db4_sens_temp06[,70:74] >= day06[2]) + 0
db4_sens_temp06[,75:79] = (db4_sens_temp06[,75:79] >= day06[3]) + 0

# 80%
dat08 = round(c((365+91), (730+91), (1095+91)) * c(0.9)) + 1
db4_sens_temp08[,65:69] = (db4_sens_temp08[,65:69] >= dat08[1]) + 0
db4_sens_temp08[,70:74] = (db4_sens_temp08[,70:74] >= dat08[2]) + 0
db4_sens_temp08[,75:79] = (db4_sens_temp08[,75:79] >= dat08[3]) + 0

fit_S5_sens_04_1 = coxph(form_temp_sens1, data = db4_sens_temp04)
fit_S5_sens_04_2 = coxph(form_temp_sens2, data = db4_sens_temp04)
fit_S5_sens_04_3 = coxph(form_temp_sens3, data = db4_sens_temp04)

fit_S5_sens_06_1 = coxph(form_temp_sens1, data = db4_sens_temp06)
fit_S5_sens_06_2 = coxph(form_temp_sens2, data = db4_sens_temp06)
fit_S5_sens_06_3 = coxph(form_temp_sens3, data = db4_sens_temp06)

fit_S5_sens_08_1 = coxph(form_temp_sens1, data = db4_sens_temp08)
fit_S5_sens_08_2 = coxph(form_temp_sens2, data = db4_sens_temp08)
fit_S5_sens_08_3 = coxph(form_temp_sens3, data = db4_sens_temp08)


res_S7 = 
  rbind(
    sprintf("%.2f (%.2f, %.2f)", summary(fit_S5_sens_04_1)$conf.int[1:4,1],
            summary(fit_S5_sens_04_1)$conf.int[1:4,3], summary(fit_S5_sens_04_1)$conf.int[1:4,4]),
    
    sprintf("%.2f (%.2f, %.2f)", summary(fit_S5_sens_06_1)$conf.int[1:4,1],
            summary(fit_S5_sens_06_1)$conf.int[1:4,3], summary(fit_S5_sens_06_1)$conf.int[1:4,4]),
    
    sprintf("%.2f (%.2f, %.2f)", summary(fit_S5_sens_08_1)$conf.int[1:4,1],
            summary(fit_S5_sens_08_1)$conf.int[1:4,3], summary(fit_S5_sens_08_1)$conf.int[1:4,4]),
    
    sprintf("%.2f (%.2f, %.2f)", summary(fit_S5_sens_04_2)$conf.int[1:4,1],
            summary(fit_S5_sens_04_2)$conf.int[1:4,3], summary(fit_S5_sens_04_2)$conf.int[1:4,4]),
    
    sprintf("%.2f (%.2f, %.2f)", summary(fit_S5_sens_06_2)$conf.int[1:4,1],
            summary(fit_S5_sens_06_2)$conf.int[1:4,3], summary(fit_S5_sens_06_2)$conf.int[1:4,4]),
    
    sprintf("%.2f (%.2f, %.2f)", summary(fit_S5_sens_08_2)$conf.int[1:4,1],
            summary(fit_S5_sens_08_2)$conf.int[1:4,3], summary(fit_S5_sens_08_2)$conf.int[1:4,4]),
    
    sprintf("%.2f (%.2f, %.2f)", summary(fit_S5_sens_04_3)$conf.int[1:4,1],
            summary(fit_S5_sens_04_3)$conf.int[1:4,3], summary(fit_S5_sens_04_3)$conf.int[1:4,4]),
    
    sprintf("%.2f (%.2f, %.2f)", summary(fit_S5_sens_06_3)$conf.int[1:4,1],
            summary(fit_S5_sens_06_3)$conf.int[1:4,3], summary(fit_S5_sens_06_3)$conf.int[1:4,4]),
    
    sprintf("%.2f (%.2f, %.2f)", summary(fit_S5_sens_08_3)$conf.int[1:4,1],
            summary(fit_S5_sens_08_3)$conf.int[1:4,3], summary(fit_S5_sens_08_3)$conf.int[1:4,4])
  )
write.csv(res_S7, "spj_MGUS_Table_S7.csv", row.names = F)
colnames(db4_sens)
# --------------------- #
#### Part 2 analyis  ####
# --------------------- #
# Statin prescription
# db4_statin = db4[db4$DDD_STATIN > 0,]
# round((db4 %>% select(contains("DDD_")) %>% apply(2, function(x) sum(x > 0))) / nrow(db4) * 100, 1)

# setting
DDD_col = pal_lancet()(3)[c(1,3,2)]
DDD_lab = c("Low", "Moderate", "High")
# DDD_level_ind = paste0(c("1. ", "2. ", "3. "), 
#                        c(ifelse(level1_ths == 0, "Non-user", paste0("0-", level1_ths)),
#                          paste0(level1_ths+1, "-", level2_ths),
#                          paste0("> ", level2_ths)), " DDD-days")
res_form_wcif = Surv(C_TIME1, C_EVENT1 == 1) ~ 1

#### Table 1: baseline characteristics ####
x_var_base = c(x_var3, colnames(db4 %>% select(contains("DDD")))[-8])
db_t1 = db4 %>% mutate_at(c("C_EVENT", "STATIN_D_yes", "ATORVA_D_yes", "ROSUVA_D_yes", x_var3), as.factor)

base_res_statin = CreateTableOne(vars = c("C_EVENT1", "first_d472_age", x_var_base), strata = "STATIN_D_yes", data = db_t1)
base_res_atorva = CreateTableOne(vars = c("C_EVENT1", "first_d472_age", x_var_base), strata = "ATORVA_D_yes", data = db_t1)
base_res_rosuva = CreateTableOne(vars = c("C_EVENT1", "first_d472_age", x_var_base), strata = "ROSUVA_D_yes", data = db_t1)

write.csv(print(base_res_statin), "spj_MGUS_Table1_short.csv")
write.csv(print(base_res_statin, showAllLevels = T), "spj_MGUS_Table1_long.csv")

write.csv(print(base_res_atorva), "spj_MGUS_TableS_atorva_short.csv")
write.csv(print(base_res_atorva, showAllLevels = T), "spj_MGUS_TableS_atorva_long.csv")

write.csv(print(base_res_rosuva), "spj_MGUS_TableS_rosuva_short.csv")
write.csv(print(base_res_rosuva, showAllLevels = T), "spj_MGUS_TableS_rosuva_long.csv")

survfit(Surv(C_TIME1, C_EVENT1 != 1) ~ STATIN_D_yes, data = db4) # median follow-up time
survdiff(Surv(C_TIME1, C_EVENT1 != 0) ~ STATIN_D_yes, data = db4)

mean((db4$DDD_ATORVASTATIN > 0))
mean((db4$DDD_ROSUVASTATIN > 0))
mean((db4$DDD_SIMVASTATIN > 0))
mean((db4$DDD_FLUVASTATIN > 0))
mean((db4$DDD_PITAVASTATIN > 0))


#### Table 2: Incidence rate ####
# db4$ROSUVA_D_yes = case_when(db4$DDD_ROSUVASTATIN <= 28 ~ 0, db4$DDD_ROSUVASTATIN <= 365 ~ 1, T ~ 2) %>% factor
db_t2 = db4 %>%
  mutate(cut_STATIN = case_when(STATIN_D_yes == 0 ~ DDD_lab[1], STATIN_D_yes == 1 ~ DDD_lab[2], T ~ DDD_lab[3]),
         cut_ATORVA = case_when(ATORVA_D_yes == 0 ~ DDD_lab[1], ATORVA_D_yes == 1 ~ DDD_lab[2], T ~ DDD_lab[3]),
         cut_ROSUVA = case_when(ROSUVA_D_yes == 0 ~ DDD_lab[1], ROSUVA_D_yes == 1 ~ DDD_lab[2], T ~ DDD_lab[3])) %>%
  mutate_at(c("cut_STATIN", "cut_ATORVA", "cut_ROSUVA"), function(x) factor(x, levels = c("Low", "Moderate", "High")))
table(db4$STATIN_D_yes, db4$C_EVENT1)
table(db4$ATORVA_D_yes, db4$C_EVENT1)
table(db4$ROSUVA_D_yes, db4$C_EVENT1)
IR_res = 
  db_t2 %>% rename(cut = cut_STATIN) %>% group_by(cut) %>% mutate(C_TIME1 = C_TIME1 - 3) %>%
  summarize(nn = n(), n_event = sum(C_EVENT1 == 1), py = sum(C_TIME1)/12) %>%
  mutate(IR = n_event / py * 100000, drug = "STATIN") %>%
  bind_rows(
    db_t2 %>% rename(cut = cut_ATORVA) %>% group_by(cut) %>% mutate(C_TIME1 = C_TIME1 - 3) %>% 
      summarize(nn = n(), n_event = sum(C_EVENT1 == 1), py = sum(C_TIME1)/12) %>%
      mutate(IR = n_event / py * 100000, drug = "ATORVA"),
    db_t2 %>% rename(cut = cut_ROSUVA) %>% group_by(cut) %>% mutate(C_TIME1 = C_TIME1 - 3) %>%
      summarize(nn = n(), n_event = sum(C_EVENT1 == 1), py = sum(C_TIME1)/12) %>%
      mutate(IR = n_event / py * 100000, drug = "ROSUVA")  
  ) 
IR_res$prop = IR_res$n_event / IR_res$nn * 100
IR_res$IR_LCI = qchisq(0.025, df = IR_res$n_event * 2) / 2 / IR_res$py * 100000
IR_res$IR_UCI = qchisq(0.975, df = IR_res$n_event * 2 + 2) / 2 / IR_res$py * 100000
IR_res$IR_res = sprintf("%.1f (%.1f, %.1f)", IR_res$IR, IR_res$IR_LCI, IR_res$IR_UCI)

glm(as.formula(paste0("I(C_EVENT1 == 1) ~ ROSUVA_D_yes + offset(log(C_TIME1))+", paste0(x_var4, collapse = "+"))), data = db4, family = "quasipoisson") %>% summary



t2_form = paste0("Surv(C_TIME1, C_EVENT1 == 1) ~ ", c("cut_STATIN", "cut_ATORVA", "cut_ROSUVA"), "+", paste0(x_var4, collapse = "+"))
t2_fit1 = coxph(as.formula(t2_form[1]), data = db_t2)
t2_fit2 = coxph(as.formula(t2_form[2]), data = db_t2)
t2_fit3 = coxph(as.formula(t2_form[3]), data = db_t2)

IR_res$adj_HR = c(
  1, apply(summary(t2_fit1)$conf.int[1:2,c(1,3,4)], 1, function(x){sprintf("%.2f (%.2f, %.2f)", x[1], x[2], x[3])}),
  1, apply(summary(t2_fit2)$conf.int[1:2,c(1,3,4)], 1, function(x){sprintf("%.2f (%.2f, %.2f)", x[1], x[2], x[3])}),
  1, apply(summary(t2_fit3)$conf.int[1:2,c(1,3,4)], 1, function(x){sprintf("%.2f (%.2f, %.2f)", x[1], x[2], x[3])})
)
IR_res
write.csv(IR_res %>% select(drug, cut, nn, n_event, prop, py, IR, adj_HR), "spj_MGUS_Table2.csv", row.names = F)

#### Table S4. Weighted baseline characteristics + Figure S3. Love plot ####
# table(db4$STATIN_yes, db4$mi_yes)
db4_t2 = db4 %>% mutate_at(colnames(db4)[str_detect(colnames(db4), "yes")], as.factor)
# exp_form_P1 = paste0("STATIN_D_yes ~ ", paste0(x_var3, collapse = "+"), "+ARB_P_yes*EZETIMIBE_P_yes") %>% as.formula
# exp_form_P2 = paste0("ATORVA_D_yes ~ ", paste0(x_var3, collapse = "+"), "+ageg4*ARB_P_yes") %>% as.formula
# exp_form_P3 = paste0("ROSUVA_D_yes ~ ", paste0(x_var3, collapse = "+"), "+ageg4*ARB_P_yes") %>% as.formula

exp_form_P1 = paste0("STATIN_D_yes ~ ", paste0(x_var4, collapse = "+")) %>% as.formula
exp_form_P2 = paste0("ATORVA_D_yes ~ ", paste0(x_var4, collapse = "+")) %>% as.formula
exp_form_P3 = paste0("ROSUVA_D_yes ~ ", paste0(x_var4, collapse = "+"), "+ageg4*ARB_P_yes") %>% as.formula

my_var_lab = c("Sex", "Age", "MI", "CHF", "PVD", "CVD", "Dementia", "CPD", "RHE", "PUD",
               "LD", "Diabetes", "HP", "RD", "Cancer", "ARB", "Aspirin", "Metformin")

DDD_love_lab = c("Low vs Moderate", "Low vs High", "Moderate vs High")

love_P1 = get_loveplot(exp_form_P1, data = db4_t2, var_lab = my_var_lab, title = "Statin", col_lab = DDD_love_lab, trunc_cut = 0.001)
love_P2 = get_loveplot(exp_form_P2, data = db4_t2, var_lab = my_var_lab, title = "Atorvastatin", col_lab = DDD_love_lab, trunc_cut = 0.001)
love_P3 = get_loveplot(exp_form_P3, data = db4_t2, var_lab = my_var_lab, title = "Rosuvastatin", col_lab = DDD_love_lab, trunc_cut = 0.001)

round(ExtractSmd(love_P1$adj_dat), 3) %>% max
round(ExtractSmd(love_P2$adj_dat), 3) %>% max
round(ExtractSmd(love_P3$adj_dat), 3) %>% max

write.csv(print(love_P1$adj_dat), "spj_MGUS_TableS6_short.csv", row.names = F)
write.csv(print(love_P1$adj_dat, showAllLevels = T), "spj_MGUS_TableS6_long.csv", row.names = F)

ggarrange(love_P1$plot, love_P2$plot, love_P3$plot, ncol = 3, nrow = 1, common.legend = T, legend = "bottom")
ggsave("spj_MGUS_FigureS3.png", device = "png", dpi = 300, width = 14, height = 6)
ggsave("spj_MGUS_FigureS3.pdf", device = "pdf", dpi = 300, width = 14, height = 6)

#### Figure 1. Weighted cumulative incidence ####
# db4$STATIN_D_yes = case_when(db4$DDD_STATIN <= 30 ~ 0, db4$DDD_STATIN <= 91 ~ 1, T ~ 2) %>% factor
# db4$ATORVA_D_yes = case_when(db4$DDD_ATORVASTATIN <= 30 ~ 0, db4$DDD_ATORVASTATIN <= 91 ~ 1, T ~ 2) %>% factor
# db4$ROSUVA_D_yes = case_when(db4$DDD_ROSUVASTATIN <= 0 ~ 0, db4$DDD_ROSUVASTATIN <= 182 ~ 1, T ~ 2) %>% factor
STATIN_sw = get_sw(exp_form_P1, data = db4, trunc_cut = 0.001)
ATORVA_sw = get_sw(exp_form_P2, data = db4, trunc_cut = 0.001)
ROSUVA_sw = get_sw(exp_form_P3, data = db4, trunc_cut = 0.001)

fit_cs1 = coxph(Surv(C_TIME1, C_EVENT1 == 1) ~ STATIN_D_yes, data = db4, weights = STATIN_sw$weight, robust = T)
fit_cs2 = coxph(Surv(C_TIME1, C_EVENT1 == 1) ~ ATORVA_D_yes + ageg4 + chf_yes + pvd_yes + cvd_yes + 
                  ld_yes + db_yes + ARB_P_yes + ASPIRIN_P_yes + METFORMIN_P_yes, data = db4, weights = ATORVA_sw$weight, robust = T)
fit_cs3 = coxph(Surv(C_TIME1, C_EVENT1 == 1) ~ ROSUVA_D_yes, data = db4, weights = ROSUVA_sw$weight, robust = T)

fg_dat1 = crprep("C_TIME1", "C_EVENT1", data = db4, trans = 1, cens = 0, id = "JID", 
                 keep = c("STATIN_D_yes", "ATORVA_D_yes", "ROSUVA_D_yes", x_var3)) %>%
  left_join(data.frame(JID = db4$JID, STATIN_sw = STATIN_sw$weight, ATORVA_sw = ATORVA_sw$weight, ROSUVA_sw = ROSUVA_sw$weight), by = "JID") %>%
  mutate(STATIN_sw = STATIN_sw * weight.cens,
         ATORVA_sw = ATORVA_sw * weight.cens,
         ROSUVA_sw = ROSUVA_sw * weight.cens)

fit_fg1 = coxph(Surv(Tstart, Tstop, status == 1) ~ STATIN_D_yes, data = fg_dat1, weights = STATIN_sw, robust = T)
fit_fg2 = coxph(Surv(Tstart, Tstop, status == 1) ~ ATORVA_D_yes + ageg4 + chf_yes + pvd_yes + cvd_yes + 
                  ld_yes + db_yes + ARB_P_yes + ASPIRIN_P_yes + METFORMIN_P_yes, data = fg_dat1, weights = ATORVA_sw, robust = T)
fit_fg3 = coxph(Surv(Tstart, Tstop, status == 1) ~ ROSUVA_D_yes, data = fg_dat1, weights = ROSUVA_sw, robust = T)

fit_res_csfg = 
  cbind.data.frame(
    model = rep(c("Cause-specific", "Fine-Gray"), each = 3),
    exposure = rep(c("STATIN", "ATORVASTATIN", "ROSUVASTATIN"), 2),
    rbind(
      apply(summary(fit_cs1)$conf.int[1:2,c(1,3,4)], 1, function(x) sprintf("%.2f (%.2f, %.2f)", x[1], x[2], x[3])),
      apply(summary(fit_cs2)$conf.int[1:2,c(1,3,4)], 1, function(x) sprintf("%.2f (%.2f, %.2f)", x[1], x[2], x[3])),
      apply(summary(fit_cs3)$conf.int[1:2,c(1,3,4)], 1, function(x) sprintf("%.2f (%.2f, %.2f)", x[1], x[2], x[3])),
      apply(summary(fit_fg1)$conf.int[1:2,c(1,3,4)], 1, function(x) sprintf("%.2f (%.2f, %.2f)", x[1], x[2], x[3])),
      apply(summary(fit_fg2)$conf.int[1:2,c(1,3,4)], 1, function(x) sprintf("%.2f (%.2f, %.2f)", x[1], x[2], x[3])),
      apply(summary(fit_fg3)$conf.int[1:2,c(1,3,4)], 1, function(x) sprintf("%.2f (%.2f, %.2f)", x[1], x[2], x[3]))
    )
  )
colnames(fit_res_csfg)[3:4] = c("Low vs Moderate", "Low vs High")
write.csv(fit_res_csfg, "spj_MGUS_TableS6.csv", row.names = F)

exp_form_P3 = paste0("ROSUVA_D_yes ~ ", paste0(x_var4, collapse = "+"), "+ageg4*ARB_P_yes+ATORVA_D_yes") %>% as.formula
love_P3 = get_loveplot(exp_form_P3, data = db4_t2, var_lab = c(my_var_lab, "Atorva"), title = "Rosuvastatin", col_lab = DDD_love_lab, trunc_cut = 0.001)
round(ExtractSmd(love_P3$adj_dat), 3) %>% max

ROSUVA_sw = get_sw(exp_form_P3, data = db4, trunc_cut = 0.001)
fit_cs3 = coxph(Surv(C_TIME1, C_EVENT1 == 1) ~ ROSUVA_D_yes+SEX_TP_CD+chf_yes+pud_yes+ATORVA_D_yes, data = db4, weights = ROSUVA_sw$weight, robust = T)
summary(fit_cs3)$conf.int

#### Figure 1. Weighted cumulative incidence ####
fit_1A = survfit(Surv(Tstart, Tstop, status == 1) ~ STATIN_D_yes, data = fg_dat1, weights = STATIN_sw)
fit_1B = survfit(Surv(Tstart, Tstop, status == 1) ~ ATORVA_D_yes, data = fg_dat1, weights = ATORVA_sw)
fit_1C = survfit(Surv(Tstart, Tstop, status == 1) ~ ROSUVA_D_yes, data = fg_dat1, weights = ROSUVA_sw)

gray_1A = coxph(Surv(Tstart, Tstop, status == 1) ~ STATIN_D_yes, data = fg_dat1, weights = STATIN_sw, robust = T)
gray_1B = coxph(Surv(Tstart, Tstop, status == 1) ~ ATORVA_D_yes, data = fg_dat1, weights = ATORVA_sw, robust = T)
gray_1C = coxph(Surv(Tstart, Tstop, status == 1) ~ ROSUVA_D_yes, data = fg_dat1, weights = ROSUVA_sw, robust = T)

# Gray t-test
summary(gray_1A)$robscore
summary(gray_1B)$robscore
summary(gray_1C)$robscore

fit_1A$n.risk[fit_1A$time == 0] = as.numeric(table(db4$STATIN_D_yes))
fit_1B$n.risk[fit_1B$time == 0] = as.numeric(table(db4$ATORVA_D_yes))
fit_1C$n.risk[fit_1C$time == 0] = as.numeric(table(db4$ROSUVA_D_yes))

plot_1A = 
  fit_1A %>%
  ggsurvplot(risk.table = T, legend = "none", conf.int = F, censor = F, xlab = "Time (Month)", ylab = "Survival probability", 
             legend.labs = c("Low", "Moderate", "High"),
             xlim = c(0, 72), break.x.by = 12, risk.table.height = 0.3, fontsize = 4,
             legend.title = "", palette = DDD_col, fun = 'event', title = "(A)")

plot_1B = 
  fit_1B %>%
  ggsurvplot(risk.table = T, legend = "none", conf.int = F, censor = F, xlab = "Time (Month)", ylab = "Survival probability", 
             legend.labs = c("Low", "Moderate", "High"),
             xlim = c(0, 72), break.x.by = 12, risk.table.height = 0.3, fontsize = 4,
             legend.title = "", palette = DDD_col, fun = 'event', title = "(B)")

plot_1C = 
  fit_1C %>%
  ggsurvplot(risk.table = T, legend = "none", conf.int = F, censor = F, xlab = "Time (Month)", ylab = "Survival probability", 
             legend.labs = c("Low", "Moderate", "High"),
             xlim = c(0, 72), break.x.by = 12, risk.table.height = 0.3, fontsize = 4,
             legend.title = "", palette = DDD_col, fun = 'event', title = "(C)")

plot_1A$plot = 
  plot_1A$plot + 
  scale_y_continuous(breaks = seq(0, 0.06, 0.015), labels = seq(0, 0.06, 0.015)*100, limits = c(0, 0.07)) +
  labs(y = "Cumulative incidence (%)") +
  theme(axis.title.y = element_text(vjust = -7))

plot_1B$plot = 
  plot_1B$plot + 
  scale_y_continuous(breaks = seq(0, 0.06, 0.015), labels = seq(0, 0.06, 0.015)*100, limits = c(0, 0.07)) +
  labs(y = "Cumulative incidence (%)") +
  theme(axis.title.y = element_text(vjust = -7))

plot_1C$plot = 
  plot_1C$plot + 
  scale_y_continuous(breaks = seq(0, 0.06, 0.015), labels = seq(0, 0.06, 0.015)*100, limits = c(0, 0.07)) +
  labs(y = "Cumulative incidence (%)") +
  theme(axis.title.y = element_text(vjust = -7))

fit_res = 
  cbind.data.frame(
    Level = c("Moderate", "High"),
    `HR (95% CI)` = 
      c(apply(summary(fit_cs1)$conf.int[1:2,c(1,3,4)], 1, function(x) sprintf("%.2f (%.2f, %.2f)", x[1], x[2], x[3])),
        apply(summary(fit_cs2)$conf.int[1:2,c(1,3,4)], 1, function(x) sprintf("%.2f (%.2f, %.2f)", x[1], x[2], x[3])),
        apply(summary(fit_cs3)$conf.int[1:2,c(1,3,4)], 1, function(x) sprintf("%.2f (%.2f, %.2f)", x[1], x[2], x[3]))),
    `P-value` = 
      c(sprintf("%.3f", summary(fit_cs1)$coefficients[1:2,6]),
        sprintf("%.3f", summary(fit_cs2)$coefficients[1:2,6]),
        sprintf("%.3f", summary(fit_cs3)$coefficients[1:2,6]))
  )


plot_1A$plot = plot_1A$plot + annotate(geom = "table2", x = 18, y = 0.065, label = list(fit_res[1:2,])) +
  annotate(geom = "text", x = 60, y = 0.0075, label = "P = 0.335")
plot_1B$plot = plot_1B$plot + annotate(geom = "table2", x = 18, y = 0.065, label = list(fit_res[3:4,])) +
  annotate(geom = "text", x = 60, y = 0.0075, label = "P = 0.682")
plot_1C$plot = plot_1C$plot + annotate(geom = "table2", x = 18, y = 0.065, label = list(fit_res[5:6,])) +
  annotate(geom = "text", x = 60, y = 0.0075, label = "P = 0.002")

ggsave("spj_MGUS_Figure1A.png", plot_1A, device = "png", dpi = 600, width = 6, height = 5)
ggsave("spj_MGUS_Figure1B.png", plot_1B, device = "png", dpi = 600, width = 6, height = 5)
ggsave("spj_MGUS_Figure1C.png", plot_1C, device = "png", dpi = 600, width = 6, height = 5)

ggsave("spj_MGUS_Figure1A.pdf", plot_1A, device = "pdf", dpi = 600, width = 6, height = 5)
ggsave("spj_MGUS_Figure1B.pdf", plot_1B, device = "pdf", dpi = 600, width = 6, height = 5)
ggsave("spj_MGUS_Figure1C.pdf", plot_1C, device = "pdf", dpi = 600, width = 6, height = 5)


# -- sensitivity analysis -- #
exp_form_P1 = paste0("STATIN_D_yes ~ ", paste0(x_var4, collapse = "+")) %>% as.formula
exp_form_P2 = paste0("ATORVA_D_yes ~ ", paste0(x_var4, collapse = "+")) %>% as.formula
# exp_form_P3 = paste0("ROSUVA_D_yes ~ ", paste0(x_var4, collapse = "+")) %>% as.formula
exp_form_P3 = paste0("ROSUVA_D_yes ~ ", paste0(x_var4, collapse = "+"), "+ageg4*ARB_P_yes") %>% as.formula

my_var_lab = c("Sex", "Age", "MI", "CHF", "PVD", "CVD", "Dementia", "CPD", "RHE", "PUD",
               "LD", "Diabetes", "HP", "RD", "Cancer", "ARB", "Aspirin", "Metformin")

DDD_love_lab = c("Low vs Moderate", "Low vs High", "Moderate vs High")

sens_trunc = c(0, 0.001, 0.01, 0.1)
sens_weight_res = list()
sens_res = list()
for(i in 1:length(sens_trunc)){
  
  love_P1 = get_loveplot(exp_form_P1, data = db4_t2, var_lab = my_var_lab, title = "Statin", col_lab = DDD_love_lab, trunc_cut = sens_trunc[i])
  love_P2 = get_loveplot(exp_form_P2, data = db4_t2, var_lab = my_var_lab, title = "Atorvastatin", col_lab = DDD_love_lab, trunc_cut = sens_trunc[i])
  love_P3 = get_loveplot(exp_form_P3, data = db4_t2, var_lab = my_var_lab, title = "Rosuvastatin", col_lab = DDD_love_lab, trunc_cut = sens_trunc[i])
  
  cov_P1 = rownames(ExtractSmd(love_P1$adj_dat))[apply(ExtractSmd(love_P1$adj_dat)[,2:4], 1, max) > 0.1]
  cov_P2 = rownames(ExtractSmd(love_P2$adj_dat))[apply(ExtractSmd(love_P2$adj_dat)[,2:4], 1, max) > 0.1]
  cov_P3 = rownames(ExtractSmd(love_P3$adj_dat))[apply(ExtractSmd(love_P3$adj_dat)[,2:4], 1, max) > 0.1]
  
  STATIN_sw = get_sw(exp_form_P1, data = db4, trunc_cut = sens_trunc[i])
  ATORVA_sw = get_sw(exp_form_P2, data = db4, trunc_cut = sens_trunc[i])
  ROSUVA_sw = get_sw(exp_form_P3, data = db4, trunc_cut = sens_trunc[i])
  
  if(length(cov_P1) > 0) {cov_PP1 = "+"} else {cov_PP1 = ""}
  if(length(cov_P2) > 0) {cov_PP2 = "+"} else {cov_PP2 = ""}
  if(length(cov_P3) > 0) {cov_PP3 = "+"} else {cov_PP3 = ""}
  
  new_form_P1 = paste0("Surv(C_TIME1, C_EVENT1 == 1) ~ STATIN_D_yes", cov_PP1, paste0(cov_P1, collapse = "+")) %>% as.formula
  new_form_P2 = paste0("Surv(C_TIME1, C_EVENT1 == 1) ~ ATORVA_D_yes", cov_PP2, paste0(cov_P2, collapse = "+")) %>% as.formula
  new_form_P3 = paste0("Surv(C_TIME1, C_EVENT1 == 1) ~ ROSUVA_D_yes", cov_PP3, paste0(cov_P3, collapse = "+")) %>% as.formula
  
  fit_cs1 = coxph(new_form_P1, data = db4, weights = STATIN_sw$weight, robust = T)
  fit_cs2 = coxph(new_form_P2, data = db4, weights = ATORVA_sw$weight, robust = T)
  fit_cs3 = coxph(new_form_P3, data = db4, weights = ROSUVA_sw$weight, robust = T)
  
  sens_weight_res[[i]] = data.frame(Treatment = c("Statin", "Atorva", "Rosuva"),
                                    Mean = c(mean(STATIN_sw$weight), mean(ATORVA_sw$weight), mean(ROSUVA_sw$weight)),
                                    SD = c(sd(STATIN_sw$weight), sd(ATORVA_sw$weight), sd(ROSUVA_sw$weight)),
                                    Min = c(min(STATIN_sw$weight), min(ATORVA_sw$weight), min(ROSUVA_sw$weight)),
                                    Max = c(max(STATIN_sw$weight), max(ATORVA_sw$weight), max(ROSUVA_sw$weight)))
  
  sens_res[[i]] = 
    cbind.data.frame(
      Level = c("Moderate", "High"),
      `HR (95% CI)` = 
        c(apply(summary(fit_cs1)$conf.int[1:2,c(1,3,4)], 1, function(x) sprintf("%.2f (%.2f, %.2f)", x[1], x[2], x[3])),
          apply(summary(fit_cs2)$conf.int[1:2,c(1,3,4)], 1, function(x) sprintf("%.2f (%.2f, %.2f)", x[1], x[2], x[3])),
          apply(summary(fit_cs3)$conf.int[1:2,c(1,3,4)], 1, function(x) sprintf("%.2f (%.2f, %.2f)", x[1], x[2], x[3]))),
      `P-value` = 
        c(sprintf("%.3f", summary(fit_cs1)$coefficients[1:2,6]),
          sprintf("%.3f", summary(fit_cs2)$coefficients[1:2,6]),
          sprintf("%.3f", summary(fit_cs3)$coefficients[1:2,6]))
    )
  
}

bind_cols(
  do.call("bind_rows", 
          lapply(sens_weight_res, function(x){
            x$res1 = sprintf("%.2f (%.2f)", x$Mean, x$SD)
            x$res2 = sprintf("%.2f / %.2f", x$Min, x$Max)
            x %>% select(Treatment, res1, res2)
          })
  ),
  do.call("bind_rows", 
          lapply(sens_res, function(x){
            x$trt = factor(rep(c("Statin", "Atorva", "Rosuva"), each = 2), levels = c("Statin", "Atorva", "Rosuva"))
            x$Level = factor(x$Level, levels = c("Moderate", "High"))
            x %>% select(-`P-value`) %>% spread(key = "Level", value = "HR (95% CI)") %>%
              mutate(trt = as.character(trt), Moderate = as.character(Moderate), High = as.character(High))
          })
  )
) %>%
  select(-trt) %>%
  write.csv("spj_MGUS_TableS7.csv", row.names = F)



