# -- library -- #
library(tidyverse)
library(survival)
library(nnet)
library(cmprsk)
library(EValue)

# -- function -- #
source('HIRA_statin_function.R')

# -- data -- #
source('HIRA_statin_cohort.R')

LV <- c("Low (<=30)", "Moderate (31-365)", "High (>=366)")
for (v in c("statin_grp", "atorva_grp", "rosuva_grp")){
  dat[[v]] <- factor(dat[[v]], levels = LV)
}
  
dat$sex <- factor(dat$sex, levels = c("Male", "Female"))
dat$age_group <- factor(dat$age_group, levels = c("<60", "60-69", ">=70"))
dat$male <- as.integer(dat$sex == "Male")
dat$ev <- factor(dat$status, levels = 0:2, labels = c("censor", "MM", "death"))

MAXFU <- 6
dat$status = case_when(dat$time_years >= MAXFU ~ 0, T ~ dat$status)
dat$time_years <- pmin(dat$time_years, MAXFU)

COVS_DEMO <- c("age", "male")
COVS_COMO <- c("mi", "chf", "pvd", "cvd", "dementia", "cpd", "rheum",
               "pud", "liver", "dm", "hemiplegia", "renal", "cancer")
COVS_MED  <- c("acei_arb", "aspirin", "ezetimibe", "metformin")
COVS_ALL  <- c(COVS_DEMO, COVS_COMO, COVS_MED)   # IPTW / 다변량 보정 공변량

# -- analysis -- #
rev_km <- survfit(Surv(time_years, as.integer(status == 0)) ~ 1, data = dat)
mfu <- quantile(rev_km, 0.5)
cat(sprintf("Median follow-up (reverse KM): %.1f months (95%% CI, %.1f-%.1f)\n", mfu$quantile * 12, mfu$lower * 12, mfu$upper * 12))
cat(sprintf("Median observed follow-up : %.1f months\n", median(dat$time_years) * 12))

ev_tab <- data.frame(item = c("Censoring", "Progression to MM", "Death"),
                     n = as.numeric(table(dat$ev)),
                     pct = round(100 * as.numeric(prop.table(table(dat$ev))), 1))
print(ev_tab)

# baseline characteristics
desc_overall <- rbind(
  data.frame(variable = "Age, mean (SD)",
             value = sprintf("%.1f (%.1f)", mean(dat$age), sd(dat$age))),
  data.frame(variable = paste0("Age ", levels(dat$age_group)),
             value = sprintf("%d (%.1f%%)", as.numeric(table(dat$age_group)),
                             100 * as.numeric(prop.table(table(dat$age_group))))),
  data.frame(variable = c("Male", "Female"),
             value = sprintf("%d (%.1f%%)", as.numeric(table(dat$sex)),
                             100 * as.numeric(prop.table(table(dat$sex))))),
  data.frame(variable = c(COVS_COMO, "statin", COVS_MED),
             value = sapply(c(COVS_COMO, "statin", COVS_MED), function(v)
               sprintf("%d (%.1f%%)", sum(dat[[v]]), 100 * mean(dat[[v]]))))
)
print(desc_overall, row.names = FALSE)

# cumulative incidence
aj <- survfit(Surv(time_years, ev) ~ 1, data = dat)
s_aj <- summary(aj, times = c(1, 2, 4, 6))
pm <- s_aj$pstate; colnames(pm) <- aj$states
cif_overall <- data.frame(
  years = s_aj$time,
  CIF_MM_pct    = round(100 * pm[, "MM"], 2),
  CIF_death_pct = round(100 * pm[, "death"], 2))
cat("\n-- Cumulative incidence (%) --\n"); print(cif_overall, row.names = FALSE)

km_prog <- survfit(Surv(time_years, prog) ~ 1, data = subset(dat, prog == 1))
q <- quantile(km_prog, 0.5)
cat(sprintf("\nMedian time to progression (progressors only): %.1f months (95%% CI, %.1f-%.1f)\n", q$quantile * 12, q$lower * 12, q$upper * 12))

# baseline characteristics - statin
t1 <- table1(dat, "statin_grp",
             cont_vars = c("age", "atorva_cddd", "rosuva_cddd"),
             bin_vars  = c("male", "prog", "death", COVS_COMO, COVS_MED))
names(t1) <- c("variable", paste0(LV, " (n=", as.numeric(table(dat$statin_grp)), ")"), "P value")
print(t1, row.names = FALSE)

# screening
DRUGS <- c("statin", "acei_arb", "aspirin", "ezetimibe", "metformin")

etab7 <- do.call(rbind, lapply(DRUGS, screen_one, data = dat))
print(etab7, row.names = FALSE)

# sensitivity
etab8 <- list()
for (w in 1:3) {
  for (thr in c(0.70, 0.80, 0.90)) {
    d <- dat
    cutoff <- 365 * thr
    for (dr in DRUGS) d[[dr]] <- as.integer(d[[paste0(dr, "_cddd_w", w)]] >= cutoff)
    res <- do.call(rbind, lapply(DRUGS, function(dr) {
      r <- screen_one(dr, d); data.frame(drug = dr, HR = r$multivariate) }))
    etab8[[length(etab8) + 1]] <- data.frame(
      window = paste0(w, " year"), threshold = paste0(thr * 100, "%"),
      as.data.frame(t(setNames(res$HR, res$drug))), check.names = FALSE)
  }
}
etab8 <- do.call(rbind, etab8)
print(etab8, row.names = FALSE)

# IPT weighting
iptw <- lapply(c(statin = "statin_grp", atorva = "atorva_grp", rosuva = "rosuva_grp"),
               function(v) make_iptw(dat, v, COVS_ALL))
dat$w_statin <- iptw$statin$w
dat$w_atorva <- iptw$atorva$w
dat$w_rosuva <- iptw$rosuva$w

cat("Stabilized weight summary\n")
print(round(sapply(dat[, c("w_statin", "w_atorva", "w_rosuva")], summary), 3))

# SMD (보정 전 / 후) 및 double adjustment 대상 선정
bal <- lapply(c(statin = "statin_grp", atorva = "atorva_grp", rosuva = "rosuva_grp"),
              function(v) {
                b0 <- smd_table(dat, v, COVS_ALL)
                b1 <- smd_table(dat, v, COVS_ALL, w = dat[[paste0("w_", sub("_grp", "", v))]])
                data.frame(variable = COVS_ALL,
                           SMD_unweighted = b0$max_abs_SMD, SMD_weighted = b1$max_abs_SMD)
              })
for (nm in names(bal)) {
  cat("\n--", nm, ": SMD >0.1 after weighting (double adjustment 대상) --\n")
  print(subset(bal[[nm]], SMD_weighted > 0.1), row.names = FALSE)
}
dbl_adj <- lapply(bal, function(b) b$variable[b$SMD_weighted > 0.1])

# love plot
png(file.path(OUT, "eFigure3_love_plot.png"), width = 1500, height = 600, res = 130)
par(mfrow = c(1, 3), mar = c(4, 8, 3, 1))
for (nm in names(bal)) {
  b <- bal[[nm]]; o <- order(b$SMD_unweighted)
  plot(b$SMD_unweighted[o], seq_along(o), pch = 1, col = "grey40",
       xlim = c(0, max(0.3, max(b$SMD_unweighted))), yaxt = "n", ylab = "",
       xlab = "Standardized mean difference", main = nm)
  points(b$SMD_weighted[o], seq_along(o), pch = 16, col = "firebrick")
  axis(2, at = seq_along(o), labels = b$variable[o], las = 1, cex.axis = 0.7)
  abline(v = 0.1, lty = 2)
  legend("bottomright", c("Unweighted", "IPTW"), pch = c(1, 16),
         col = c("grey40", "firebrick"), bty = "n", cex = 0.8)
}
dev.off()

# incidence rate
inc_tab <- do.call(rbind, lapply(
  c("Any statin" = "statin_grp", "Atorvastatin" = "atorva_grp",
    "Rosuvastatin" = "rosuva_grp"),
  function(v) {
    do.call(rbind, lapply(LV, function(l) {
      s <- dat[dat[[v]] == l, ]
      py <- sum(s$time_years); e <- sum(s$prog)
      r <- inc_rate(e, py)
      data.frame(drug = v, group = l, patients = nrow(s),
                 events = e, pct = round(100 * e / nrow(s), 1),
                 person_years = round(py, 1),
                 incidence = sprintf("%.1f (%.1f-%.1f)", r[1], r[2], r[3]))
    }))
  }))
inc_tab$drug <- rep(c("Any statin", "Atorvastatin", "Rosuvastatin"), each = 3)
print(inc_tab, row.names = FALSE)

# IPTW analysis
par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3, 1))
gray_res <- list()
for (k in seq_along(c("statin", "atorva", "rosuva"))) {
  nm  <- c("statin", "atorva", "rosuva")[k]
  ttl <- c("A. Any statin", "B. Atorvastatin", "C. Rosuvastatin")[k]
  v   <- paste0(nm, "_grp"); wv <- paste0("w_", nm)
  
  fit <- survfit(as.formula(paste("Surv(time_years, ev) ~", v)),
                 data = dat, weights = dat[[wv]])
  plot(fit[, "MM"], col = c("grey30", "steelblue", "firebrick"), lwd = 2,
       xlab = "Years since landmark", ylab = "Cumulative incidence of MM",
       ylim = c(0, 0.10), main = ttl)
  legend("topleft", LV, col = c("grey30", "steelblue", "firebrick"),
         lwd = 2, bty = "n", cex = 0.8)
  
  gt <- cmprsk::cuminc(ftime = dat$time_years, fstatus = dat$status,
                       group = dat[[v]], cencode = 0)
  p <- gt$Tests["1", "pv"]
  gray_res[[nm]] <- data.frame(drug = nm, gray_p = signif(p, 3))
  mtext(sprintf("Gray test P = %.3f", p), side = 3, line = -1.5, cex = 0.7)
}

wcox <- function(nm, extra = NULL, data = dat) {
  v <- paste0(nm, "_grp"); wv <- paste0("w_", nm)
  rhs <- paste(c(v, dbl_adj[[nm]], extra), collapse = " + ")
  fit <- coxph(as.formula(paste("Surv(time_years, prog) ~", rhs)),
               data = data, weights = data[[wv]], id = id, robust = TRUE)
  r <- hr_ci(fit, paste0(v, LV[-1]))
  r$term <- LV[-1]; r$drug <- nm
  r[, c("drug", "term", "HR", "lower", "upper", "label")]
}
cox_main <- do.call(rbind, lapply(c("statin", "atorva", "rosuva"), wcox))
print(cox_main, row.names = FALSE)
sav(cox_main, "Table2_weighted_HR_cause_specific.csv")

# atorva adjust
cox_ros_adj <- wcox("rosuva", extra = "atorva_grp")
print(cox_ros_adj, row.names = FALSE)

# subdistribution
fg_model <- function(nm, data = dat) {
  v <- paste0(nm, "_grp"); wv <- paste0("w_", nm)
  d <- data; d$w <- d[[wv]]
  fgd <- finegray(as.formula(paste("Surv(time_years, ev) ~", v, "+ w + id")),
                  data = d, etype = "MM", weights = w)
  fit <- coxph(as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~", v)),
               data = fgd, weights = fgwt, id = id, robust = TRUE)
  r <- hr_ci(fit, paste0(v, LV[-1])); r$term <- LV[-1]; r$drug <- nm
  r[, c("drug", "term", "HR", "lower", "upper", "label")]
}
fg_res <- do.call(rbind, lapply(c("statin", "atorva", "rosuva"), fg_model))
print(fg_res, row.names = FALSE)

# weight truncation sensitivity
trunc_sens <- do.call(rbind, lapply(
  list(list("No truncation", NULL), list("99th %ile", 0.99),
       list("99.5th %ile", 0.995), list("99.9th %ile", 0.999),
       list("99.95th %ile (main)", 0.9995)),
  function(s) {
    ww <- make_iptw(dat, "rosuva_grp", COVS_ALL, trunc_q = s[[2]])$w
    d <- dat; d$w_rosuva <- ww
    r <- wcox("rosuva", data = d)
    data.frame(method = s[[1]], max_weight = round(max(ww), 2),
               moderate = r$label[1], high = r$label[2])
  }))
print(trunc_sens, row.names = FALSE)

# E-value
ev_tab2 <- do.call(rbind, lapply(seq_len(nrow(cox_main)), function(i) {
  r <- cox_main[i, ]
  data.frame(drug = r$drug, group = r$term, HR = r$HR,
             E_value_point = round(evalue(r$HR), 2),
             E_value_CI    = round(evalue_ci(r$HR, r$lower, r$upper), 2))
}))
print(ev_tab2, row.names = FALSE)
