source("code/levellogger_functions.R")
source("code/levellogger_packages.R")
source("code/levellogger_plan.R")

data <- readd(daily_water_balance)

data[, pet_cm_d := -pet_cm_d]

raw_mod <- 
  quantreg::rq(raw_et_cm_d ~ pet_cm_d, 
               data = data[(dry_day)])

corrected_mod <- 
  quantreg::rq(corrected_et_cm_d ~ pet_cm_d, 
               data = data[(dry_day)])

raw_t <- 
  t.test(data[(dry_day)]$raw_et_cm_d - coef(raw_mod)[[1]], data[(dry_day)]$pet_cm_d, paired = TRUE)
corr_t <- 
  t.test(data[(dry_day)]$corrected_et_cm_d - coef(corrected_mod)[[1]], data[(dry_day)]$pet_cm_d, paired = TRUE)

p.adjust(c(raw_t$p.value, corr_t$p.value), "bonferroni")
