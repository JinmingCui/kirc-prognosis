# =============================================================
# 02_nomogram.R
# Nomogram for overall survival in TCGA-KIRC
# (AQP1 group + age + pathologic stage)
#
# Run 01_aqp1_survival.R first: this script uses df_multi and
# cox_multi created there.
# =============================================================


# ---- 0. Packages --------------------------------------------
# install.packages("rms")   # run once
library(rms)


# ---- 1. Fit the Cox model with rms --------------------------
# rms needs a plain data.frame and a datadist object
nom_data <- as.data.frame(df_multi)
dd <- datadist(nom_data)
options(datadist = "dd")

# Sex was not included: it was not significant and did not
# improve the C-index (see 01_aqp1_survival.R)
fit_cph <- cph(Surv(OS.year, OS) ~ AQP1_group + age + stage,
               data = nom_data, x = TRUE, y = TRUE, surv = TRUE)
fit_cph


# ---- 2. Nomogram --------------------------------------------
surv_fun <- Survival(fit_cph)

nom <- nomogram(fit_cph,
                fun = list(function(x) surv_fun(1, x),
                           function(x) surv_fun(3, x),
                           function(x) surv_fun(5, x)),
                funlabel = c("1-year survival", "3-year survival", "5-year survival"),
                lp = FALSE)

png("figures/nomogram.png", width = 11, height = 6, units = "in", res = 300)
plot(nom)
dev.off()


# ---- 3. Predicted survival for an example patient -----------
# 60 years old, Stage III, AQP1 low
patient_A <- data.frame(
  AQP1_group = factor("Low", levels = c("Low", "High")),
  age = 60,
  stage = factor("Stage III", levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
)

lp_A <- predict(fit_cph, newdata = patient_A)
surv_fun(c(1, 3, 5), lp_A)   # 1-, 3-, 5-year survival: 0.87, 0.68, 0.52


# ---- 4. Calibration at 5 years (bootstrap, B = 200) ---------
fit_cph5 <- cph(Surv(OS.year, OS) ~ AQP1_group + age + stage,
                data = nom_data, x = TRUE, y = TRUE, surv = TRUE,
                time.inc = 5)

set.seed(2026)
cal5 <- calibrate(fit_cph5, cmethod = "KM", method = "boot",
                  u = 5, m = 105, B = 200)

png("figures/calibration_5year.png", width = 7, height = 6, units = "in", res = 300)
plot(cal5,
     xlab = "Predicted 5-year survival",
     ylab = "Observed 5-year survival")
dev.off()


# ---- 5. Does AQP1 add to age + stage? -----------------------
cox_clin <- coxph(Surv(OS.year, OS) ~ age + stage, data = df_multi)

summary(cox_clin)$concordance    # C-index without AQP1: 0.761
summary(cox_multi)$concordance   # C-index with AQP1:    0.770
anova(cox_clin, cox_multi)       # likelihood ratio test: p = 0.0004
