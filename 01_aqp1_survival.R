# =============================================================
# 01_aqp1_survival.R
# AQP1 expression and overall survival in TCGA-KIRC
# (clear cell renal cell carcinoma)
#
# Data: UCSC Xena, GDC TCGA-KIRC
#   - TCGA-KIRC.star_tpm.tsv.gz   gene expression, log2(TPM+1)
#   - TCGA-KIRC.survival.tsv.gz   overall survival
#   - TCGA-KIRC.clinical.tsv.gz   clinical phenotype
# =============================================================


# ---- 0. Packages --------------------------------------------
library(tidyverse)
library(survival)
library(survminer)


# ---- 1. Read data -------------------------------------------
expr <- read_tsv("TCGA-KIRC.star_tpm.tsv.gz")
surv <- read_tsv("TCGA-KIRC.survival.tsv.gz")
clin <- read_tsv("TCGA-KIRC.clinical.tsv.gz")


# ---- 2. Extract AQP1 (Ensembl ID ENSG00000240583) -----------
# One row (AQP1) x 610 sample columns -> one row per sample
aqp1 <- expr %>%
  filter(str_starts(Ensembl_ID, "ENSG00000240583")) %>%
  pivot_longer(-Ensembl_ID, names_to = "sample", values_to = "AQP1") %>%
  select(sample, AQP1)


# ---- 3. Keep primary tumor samples only ---------------------
# Characters 14-15 of a TCGA barcode give the sample type:
# 01 = primary tumor, 11 = solid tissue normal, 05 = additional new primary
table(substr(aqp1$sample, 14, 15))

aqp1_tumor <- aqp1 %>%
  filter(substr(sample, 14, 15) == "01")


# ---- 4. Merge with survival data ----------------------------
df <- inner_join(aqp1_tumor, surv, by = "sample")


# ---- 5. One sample per patient ------------------------------
# 4 patients have two tumor samples (01A and 01B); keep the first (01A)
df_clean <- df %>%
  distinct(`_PATIENT`, .keep_all = TRUE)

nrow(df_clean)          # 529 patients
table(df_clean$OS)      # 0 = censored, 1 = death


# ---- 6. Group patients by median AQP1 -----------------------
# Low is the reference group for Cox regression
df_clean <- df_clean %>%
  mutate(AQP1_group = ifelse(AQP1 > median(AQP1), "High", "Low"),
         AQP1_group = factor(AQP1_group, levels = c("Low", "High")),
         OS.year    = OS.time / 365)


# ---- 7. Kaplan-Meier curve + log-rank test ------------------
fit <- survfit(Surv(OS.year, OS) ~ AQP1_group, data = df_clean)
fit

p_km <- ggsurvplot(fit,
                   data = df_clean,
                   pval = TRUE,
                   conf.int = TRUE,
                   risk.table = TRUE,
                   legend.labs = c("AQP1 Low", "AQP1 High"),
                   legend.title = "",
                   xlab = "Time (years)",
                   palette = c("#2E86C1", "#E67E22"))

dir.create("figures", showWarnings = FALSE)

png("figures/aqp1_km.png", width = 8, height = 7, units = "in", res = 300)
print(p_km)
dev.off()


# ---- 8. Univariate Cox regression ---------------------------
cox_uni <- coxph(Surv(OS.year, OS) ~ AQP1_group, data = df_clean)
summary(cox_uni)        # HR = 0.36 (0.26-0.50)


# ---- 9. Add age and pathologic stage ------------------------
clin_sub <- clin %>%
  select(sample,
         age   = age_at_index.demographic,
         stage = ajcc_pathologic_stage.diagnoses,
         gender = gender.demographic)

df_multi <- df_clean %>%
  left_join(clin_sub, by = "sample") %>%
  mutate(stage = factor(stage,
                        levels = c("Stage I", "Stage II", "Stage III", "Stage IV")),
         gender = factor(gender, levels = c("female", "male"))) %>%
  filter(!is.na(age), !is.na(stage))

nrow(df_multi)          # 526 patients (3 missing stage)


# ---- 10. Multivariable Cox regression + forest plot ---------
cox_multi <- coxph(Surv(OS.year, OS) ~ AQP1_group + age + stage,
                   data = df_multi)
# Sensitivity analysis: also adjust for sex
cox_multi2 <- coxph(Surv(OS.year, OS) ~ AQP1_group + age + stage + gender,
                    data = df_multi)
summary(cox_multi2)     # AQP1 High: HR = 0.51 (0.36-0.72)

summary(cox_multi)      # AQP1 High: HR = 0.55 (0.39-0.77)

# ggforest() needs a plain data.frame, not a tibble
p_forest <- ggforest(cox_multi, data = as.data.frame(df_multi))

png("figures/forest_plot.png", width = 10, height = 5.5, units = "in", res = 300)
print(p_forest)
dev.off()


# ---- 11. AQP1 expression by stage ---------------------------
df_multi %>%
  group_by(stage) %>%
  summarise(n = n(), median_AQP1 = median(AQP1))

p_stage <- ggplot(df_multi, aes(x = stage, y = AQP1, fill = stage)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 0.8, alpha = 0.4) +
  stat_compare_means(method = "kruskal.test") +
  labs(x = "Pathologic stage", y = "AQP1 expression (log2 TPM+1)") +
  theme_classic() +
  theme(legend.position = "none")

png("figures/aqp1_by_stage.png", width = 7, height = 5, units = "in", res = 300)
print(p_stage)
dev.off()

