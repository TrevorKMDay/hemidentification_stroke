library(tidyverse)

# partial r-squared
library(rsq)

# Setup

setwd("~/MyDrive/Projects/hemisphere_fingerprinting_stroke/code/")
source("../../hemisphere_fingerprinting/code/bootstrap_mcc.R")
source("acc_se.R")

centroids <- read_csv("datasets/centroids_all.csv", show_col_types = FALSE)
ld1ld2_0 <- read_csv("datasets/ld1ld2_results.csv", show_col_types = FALSE)

# Count missing peaks

activ <- ld1ld2_0 %>%
  select(sub, PatientStatus, ends_with("_activ")) %>%
  distinct() %>%
  mutate(
    missing_activ =  is.na(lh_lang_activ)
  )

table(activ$PatientStatus)
table(activ$PatientStatus, activ$missing_activ)

# Ld correlation

lds_only <- ld1ld2_0 %>%
  select(sub, hemi, PatientStatus, LD1, LD2)

rcorr(lds_only$LD1, lds_only$LD2)
rcorr(abs(lds_only$LD1), lds_only$LD2)

lds_by_group <- lds_only %>%
  group_by(PatientStatus, hemi) %>%
  nest() %>%
  mutate(
    rcorr = map(data, ~rcorr(abs(.x$LD1), .x$LD2))
  )

hc_lh_centroid <- centroids %>%
  filter(
    class == "HC_LH",
    status == "train_and_test"
  ) %>%
  pivot_longer(ends_with("_m")) %>%
  pull(value)

hc_rh_centroid <- centroids %>%
  filter(
    class == "HC_RH",
    status == "train_and_test"
  ) %>%
  pivot_longer(ends_with("_m")) %>%
  pull(value)

# Distance calculation

distance <- function(a, b) {

  euclid_dist <- sqrt((a[1] - b[1])^2 + (a[2] - b[2])^2)
  return(euclid_dist)

}

# Add distance to data

ld1ld2_all <- ld1ld2_0 %>%
  mutate(

    # Distance from HC LH centroid to each point - might as well calculate it
    # for everyone
    lhc_dist = map2_dbl(LD1, LD2, ~distance(hc_lh_centroid, c(.x, .y))),

    rhc_dist = map2_dbl(LD1, LD2, ~distance(hc_rh_centroid, c(.x, .y))),

  ) 

ld1ld2 <- ld1ld2_all %>%
  filter_out(
    status %in% c("test_long", "test_rhlesion")
  )

ggplot(ld1ld2, aes(x = lhc_dist, color = class)) +
  geom_density() +
  theme_bw() +
  labs(title = "Euclid. dist. from HC LH centroid to hemi")

ggplot(ld1ld2, aes(x = rhc_dist, color = class)) +
  geom_density() +
  theme_bw() +
  labs(title = "Euclid. dist. from HC LH centroid to hemi")

ggplot(ld1ld2, aes(x = lhc_dist, y = rhc_dist, color = class)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  facet_wrap(vars(class)) +
  theme_bw()

ld1ld2_rhlesion <- ld1ld2_all %>%
  filter(
    status %in% c("train_and_test", "test_rhlesion")
  )

ld1ld2_rhonly <- ld1ld2_all %>%
  filter(
    hemi == "RH",
    status %in% c("train_and_test", "test_rhlesion")
  )

ld1ld2 %>%
  select(ends_with("activ")) %>%
  corrr::correlate()

# Basic modeling

table(ld1ld2_rhlesion$status)

lm_rhlesion_ld1 <- lm(LD1 ~ EHI_value + status, data = ld1ld2_rhlesion)
summary(lm_rhlesion_ld1)

lm_rhlesion_ld1 <- lm(LD2 ~ EHI_value + status, data = ld1ld2_rhlesion)
summary(lm_rhlesion_ld1)

# Which is better predicted by handedness?

lm_ld1_ehi <- lm(LD1 ~ PatientStatus + EHI_value, data = ld1ld2_rhonly)
lm_dist_ehi <- lm(lhc_dist ~ PatientStatus + EHI_value, data = ld1ld2_rhonly)

summary(lm_ld1_ehi)
summary(lm_dist_ehi)

AIC(lm_ld1_ehi, lm_dist_ehi)
map_dbl(list(lm_ld1_ehi, lm_dist_ehi), ~summary(.x)$r.squared)

rsq.partial(lm_ld1_ehi)
rsq.partial(lm_dist_ehi)

lm_ld2_lv <- lm(LD2 ~ logLV, data = ld1ld2_rhonly)
lm_dist_lv <- lm(lhc_dist ~ logLV, data = ld1ld2_rhonly)
lm_dist2_lv <- lm(rhc_dist ~ logLV, data = ld1ld2_rhonly)

summary(lm_ld2_lv)
summary(lm_dist_lv)
summary(lm_dist2_lv)

rsq.partial(lm_ld2_lv)
rsq.partial(lm_dist_lv)
rsq.partial(lm_dist2_lv)

AIC(lm_ld2_lv, lm_dist_lv, lm_dist2_lv)
map_dbl(list(lm_ld2_lv, lm_dist_lv, lm_dist2_lv), ~summary(.x)$r.squared)

# Tests ====

library(corrr)

basic_corr <- ld1ld2 %>%
  select(age_at_stroke, mri_age, mri_chronicity_mo, EHI_value, WAB_value, 
          logLV, pct_lang_mean, residual_pct_lang, 
          matches("[lr]h_.*_activ")) %>%
  correlate() %>%
  stretch()

ggplot(basic_corr, aes(x = x, y = y)) +
  geom_tile(aes(fill = r)) +
  geom_text(aes(label = round(r, 2)), size = 3) +
  scale_fill_gradient2(limits = c(-0.5, 0.5), oob = scales::squish) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## EHI =====

lm_ld1_lh_ehi <- lm(LD1 ~ EHI_value, data = filter(ld1ld2, hemi == "LH"))

lm_ld1_rh_ehi <- lm(LD1 ~ PatientStatus*EHI_value, 
                    data = filter(ld1ld2, hemi == "RH"))

summary(lm_ld1_lh_ehi)
summary(lm_ld1_rh_ehi)

lm_ld1_lh_demo <- lm(LD1 ~ EHI_value + Sex + mri_age, 
                      data = filter(ld1ld2, hemi == "LH"))

summary(lm_ld1_lh_demo)

lm_ld1_rh_hcdemo <- lm(LD1 ~ EHI_value + Sex + mri_age, 
                       data = filter(ld1ld2, class == "HC_RH"))

summary(lm_ld1_rh_hcdemo)

lm_ld1_rh_lhsdemo <- lm(LD1 ~ EHI_value + Sex + logLV + mri_age + 
                          age_at_stroke + residual_pct_lang, 
                        data = filter(ld1ld2, class == "LHS_RH"))

summary(lm_ld1_rh_lhsdemo)

summary(step(lm_ld1_rh_lhsdemo, direction = "ba"))

## Activation
lm_ld1_lh_activ <- lm(LD1 ~ EHI_value + lh_lang_activ + lh_semexec_activ, 
                      data = filter(ld1ld2, hemi == "LH"))

summary(lm_ld1_lh_activ)

lm_ld1_lh_activ1 <- lm(LD1 ~ EHI_value + lh_lang_activ, 
                      data = filter(ld1ld2, hemi == "LH"))

summary(lm_ld1_lh_activ1)

lm_ld1_hcrh_activ <- lm(LD1 ~ EHI_value + rh_lang_activ + rh_semexec_activ, 
                        data = filter(ld1ld2, class == "HC_RH"))

summary(lm_ld1_hcrh_activ)

lm_ld1_hcrh_activ1 <- lm(LD1 ~ EHI_value + rh_lang_activ, 
                        data = filter(ld1ld2, class == "HC_RH"))

summary(lm_ld1_hcrh_activ1)

lm_ld1_lhsrh_activ <- lm(LD1 ~ EHI_value + rh_lang_activ + rh_semexec_activ, 
                        data = filter(ld1ld2, class == "LHS_RH"))

summary(lm_ld1_lhsrh_activ)

lm_ld1_lhsrh_activ1 <- lm(LD1 ~ EHI_value + rh_lang_activ, 
                          data = filter(ld1ld2, class == "LHS_RH"))

summary(lm_ld1_lhsrh_activ1)


## LD2 

ld1ld2_lhsonly <- filter(ld1ld2, PatientStatus == "LHS") %>%
  mutate(
    log_chronicity = log10(mri_chronicity_mo)
  )

lm_ld2_logLV <- lm(LD2 ~ EHI_value + logLV, data = ld1ld2_lhsonly)

summary(lm_ld2_logLV)

lm_ld2_stroke <- lm(LD2 ~ logLV + Sex + mri_age + age_at_stroke + 
                      residual_pct_lang, 
                    data = ld1ld2_lhsonly)

summary(lm_ld2_stroke)

summary(step(lm_ld2_stroke, direction = "backward"))

table(is.na(ld1ld2_lhsonly$rh_lang_activ), is.na(ld1ld2_lhsonly$rh_semexec_activ))

lm_ld2_stroke_activ <- lm(LD2 ~ logLV + rh_lang_activ + rh_semexec_activ, 
                          data = ld1ld2_lhsonly)

summary(lm_ld2_stroke_activ)


lm_ld2_stroke_activ1 <- lm(LD2 ~ logLV + rh_lang_activ,data = ld1ld2_lhsonly)

summary(lm_ld2_stroke_activ1)

# Behavioral prediction ====
 

# Larger lesions lead to lower WA values
wab_prediction <- lm(WAB_value ~ EHI_value + logLV + residual_pct_lang + 
                        mri_age + log_chronicity,
                      data = ld1ld2_lhsonly) 

summary(wab_prediction)

# lavaan 

library(lavaan)
library(lavaanPlot)

results_lhs_lavaan <- ld1ld2_lhsonly %>%
  mutate(
    # These are scored out of [-100, 100] and [0, 100], normalize to [x, 1.0]
    EHI_value = EHI_value / 100,
    WAB_value = WAB_value / 100,
    logChron = log10(mri_chronicity_mo)
  )

ld1_spec <- "

    # LD1 is associated with handedness

    LD1 ~ coef_LV_LD*logLV + coef_LD_EHI*EHI_value
    WAB_value ~ coef_W_LD*LD1 + coef_W_LV*logLV

    # Indirect effect 
    mediation := coef_LV_LD*coef_W_LD

    total := (coef_LV_LD*coef_W_LD) + coef_W_LV

"

lm(LD1 ~ logLV + EHI_value, data = results_lhs_lavaan) %>%
  summary()

sem1 <- sem(ld1_spec, data = results_lhs_lavaan, se = "bootstrap")
summary(sem1)

lp1 <- lavaanPlot(sem1, coefs = TRUE, stars = c("regress", "latent", "covs"))
lp1

save_png(lp1, "sem_LD1.png", width = 500)

standardizedSolution(sem1)

# Interpret LD1-WAB
# Multiply by 100 to reverse WAB transformation above
effect_LD1_on_WAB <- coef(sem1)[names(coef(sem1)) == "coef_W_LD"] * 100
ld1_sd <- sd(results_lhs_lavaan$LD1)

one_SD_delta <- ld1_sd * effect_LD1_on_WAB

ld2_spec <- "

    LD2 ~ a*logLV
    # reserve c for EHI
    WAB_value ~ b*LD2 + d*logLV + e*logChron

    # Indirect effect 
    mediation := a*b

    total := (a*b) + d + e

"

sem2 <- sem(ld2_spec, data = results_lhs_lavaan, se = "bootstrap")
summary(sem2)

arrange(modificationIndices(sem2), desc(mi))

lp2 <- lavaanPlot(sem2, coefs = TRUE, stars = c("regress", "latent", "covs"))
lp2

save_png(lp2, "sem_LD2.png", width = 500)
