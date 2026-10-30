library(tidyverse)
library(reticulate)

setwd("~/MyDrive/Projects/hemisphere_fingerprinting_stroke/code/")

# Load HCP as comparison ======

hcp <- read_csv("../../hemisphere_fingerprinting/code/modeling/analyze_results/all_hemis_LDs.csv",
                show_col_types = FALSE) %>%
  select(-`...1`)

ggplot(hcp, aes(x = LD1, fill = class)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  labs(x = "Hemisphere axis")

# Stroke =====

# Load model to get % var explained
# py_install(c("pandas", "scipy", "matplotlib"), upgrade = TRUE)
m1 <- py_load_object("results/all/HC_group-all_method-lda_outcome-class_results.pickle")
m2 <- py_load_object("results/all/HC_group-all_method-lda_outcome-hemi_results.pickle")

# Demographics
demographics <- read_csv("../data/demographics.csv",
                          show_col_types = FALSE) %>%
  mutate(
    WAB_value = case_when(
      PatientStatus == "HC" ~ NA,
      .default = WAB_value
    )
  )

lesion_pct <- read_csv("../data/lesioned_pct.csv", show_col_types = FALSE)

pred_files <- list.files("predicted/", full.names = TRUE) %>%
  str_subset("Lpct", negate = TRUE)

predicted0 <- tibble(f = pred_files) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE),
    data = map(data, ~select(.x, -any_of("group")))
  )

predicted <- predicted0 %>%
  mutate(
    f = str_remove(f, "predicted/*")
  ) %>%
  separate_wider_delim(f, delim = "_", names = c("input", "model", NA),
                        cols_remove = FALSE) %>%
  unnest(data) %>%
  mutate(
    pred_patient = str_extract(class, "^[HCLS]*"),
    LD2 = jitter(replace_na(LD2, 0), amount = 1)
  ) %>%
  arrange(sub, input)

ggplot(predicted, aes(x = LD1, y = LD2, color = pred_patient)) +
  geom_point(alpha = 0.5) +
  facet_grid(cols = vars(input), rows = vars(model)) +
  theme_bw() +
  labs(x = "LD1 (hemi)", y = "LD2 (patient status)",
      caption = "LD2 for hemi model jittered between [-1, 1]")

ggsave("plots/all_inputs_both_models.png", width = 8.4, height = 4.0)

status_prediction <- predicted %>%
  filter(
    input == "input-all",
    model == "model-class",
    hemi == "RH"
  ) %>%
  separate_wider_delim(predicted, delim = "_",
                        names = c("pred_status", "pred_hemi"))

status_contigency <- table(status_prediction$PatientStatus,
  status_prediction$pred_status)

# HC change over age

hc <- predicted %>%
  filter(
    input == "input-HC"
  ) %>%
  select(sub, model, hemi, LD1) %>%
  pivot_wider(names_from = hemi, values_from = LD1) %>%
  mutate(
    distance = abs(L) + abs(R)
  ) %>%
  left_join(demographics)

hc %>%
  select(sub, model, distance) %>%
  pivot_wider(names_from = model, values_from = distance) %>%
  corrr::correlate()

ggplot(hc, aes(x = EHI_age, y = distance)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(vars(model)) +
  theme_bw()

lm(distance ~ EHI_age, data = filter(hc, model == "model-hemi")) %>%
  summary()

lm(distance ~ EHI_age, data = filter(hc, model == "model-class")) %>%
  summary()

library(caret)

sens <- sensitivity(as.factor(status_prediction$pred_status),
                    as.factor(status_prediction$PatientStatus))

spec <- specificity(as.factor(status_prediction$pred_status),
                    as.factor(status_prediction$PatientStatus))

mltools::mcc(preds = status_prediction$pred_status,
             actuals = status_prediction$PatientStatus)

# Model of interest

HCandLHS <- predicted %>%
  filter(
    model == "model-class",
    input == "input-all",
    hemi == "RH"
  ) %>%
  select(-PatientStatus, -LesionVolume) %>%
  left_join(demographics, by = join_by(sub)) %>%
  mutate(
    pred_group = str_extract(predicted, "^[LHSC]*")
  )

overall_accuracy <- sum(HCandLHS$class == HCandLHS$predicted) / nrow(HCandLHS)
patient_accuracy <- sum(str_detect(HCandLHS$class, HCandLHS$pred_group)) / nrow(HCandLHS)
hemi_accuracy <- sum(str_detect(HCandLHS$class, str_glue("{HCandLHS$hemi}H"))) / nrow(HCandLHS)


ggplot(HCandLHS, aes(x = EHI_age, y = LD2, color = PatientStatus)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = "Age (y)", title = "LD2 against age")

effsize::cohen.d(LD1 ~ PatientStatus, data = HCandLHS)
t.test(LD1 ~ PatientStatus, data = HCandLHS)

lm(LD1 ~ EHI_age, data = filter(HCandLHS, hemi == "RH", PatientStatus == "HC")) %>%
  summary()

lm(LD1 ~ log10(LesionVolume), data = HCandLHS) %>%
  summary()

lm(LD1 ~ WAB_value, data = HCandLHS) %>%
  summary()

lm(scale(WAB_value) ~ scale(log10(LesionVolume)), data = HCandLHS) %>%
  summary()

lm(scale(LD1) ~ scale(WAB_value), data = HCandLHS) %>%
  summary()

lm(LD1 ~ EHI_age*patient, data = HCandLHS) %>%
  summary()

ggplot(HCandLHS, aes(x = EHI_value, y = LD2, color = PatientStatus)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = "EHI", title = "LD2 against EHI")

lm(LD2 ~ EHI_value*PatientStatus, data = HCandLHS) %>%
  summary()

ggplot(HCandLHS, aes(x = WAB_value, y = LD2, color = PatientStatus)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = "WAB", title = "LD2 against WAB")

ggplot(HCandLHS, aes(x = WAB_value, y = LD1, color = PatientStatus)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = "WAB", title = "LD1 against WAB")

lm(LD1 ~ WAB_value, data = HCandLHS) %>%
  summary()

lm(LD2 ~ WAB_value, data = HCandLHS) %>%
  summary()

HC_with_WAB <- HCandLHS %>%
  filter(
    patient == "HC",
    !is.na(WAB_value)
  )

# Only RH cell

centroids <- HCandLHS %>%
  group_by(PatientStatus) %>%
  summarize(
    n = n(),
    mLD1 = mean(LD1),
    mLD2 = mean(LD2),
    sdLD1 = sd(LD1),
    sdLD2 = sd(LD2),
  ) %>%
  mutate(
    seLD1 = sdLD1 / sqrt(n),
    seLD2 = sdLD2 / sqrt(n)
  )

ggplot(HCandLHS, aes(color = PatientStatus)) +
  geom_point(aes(x = LD1, y = LD2), alpha = 0.5) +
  geom_smooth(aes(x = LD1, y = LD2), method = "lm") +
  geom_point(
    data = centroids,
    aes(x = mLD1, y = mLD2),
    shape = 15, size = 0.5
  ) +
  geom_errorbar(
    data = centroids,
    aes(y = mLD2, xmin = mLD1 - seLD1, xmax = mLD1 + seLD1),
    orientation = "y",width = 0.1
  ) +
  geom_errorbar(
    data = centroids,
    aes(x = mLD1, ymin = mLD2 - seLD2, ymax = mLD2 + seLD2),
    width = 0.1
  ) +
  theme_bw() +
  labs(x = "LD1", title = "LD1 against age")

# distance

HCandLHS_2 <- HCandLHS %>%
  mutate(

    LD1_dist = LD1 - centroids$mLD1[1],
    LD2_dist = LD2 - centroids$mLD2[1],

    # Distance from HC LD1 average; (-) means more ~LH and (+) means more RH
    LD1p = LD1 - centroids$mLD1[1],

    euclid_dist = sqrt(LD1_dist^2 + LD2_dist^2),

    # Weight according to the respective var explained ratios of the LDs
    euclid_dist2 = sqrt(m1$explained_variance_ratio_[1] * LD1_dist^2 +
                          m1$explained_variance_ratio_[2] * LD2_dist^2),

  )

ggplot(filter(HCandLHS_2, patient == "LHS"),
        aes(x = WAB_value, y = LD1p, color = Sex)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, color = "red") +
  scale_x_continuous(transform = "reverse") +
  theme_bw() +
  labs(x = "Worse aphasia", y = "LD1 - avg. HC LD1",
        title = "Worse aphasia leads to RH reinforcement")

ggplot(filter(HCandLHS_2, patient == "LHS"),
        aes(x = LesionVolume, y = LD1p, color = Sex)) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 0, color = "red") +
  theme_bw() +
  labs(x = "Lesion volume", y = "LD1 - avg. HC LD1",
        title = "Larger lesions lead to RH reinforcement")

lm_lesionvolume <- lm(LD1p ~ LesionVolume, data = HCandLHS_2) %>%
  broom::augment() %>%
  left_join(HCandLHS_2)

ggplot(lm_lesionvolume, aes(x = .std.resid, y = WAB_value)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()

ggplot(filter(HCandLHS_2, patient == "LHS"),
       aes(x = euclid_dist, y = WAB_value)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = "Distance from HC centroid", title = "LD2 against age")

ggplot(filter(HCandLHS_2, patient == "LHS"),
       aes(x = euclid_dist, y = WAB_value)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = "Distance from HC centroid", title = "LD2 against age")

ggplot(HCandLHS_2, aes(x = EHI_age, y = euclid_dist2, color = patient)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = "Age", title = "LD2 against age")

lm(euclid_dist ~ WAB_value, data = HCandLHS_2) %>%
  summary()

lm(euclid_dist2 ~ WAB_value, data = HCandLHS_2) %>%
  summary()

ggplot(filter(HCandLHS_2, patient == "LHS"),
        aes(x = WAB_chronicity / 12, y = euclid_dist2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = "Chronicity (years)", title = "LD2 against age")

ggplot(filter(HCandLHS_2, patient == "LHS"),
        aes(x = WAB_chronicity / 12, y = LD1p)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = "Chronicity (years)", title = "LD2 against age")

ggplot(filter(HCandLHS_2, patient == "LHS"),
       aes(x = log10(WAB_chronicity), y = LD1p)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw() +
  labs(x = "Chronicity (years)", title = "LD2 against age")

lm(LD1 ~ log10(WAB_chronicity), data = filter(HCandLHS_2, patient == "LHS")) %>%
  summary()

lm(WAB_value ~ log10(LesionVolume), data = HCandLHS) %>%
  summary()

lm(WAB_value ~ LD1, data = HCandLHS) %>%
  summary()

# Interactions ====

library(interactions)

HCandLHSx <- HCandLHS %>%
  mutate(
    WAB_value_Z = scale(WAB_value)[, 1],
    logLV_Z = scale(log10(LesionVolume))[, 1],
    LD1_Z = scale(LD1)[, 1]
  )

moderation_model <- lm(WAB_value_Z ~ logLV_Z*LD1_Z, data = HCandLHSx)

sim_slopes(moderation_model, pred = "logLV_Z", modx = "LD1_Z")

interact_plot(moderation_model, pred = "logLV_Z", modx = "LD1_Z",
              plot.points = TRUE)

# Lesion_pct -===

library(gamlss)
library(tidymodels)

lesioned <- HCandLHS %>%
  filter(
    PatientStatus == "LHS"
  ) %>%
  select(sub, LD1, LD2, LesionVolume) %>%
  left_join(lesion_pct) %>%
  group_by(ROI) %>%
  nest() %>%
  mutate(
    model = map(data, ~glmmTMB(LD1 ~ lesioned_pct,
                               data = .x,
                               ziformula = ~1) )
  )

# https://glmmtmb.github.io/glmmTMB/articles/glmmTMB.pdf
glmmTMB(LD1 ~ lesioned_pct,
        data = lesioned$data[[1]],
        ziformula = ~1) %>%
  tidy()

lesioned_results <- lesioned %>%
  mutate(
    broom = map(model, tidy)
  ) %>%
  select(-data, -model) %>%
  unnest(broom) %>%
  filter(
    # term == "lesioned_pct"
  )

# Language ROIs ====

rois_andrew <- readxl::read_excel("../data/BN_Atlas_274_combined.xlsx")

rois_L <- rois_andrew %>%
  filter(
    Hemisphere == "Left"
  ) %>%
  filter_out(
    Lobe %in% c("Subcortical", "Cerebellum")
  ) %>%
  select(LabelNoHemisphere, ChangNetwork) %>%
  mutate(
    LabelNoHemisphere = str_remove_all(LabelNoHemisphere, "[/+.]"),
    ChangNetwork = if_else(ChangNetwork != "NL", "LangN", "NL")
  )

lesion_pct2 <- left_join(lesion_pct, rois_L,
                         by = join_by(ROI == LabelNoHemisphere)) %>%
  mutate(
    lesioned = lesioned_pct > 0.1
  )

ggplot(filter(lesion_pct2, lesioned_pct > 0),
       aes(x = lesioned_pct, fill = ChangNetwork)) +
  geom_density()

lesions_count <- lesion_pct2 %>%
  summarize(
    .by = c(sub, ChangNetwork, lesioned),
    count = n()
  ) %>%
  complete(sub, ChangNetwork, lesioned) %>%
  mutate(
    count = replace_na(count, 0)
  )

lesioncount <- HCandLHS %>%
  filter(
    PatientStatus == "LHS"
  ) %>%
  select(sub, LD1, LD2, LesionVolume) %>%
  left_join(
    lesions_count %>%
      filter(lesioned) %>%
      select(-lesioned)
  ) %>%
  mutate(
    lesioned_regions = replace_na(count, 0)
  )

ggplot(lesioncount, aes(x = lesioned_regions, y = LD1, color = ChangNetwork)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()

lesioncount_wide <- lesioncount %>%
  pivot_wider(id_cols = c(sub, LD1, LD2),
              names_from = ChangNetwork, values_from = lesioned_regions) %>%
  mutate(
    across(c(NL, LangN), ~replace_na(.x, 0))
  ) %>%
  arrange(LangN, NL)

lm(LD1 ~ NL + LangN, data = lesioncount_wide) %>%
  summary()

lm(LD2 ~ NL + LangN, data = lesioncount_wide) %>%
  summary()

ggplot(lesioncount_wide, aes(x = LD1, y = LD2)) +
  geom_point(aes(color = NL)) +
  viridis::scale_color_viridis() +
  theme_bw()
