library(tidyverse)
library(ggnewscale)
library(viridis)
library(googlesheets4)

library(lavaan)
library(lavaanPlot)

library(lme4)
library(lmerTest)
library(glmnet)

# Setup ====

setwd("~/MyDrive/Projects/hemisphere_fingerprinting_stroke/code/")
source("../../hemisphere_fingerprinting/code/bootstrap_mcc.R")
source("acc_se.R")

tabulate_results <- function(df, ground, name, 
                             predicted = str_glue("{ground}_predicted")) {

  results <- df %>%
    group_by(!!sym(ground), !!sym(predicted))%>%
    summarize(
      n = n()
    ) %>%
    group_by(!!sym(ground)) %>%
    mutate(
      denom = sum(n),
      analysis = name,
      p = n / denom
    ) %>%
    rename(
      ground = !!sym(ground),
      predicted = !!sym(predicted)
    )

  return(results)

}

# Parallel setup ====

library(parallel)
library(doParallel)
library(foreach)

cores <- detectCores()
cl <- makeCluster(cores[1] - 1) 
registerDoParallel(cl)

# Demographics ====

annas_excluded_lh <- read_csv("lhs_only_in_Annas_projects.csv", 
                              show_col_types = FALSE)

demo0 <- read_csv("datasets/demographics.csv", show_col_types = FALSE) %>%
  filter_out(
    sub %in% annas_excluded_lh$sub
  )

patients <- unique(demo0$sub[demo0$PatientStatus == "LHS"])
conrols <- unique(demo0$sub[demo0$PatientStatus == "HC"])

# RH activation
rh_activ0 <- readxl::read_xlsx("../data/RHactivation_for_trevor.xlsx")

# Tidy up this activation table - give shorter names
rh_activ <- rh_activ0 %>%
  mutate(
    sub = str_glue("sub-{RegistryCode}")
  ) %>%
  select(sub, everything(), -Group, -RegistryCode) %>%
  select_all(~str_replace(tolower(.), "mean_activation", "activ")) %>%
  select_all(~str_replace(., "language", "lang"))

Hmisc::rcorr(as.matrix(select(rh_activ, -sub)))

# Lesion location 

lesion_loc0 <- read_csv("datasets/lesion_pct_by_division.csv", 
                        show_col_types = FALSE)

lesion_loc1 <- lesion_loc0 %>%
  filter(
    LR == "L"
  ) %>%
  pivot_wider(names_from = division, values_from = total_lesion_mm3,
              id_cols = c(sub, date)) %>%
  mutate(
    total_mm3 = replace_values(rowSums(select(., -sub, -date)), 0 ~ 1),
    # Pct of lesion lying outside language areas
    pct_lang = (Semantic + Language) / total_mm3,
    total_cm3 = total_mm3 * 0.001
  ) %>%
  left_join(
    distinct(select(demo0, sub, LesionVolume, logLV))
  ) %>%
  select(sub, date, 
    # my values:
    total_cm3, pct_lang, 
    # crl values
    LesionVolume, logLV
  )  %>%
  filter_out(
    sub %in% annas_excluded_lh$sub
  ) 

lesion_loc2 <- lesion_loc1 %>%
  group_by(sub) %>%
  summarize(
    total_cm3_mean = mean(total_cm3),
    pct_lang_mean = mean(pct_lang),
    LesionVolume = unique(LesionVolume),
    logLV = unique(logLV)
  )

patients %in% lesion_loc2$sub

# Check for missing lesion info
missing_lesion_info <- demo0 %>%
  filter(
    PatientStatus == "LHS"
  ) %>%
  left_join(lesion_loc) %>%
  filter(
    when_any(
      is.na(total_cm3_mean),
      is.na(pct_lang_mean),
      is.na(LesionVolume),
      is.na(logLV)
    )
  )

# My calculated value agains the canonical CRL value
# r = .978
cor.test(lesion_loc2$total_cm3_mean, lesion_loc2$LesionVolume)

# Is "percent of lesion in language network" correlated with volume?
# r = .301
cor.test(lesion_loc2$pct_lang_mean, lesion_loc2$LesionVolume)

residualize_pct_lang <- lm(pct_lang_mean ~ logLV, data = lesion_loc2)
summary(residualize_pct_lang)

ggplot(lesion_loc2, aes(x = LesionVolume, y = 100 * pct_lang_mean)) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_log10(limits = c(NA, NA), breaks = 10^(1:5), 
                labels = c("10", "100", "1,000", "10,000", "100,000")) +
  coord_cartesian(ylim = c(NA, NA)) +
  theme_bw() +
  labs(x = "Total lesion volume (cm^3)", y = "% of lesion in LN")

ggsave("plots/residualized_lnpct.png", width = 6.4, height = 4)

# Add the residalized valuaes
lesion_loc <- lesion_loc2 %>%
  add_column(residual_pct_lang = residualize_pct_lang$residuals)

demo <- left_join(
    demo0,
    select(lesion_loc, sub, pct_lang_mean, total_cm3_mean, residual_pct_lang),
    by = join_by(sub)
  ) %>%
  rename(
    date = ScanDate
  ) %>%
  left_join(
    # Only one value per ppt for now 
    rh_activ, 
    by = join_by(sub)
  )

## Relationship for lefties 

lm_logLV_logchron <- lm(WAB_value ~ logLV + log(mri_chronicity_mo), 
                        data = demo) 

lm_logLV_logchron_ehi <- lm(WAB_value ~ logLV + log(mri_chronicity_mo) + 
                              EHI_value, 
                              data = demo)

lm_logLV_logchron_ehi2 <- lm(WAB_value ~ logLV*EHI_value + 
                              log(mri_chronicity_mo),
                              data = demo)

AIC(lm_logLV_logchron, lm_logLV_logchron_ehi, lm_logLV_logchron_ehi2)
anova(lm_logLV_logchron, lm_logLV_logchron_ehi, lm_logLV_logchron_ehi2)

ggplot(filter(demo, PatientStatus == "LHS"), 
        aes(x = logLV, y = WAB_value, color = EHI_value > 0)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm") +
  coord_cartesian(ylim = c(NA, 100)) +
  theme_bw()

ggsave("plots/wab_on_LV_and_EHI.png", width = 6.5, height = 4)

# Load LDA results ====

results_testtrain <- read_csv(str_glue("predicted/",
                                       "input-all_outcome-class_predictedon-all_results.csv"),
                              show_col_types = FALSE) 


results_testonly <- read_csv(str_glue("predicted/",
                             "input-testonly_outcome-class_predictedon-all_results.csv"),
                    show_col_types = FALSE) 

# List of participants HC and LHS with longitudinal data
longitudinal_ppts <- unique(
    results_testonly$sub[results_testonly$status == "test_long"]
  )

results <- bind_rows(results_testtrain, results_testonly) %>%
  arrange(
    sub, date
  ) %>%
  group_by(sub) %>%
  mutate(
    days_since_first_scan = interval(min(date), date) / days(1),
    longitudinal = sub %in% longitudinal_ppts
  )


# Format this for Andrew's SVM GUI
for_svm <- results %>%
  filter(
    PatientStatus == "LHS"
  ) %>%
  mutate(
    RegistryCode = str_remove(sub, "sub-"),
  ) %>%
  ungroup() %>%
  select(RegistryCode, LD1, LD2) 

write_csv(for_svm, "datasets/for_svm.csv")

# Pull together results

results2 <- results %>%
  separate_wider_delim(
    class_predicted, delim = "_", 
    # Name these to remind that they are extracted from the three-way class
    # prediction
    names = c("class_pred_pt", "class_pred_hemi"),
    cols_remove = FALSE, too_few = "align_end"
  ) %>%
  mutate(
    hemi = str_glue("{hemi}H"),
    ground_pred_hemi = str_glue("{hemi} => {class_pred_hemi}"),
    correct_hemi = if_else(hemi == class_pred_hemi, "correct", "incorrect")
  ) %>%
  select(sub, date, hemi, PatientStatus, status, group, class, fold, 
          LD1, LD2, class_pred_hemi, correct_hemi, class_predicted, 
          ground_pred_hemi,
          logLV_g, WAB_value_g, residual_wab, residual_wab_std, 
          residual_wab_group, days_since_first_scan, longitudinal) %>%
  left_join(
    # Many of these columns end up in both the demographics and then in the 
    # processed data. This should be 'one-to-one' but the HC participants have
    # two hemi rows
    demo, 
    relationship = "many-to-one",
    by = join_by(sub, date, PatientStatus, group)
  ) 

# Basic findings

table(results2$hemi, results2$class_pred_hemi)
table(results2$class, results2$class_predicted)

acc_se(sum(results2$hemi == results2$class_pred_hemi), nrow(results2)) 
bootstrap_mcc(results2, ground = "hemi", predicted = "class_pred_hemi")

acc_se(sum(results2$class == results2$class_predicted), nrow(results2)) 
bootstrap_mcc(results2, ground = "class", predicted = "class_predicted")

results2_plot <- tabulate_results(results2, ground = "class", 
                                  name = "class") %>%
  ungroup() %>%
  mutate(
    across(c(ground, predicted), ~str_replace(., "_", " "))
  )

ggplot(results2_plot, aes(x = predicted, y = ground, fill = p)) +
  geom_tile() +
  geom_label(aes(label = str_glue("{n}\n{round(p, 2)}")), fill = "white",
             alpha = 0.75) +
  scale_x_discrete(position = "top", limits = unique(results2_plot$ground)) +
  scale_y_discrete(limits = rev(unique(results2_plot$ground))) +
  viridis::scale_fill_viridis(limits = c(0, 1)) +
  coord_equal() +
  theme_bw() +
  labs(x = "Predicted", y = "Ground")

ggsave("plots/fig3-allin_class.png", height = 3.2)

write_csv(results2, "datasets/ld1ld2_results.csv")

# Plot 1 ======

# Remove the longitudinal participants and order the rows so that the test
# rows end up on top
results2_nolong_forplot <- results2 %>%
  filter_out(
    status == "test_long"
  ) %>%
  arrange(desc(status))

ggplot(NULL, aes(x = LD1, y = LD2)) +
  geom_point(
    data = filter(results2_nolong_forplot, PatientStatus == "HC"),
    aes(fill = hemi, shape = ground_pred_hemi),
    size = 2.5, alpha = 0.75
  ) +
  scale_shape_manual(
    values = c("LH => LH" = 21, "LH => RH" = 24, "RH => RH" = 22, 
               "RH => LH" = 25)
  ) +
  scale_fill_discrete(name = "HC Hemi.") + 
  guides(fill = guide_legend(override.aes = list(size = 3, stroke = 1))) +
  new_scale_fill() +
  geom_point(
    data = filter(results2_nolong_forplot, PatientStatus == "LHS"),
    aes(fill = logLV, shape = ground_pred_hemi, color = status),
    size = 3, stroke = 1
  ) +
  scale_fill_gradient2(
    midpoint = median(results2_nolong_forplot$logLV, na.rm = TRUE),
    low = "red", mid = "white", high = "blue"
  ) +
  scale_color_manual(
    values = c(train_and_test = "black", test_only = "red", 
                test_rhlesion = "red")
  ) +
  theme_bw() +
  labs(fill = "log(lesion vol.)", shape = "Hemi: G => P") +
  guides(shape = guide_legend(override.aes = list(size = 3, stroke = 1, 
                                                  fill = "grey"))) 

# Plot 2 ====

p2 <- ggplot(NULL, aes(x = LD1, y = LD2)) +
  theme_bw() +
  geom_point(
    data = filter(results2_nolong_forplot, PatientStatus == "HC"),
    aes(color = hemi, shape = correct_hemi),
    size = 2.5, alpha = 0.75
  ) + 
  scale_color_discrete(name = "True\nsidedness\n(HC)") +
  scale_shape_manual(
    values = c("correct" = 3, "incorrect" = 4),
    name = "Prediction\n(HC)"
  ) +
  new_scale("shape") +
  new_scale("color") +
  geom_point(
    data = filter(results2_nolong_forplot, PatientStatus == "LHS"),
    aes(fill = logLV, color = status),
    shape = 22, size = 3, stroke = 1
  ) +
  scale_fill_viridis(name = expression(LHS:log[10](LV)))

p2 +
  scale_color_manual(
    values = c(train_and_test = "black", test_only = "red", 
                test_rhlesion = "red")
  )

centroids_all <- results2 %>%
  group_by(class, status) %>%
  summarize(
    n = n(),
    across(c(LD1, LD2), .fns = c(m = mean, sd = sd), .names = "{.col}_{.fn}"),
    across(ends_with("_sd"), ~ .x / sqrt(n), 
            .names = "{str_replace(.col, 'sd', 'se')}")
  )

write_csv(centroids_all, "datasets/centroids_all.csv")

# I only really want to plot these ones
centroids <- centroids_all %>%
  filter(
    status == "train_and_test"
  )

p2 + 
  new_scale("fill") +
  geom_point(
    data = centroids,
    aes(x = LD1_m, y = LD2_m, fill = class),
    size = 4.5, shape = 23, stroke = 2.5, color = "black"
  ) +
  scale_fill_manual(values = c("HC_LH" = "#F8766D", "HC_RH" = "#00BFC4",
                               "LHS_RH" = "gold"),
                      name = "Centroids", guide = "none") +
  scale_color_manual(
    values = c(test_and_train = "black", test_only = "red", 
    test_rhlesion = "red")
  ) +
  theme(legend.position = "bottom") +
  guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

ggsave("plots/fig4-hclhs_ld1_ld2.png", width = 6.5, height = 5)

# poster_plot <- ggplot(NULL, aes(x = LD1, y = LD2)) +
#   geom_point(
#     data = filter(results2, PatientStatus == "HC"),
#     aes(color = hemi, shape = correct_hemi),
#     size = 6, alpha = 0.75, stroke = 2
#   ) + 
#   scale_color_discrete(name = "Hemi") +
#   scale_shape_manual(
#     values = c("correct" = 3, "incorrect" = 4),
#     labels = c("Correct", "Incorrect"),
#     name = "Prediction"
#   ) +
#   guides(
#     # Put this up here before creating new scales
#     #color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)
#   ) +
#   new_scale("shape") +
#   new_scale("color") +
#   geom_point(
#     # LHS squares
#     data = filter(results2, PatientStatus == "LHS"),
#     aes(fill = logLV),
#     shape = 22, size = 8, stroke = 1
#   ) +
#   scale_fill_viridis(name = expression(log[10](LV)), 
#                      guide = guide_colorbar(direction = "vertical",
#                                             position = "right")) + 
#   new_scale("fill") +
#   geom_point(
#     # Centrois
#     data = filter_out(centroids, status == "test_only"),
#     aes(x = LD1_1, y = LD2_1, fill = class),
#     size = 10, shape = 23, stroke = 2.5,
#   ) +
#   scale_fill_manual(values = c("HC_LH" = "#F8766D", "HC_RH" = "#00BFC4",
#                                "LHS_RH" = "gold"),
#                       name = "Centroids", guide = "none") +
#   theme_minimal() +
#   theme(legend.position = "bottom", text = element_text(size = 24)) +
#   labs(x = "LD1 (hemisphere axis)", y = "LD2 (stroke axis)")

# poster_plot

# ggsave(poster_plot, filename = "plots/ld1_ld2_poster.png", width = 14, 
#        height = 8)

# mean(results2$LD1[results2$PatientStatus == "HC" & results2$hemi == "RH"])

# # statement plot
# ggplot(NULL, aes(x = LD1, y = LD2)) +
#   geom_point(
#     data = filter(results2, PatientStatus == "HC", hemi == "RH"),
#     aes(color = hemi),
#     size = 2, alpha = 0.5, color = "orange"
#   ) + 
#   scale_x_continuous(breaks = c(0, 1), labels = c(0, "1 (RH)")) +
#   scale_y_continuous(breaks = c(-2, 0, 4), 
#                       labels = c("Smaller\nlesions", "", "Larger\nlesions")) +
#   new_scale("shape") +
#   new_scale("color") +
#   geom_point(
#     # LHS squares
#     data = filter(results2, PatientStatus == "LHS"),
#     aes(fill = logLV),
#     shape = 22, size = 3, stroke = 1
#   ) +
#   scale_fill_viridis() + 
#   theme_classic() +
#   theme(text = element_text(size = 10.5), legend.position = "none") +
#   labs(x = "Hemisphere axis", y = "Stroke axis")

# ggsave("plots/res_statement.png", width = 2.75, height = 1.75)
