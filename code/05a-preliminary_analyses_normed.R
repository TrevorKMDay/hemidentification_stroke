library(tidyverse)

# Setup =====

setwd("/Users/tkmd/MyDrive/Projects/hemisphere_fingerprinting_stroke/code")

source("../../hemisphere_fingerprinting/code/bootstrap_mcc.R")
source("acc_se.R")

# Parallel setup ====

library(parallel)
library(doParallel)
library(foreach)

cores <- detectCores()
cl <- makeCluster(cores[1] - 1) #not to overload your computer
registerDoParallel(cl)

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

# Predicting hemisphere =====

hemi_hconly0 <- read_csv("results/loocv/input-HC_outcome-hemi_loocv.csv",
                           show_col_types = FALSE) %>%
  select(-index)

hemi_hconly <- hemi_hconly0  %>%
  select(sub, hemi, hemi_predicted, LD1) %>%
  mutate(
    correct = if_else(hemi == hemi_predicted, 1, 0)
  )

hemi_hconly %>%
  filter(
    !correct
  ) %>%
  arrange(sub)

hemi_hconly %>%
  select(-hemi_predicted) %>%
  pivot_wider(names_from = hemi, values_from = c(LD1, correct)) %>%
  ggplot(aes(x = LD1_L, y = LD1_R, color = interaction(correct_L, correct_R))) +
    geom_point() +
    theme_bw()

table(hemi_hconly$hemi, hemi_hconly$hemi_predicted)

hemi_hconly_acc <- acc_se(sum(hemi_hconly$correct), n = nrow(hemi_hconly))

hemi_hconly_mcc <- bootstrap_mcc(hemi_hconly, ground = "hemi", 
                                 predicted = "hemi_predicted")

## Predict LHS hemisphere status from hc only model

hemi_LHSonHC_0 <- read_csv("predicted/input-hc_outcome-hemi_predictedon-lhs_results.csv",
                           show_col_types = FALSE)

hemi_LHSonHC <- hemi_LHSonHC_0 %>%
  select(sub, hemi, logLV, hemi_predicted, LD1) %>%
  mutate(
    logLV_Z = scale(logLV)[, 1],
    correct = if_else(hemi == hemi_predicted, 1, 0)
  )

hemi_LHSonHC %>%
  filter(
    correct == 0
  )

table(hemi_LHSonHC$hemi, hemi_LHSonHC$hemi_predicted)

hemi_acc <- acc_se(sum(hemi_LHSonHC$correct), n = nrow(hemi_LHSonHC))

binom_LHSonHC <- glm(correct ~ logLV, data = hemi_LHSonHC, family = binomial)
summary(binom_LHSonHC)
                          
# HC + LHS

hemi0 <- read_csv("results/loocv/input-all_outcome-hemi_loocv.csv",
                           show_col_types = FALSE) 

hemi <- hemi0  %>%
  select(sub, PatientStatus, hemi, hemi_predicted, LD1) %>%
  mutate(
    correct = if_else(hemi == hemi_predicted, 1, 0)
  )

hemi %>%
  filter(
    !correct
  ) %>%
  arrange(sub)

hemi %>%
  select(-hemi_predicted) %>%
  pivot_wider(names_from = hemi, values_from = c(LD1, correct)) %>%
  ggplot(aes(x = LD1_L, y = LD1_R, color = interaction(correct_L, correct_R))) +
    geom_point() +
    theme_bw()

hemi_acc <- acc_se(sum(hemi$correct), n = nrow(hemi))
hemi_mcc <- bootstrap_mcc(hemi, ground = "hemi", predicted = "hemi_predicted")

table(hemi$hemi, hemi$hemi_predicted)

filter(hemi, hemi != hemi_predicted) %>%
  left_join(
    distinct(select(hemi0, sub, PatientStatus, EHI_value)),
    by = join_by(sub)
  )

hemi_results <- hemi %>%
  group_by(hemi, hemi_predicted) %>%
  summarize(
    n = n()
  ) %>%
  group_by(hemi) %>%
  mutate(
    denom = sum(n),
    analysis = "Hemisphere",
    p = n / denom
  ) %>%
  rename(
    ground = hemi,
    predicted = hemi_predicted
  )

# Predicting patient status =====

patient_rhonly0 <- read_csv("results/loocv/input-allRHonly_outcome-patient_loocv.csv",
                           show_col_types = FALSE) 

patient_rhonly <- patient_rhonly0 %>%
  select(sub, LesionVolume, logLV, PatientStatus, PatientStatus_predicted,
          LD1) %>%
  mutate(
    correct = if_else(PatientStatus == PatientStatus_predicted, 1, 0)
  )

patient_rhonly_acc <- acc_se(sum(patient_rhonly$correct),
                             n = nrow(patient_rhonly))

patient_rhonly_mcc <- bootstrap_mcc(patient_rhonly, ground = "PatientStatus",
                                    predicted = "PatientStatus_predicted")

patient_contingency <- table(patient_rhonly$PatientStatus, 
                             patient_rhonly$PatientStatus_predicted)

patient_rhonly_results <- patient_rhonly %>%
  group_by(PatientStatus, PatientStatus_predicted) %>%
  summarize(
    n = n()
  ) %>%
  group_by(PatientStatus) %>%
  mutate(
    denom = sum(n),
    analysis = "Patient Status",
    p = n / denom
  ) %>%
  rename(
    ground = PatientStatus,
    predicted = PatientStatus_predicted
  )

patient_rhonly_patientsonly <- filter(patient_rhonly, PatientStatus == "LHS")
patient_rhonly_HConly <- filter(patient_rhonly, PatientStatus == "HC")

acc_se(sum(patient_rhonly_patientsonly$correct), n = nrow(patient_rhonly_patientsonly))
acc_se(sum(patient_rhonly_HConly$correct), n = nrow(patient_rhonly_HConly))

binom1 <- glm(correct ~ LesionVolume, data = patient_rhonly_patientsonly, 
              family = binomial)
summary(binom1)

# This is the model that shows that larger volumes are more likely to be
# correctly predicted
binom1_log <- glm(correct ~ logLV, data = patient_rhonly_patientsonly, 
                  family = binomial)
summary(binom1_log)
exp(coef(binom1_log))

# Check against orginal subs
original <- read_csv("../data/timestamps_260512.csv") %>%
  filter(
    Unit == "fc"
  )

patient_rhonly_old <- filter(patient_rhonly, sub %in% str_glue("sub-{original$SubjectCode}"))

binom1_old <- glm(correct ~ LesionVolume, data = patient_rhonly_old, 
                    family = binomial)
summary(binom1_old)


x_at_prob <- function(model, prob) {

  lhs <- log((prob) / (1 - prob))
  int <- coef(model)[1]
  Lcoef <- coef(model)[2]

  x <- (lhs - int) / Lcoef
  
  return(x)

}

# Set limits to up to 90% accuracy
# x_at_prob(binom1_log, 0.9)
# scale_x_continuous(limits = c(NA, x_at_prob(binom1_log, 0.9))) +

ggplot(patient_rhonly_patientsonly, aes(x = LesionVolume, y = correct)) +
  geom_jitter(width = 0, height = 0.025) +
  stat_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = TRUE, fullrange = TRUE,
              color = "blue", 
              fill = "lightblue") +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25)) +
  theme_bw() +
  labs(x = expression(log[10](LesionVolume)), 
        y = "Prob of being correctly classified")

correct1 <- ggplot(patient_rhonly_patientsonly, aes(x = logLV, y = correct)) +
  geom_jitter(width = 0, height = 0.025, alpha = 0.5) +
  stat_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = TRUE, fullrange = TRUE,
              color = "blue", 
              fill = "lightblue") +
  scale_y_continuous(limits = c(-0.1, 1.1), breaks = seq(0, 1, 0.25)) +
  theme_bw() +
  labs(x = expression(log[10](LesionVolume)), 
        y = "Prob of being correctly classified",
       title = "(A) RH only")

correct1 

ggsave("plots/correct_prediction_RHonly.png", width = 3.25, height = 4)

# both hemispheres

patient_class0 <- read_csv("results/loocv/input-all_outcome-class_loocv.csv",
                           show_col_types = FALSE) %>%
  select(-index)

patient_class <- patient_class0 %>%
  select(sub, LesionVolume, logLV, class, class_predicted,  LD1) %>%
  mutate(
    correct = if_else(class == class_predicted, 1, 0)
  )

patient_class_acc <- acc_se(sum(patient_class$class == patient_class$class_predicted),
                            n = nrow(patient_class))

patient_class_mcc <- bootstrap_mcc(patient_class, ground = "class",
                                    predicted = "class_predicted")

patient_class_mcc2 <- bootstrap_mcc(filter_out(patient_class, class == "HC_LH"), 
                                    ground = "class", 
                                    predicted = "class_predicted")

table(patient_class$class, patient_class$class_predicted)

patient_class_results <- patient_class %>%
  group_by(class, class_predicted) %>%
  summarize(
    n = n()
  ) %>%
  group_by(class) %>%
  mutate(
    denom = sum(n),
    analysis = "Patient Status",
    p = n / denom
  ) %>%
  rename(
    ground = class,
    predicted = class_predicted
  )

# binom1_class <- glm(correct ~ LesionVolume, data = patient_class, 
#   family = binomial)
# summary(binom1_class)

# exp(coef(binom1_class))

# Larger volumes 
binom1_log_class <- glm(correct ~ logLV, 
                        data = filter_out(patient_class, 
                                          str_detect(class, "^HC_")), 
                        family = binomial)

summary(binom1_log_class)

exp(coef(binom1_log_class))

binomial_plot_data <- bind_rows(patient_rhonly, patient_class) %>%
  filter(
    PatientStatus == "LHS" | class == "LHS_RH"
  ) %>%
  mutate(
    analysis = if_else(!is.na(PatientStatus_predicted), "rh_only", "class")
  ) %>%
  select(sub, logLV, correct, analysis) 

ggplot(NULL, aes(x = logLV, y = correct)) +
  geom_point(data = filter(binomial_plot_data, analysis == "rh_only"),
              position = position_jitter(height = 0.05)) +
  stat_smooth(data = binomial_plot_data, 
              aes(color = analysis, fill = analysis),
              method = "glm", 
              method.args = list(family = "binomial"), 
              se = TRUE, fullrange = TRUE, linewidth = 1.5) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25)) +
  theme_bw()

ggsave("plots/correct_prediction_2and3way.png", width = 6.5, height = 4)

correct2 <- ggplot(filter(patient_class, str_detect(class, "LHS")), 
       aes(x = logLV, y = correct)) +
  geom_jitter(width = 0, height = 0.025, alpha = 0.5) +
  stat_smooth(method = "glm", 
              method.args = list(family = "binomial"), 
              se = TRUE, fullrange = TRUE,
              color = "blue", 
              fill = "lightblue") +
  scale_y_continuous(limits = c(-0.1, 1.1), breaks = seq(0, 1, 0.25)) +
  theme_bw() +
  labs(x = expression(log[10](LesionVolume)), 
       y = "Prob of being correctly classified",
       title = "(B) Three-way prediction")

correct2

ggsave("plots/correct_prediction_class.png", width = 3.25, height = 4)

correct1 + correct2 + plot_layout(axes = "collect")

ggsave("plots/correct_prediction_hemi_class.png", width = 6.5, height = 4)

## Cumulative probability

cum_prob <- binomial_plot_data %>%
  arrange(logLV) %>%
  group_by(analysis) %>%
  mutate(
    n = row_number(),
    sum_correct = cumsum(correct),
    p_correct = sum_correct / n,
    p_correct_se = (p_correct * (1 - p_correct)) / sqrt(n)
  )

ggplot(cum_prob, aes(x = logLV, y = p_correct, color = analysis)) +
  geom_pointrange(aes(ymin = p_correct - p_correct_se, 
                      ymax = p_correct + p_correct_se)) +
  geom_line() +
  theme_bw()

## Figure 1

fig1a <- ggplot(hemi_results, aes(x = predicted, y = ground, fill = p)) +
  geom_tile() +
  geom_label(aes(label = str_glue("{n}\n{round(p, 2)}")), fill = "white", 
             alpha = 0.75) +
  scale_x_discrete(position = "top", breaks = NULL) +
  scale_y_discrete(limits = c("R", "L"), breaks = NULL) +
  viridis::scale_fill_viridis(limits = c(0, 1)) +
  coord_equal() +
  theme_bw() +
  labs(x = "Predicted", y = "Ground", title = "(A) Hemisphere")

fig1b <- ggplot(patient_rhonly_results, 
                aes(x = predicted, y = ground, fill = p)) +
  geom_tile() +
  geom_label(aes(label = str_glue("{n}\n{round(p, 2)}")), fill = "white", 
             alpha = 0.75) +
  scale_x_discrete(position = "top") +
  viridis::scale_fill_viridis(limits = c(0, 1)) +
  coord_equal() +
  theme_bw() +
  labs(x = "Predicted", y = "Ground", title = "(A) Patient status (RH only)")

patient_class_results <- patient_class_results %>%
  ungroup() %>%
  mutate(
    across(c(ground, predicted), ~str_replace(.x, "_", " "))
  )

fig1c <- ggplot(patient_class_results, 
                aes(x = predicted, y = ground, fill = p)) +
  geom_tile() +
  geom_label(aes(label = str_glue("{n}\n{round(p, 2)}")), fill = "white") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = c("LHS RH", "HC RH", "HC LH")) +
  viridis::scale_fill_viridis(limits = c(0, 1)) +
  coord_equal() +
  theme_bw() +
  labs(x = "Predicted", y = "Ground", title = "(B) Patient status + hemi.") +
  annotate(
    geom = "rect",
    xmin = 1.5, xmax = 3.5,
    ymin = 0.5, ymax = 2.5,
    fill = NA, color = "red",
    linewidth = 3
  )

library(patchwork)

(fig1b + fig1c) + 
  plot_layout(axes = "collect", guides = "collect") 

ggsave("plots/fig1-hemi_patient_class.png", width = 6.5)

# Log LV =====

logLV0 <- read_csv("results/loocv/input-lhs_outcome-lLVg_loocv.csv",
                           show_col_types = FALSE) 

sort(unique(logLV0$logLV_g))
ordered_range <- c("[-0.00922,1.51]", "(1.51,2.03]", "(2.03,2.32]", 
                   "(2.32,2.71]")

logLV <- logLV0 %>%
  select(sub, LesionVolume, logLV, logLV_g, logLV_g_predicted,  LD1) %>%
  mutate(
    correct = if_else(logLV_g == logLV_g_predicted, 1, 0),
    across(starts_with("logLV_g"), ~factor(.x, levels = ordered_range))
  )

logLV %>%
  group_by(logLV_g) %>%
  summarize(
    n = n(),
    min = min(logLV),
    max = max(logLV)
  )

# ggplot(logLV, aes(x = LesionVolume, y = logLV_g)) +
#   geom_point()

logLV_acc <- acc_se(sum(logLV$logLV_g == logLV$logLV_g_predicted),
                            n = nrow(logLV))

logLV_mcc <- bootstrap_mcc(logLV, ground = "logLV_g",
                                    predicted = "logLV_g_predicted")


# WAB

WAB0 <- read_csv("results/loocv/input-lhs_outcome-WABg_loocv.csv",
                           show_col_types = FALSE)

WAB <- WAB0 %>%
  select(sub, LesionVolume, WAB_value, WAB_value_g, WAB_value_g_predicted,  LD1) %>%
  mutate(
    correct = if_else(WAB_value_g == WAB_value_g_predicted, 1, 0),
  )

ggplot(WAB, aes(x = WAB_value, y = WAB_value_g)) +
  geom_point()

WAB_acc <- acc_se(sum(WAB$WAB_value_g == WAB$WAB_value_g_predicted), n = nrow(WAB))

WAB_mcc <- bootstrap_mcc(WAB, ground = "WAB_value_g",
                         predicted = "WAB_value_g_predicted")

wab_table <- table(WAB$WAB_value_g, WAB$WAB_value_g_predicted)
wab_cols <- margin.table(table(WAB$WAB_value_g, WAB$WAB_value_g_predicted), margin = 2)

for (i in 1:4) {
  print(round(100 * wab_table[i, ] / wab_cols[i]))
}

# rzd WAB

rWAB0 <- read_csv("results/loocv/input-lhs_outcome-rWABg_loocv.csv",
                           show_col_types = FALSE) 

rWAB <- rWAB0 %>%
  select(sub, LesionVolume, residual_wab, starts_with("resid"),  LD1) %>%
  mutate(
    correct = if_else(residual_wab == residual_wab_group_predicted, 1, 0),
    across(starts_with("residual_wab_group"), 
            ~factor(.x, levels = c("rzd_wab_worse", "rzd_wab_avg", 
                                    "rzd_wab_better")))
  )

ggplot(rWAB, aes(x = residual_wab, y = residual_wab_group)) +
  geom_point()

rWAB_acc <- acc_se(sum(rWAB$residual_wab_group == rWAB$residual_wab_group_predicted), 
                    n = nrow(rWAB))

rWAB_mcc <- bootstrap_mcc(rWAB, ground = "residual_wab_group",
                          predicted = "residual_wab_group_predicted")

rwab_table <- table(rWAB$residual_wab_group, rWAB$residual_wab_group_predicted)

rwab_cols <- margin.table(table(rWAB$residual_wab_group, 
                                rWAB$residual_wab_group_predicted), 
                          margin = 1)

for (i in 1:3) {
  print(round(100 * rwab_table[i, ] / rwab_cols[i]))
}

# LogLV and WAB results

logLV_results <- tabulate_results(logLV, ground = "logLV_g", 
                                  name = "log(LV)") %>%
  ungroup() %>%
  mutate(
    across(c(ground, predicted), ~str_replace(.x, ",", ",\n"))
  )

WAB_results <- tabulate_results(WAB, ground = "WAB_value_g", name = "WAB")  %>%
  ungroup() %>%
  mutate(
    across(c(ground, predicted), ~str_replace(.x, ",", ",\n"))
  )

label_size <- 3

fig2a <- ggplot(logLV_results, aes(x = predicted, y = ground, fill = p)) +
  geom_tile() +
  geom_label(aes(label = str_glue("{n}\n{round(p, 2)}")), fill = "white",
              size = label_size, alpha = 0.75) +
  scale_x_discrete(position = "top", limits = rev(unique(logLV_results$ground))) +
  scale_y_discrete(limits = unique(logLV_results$ground)) +
  viridis::scale_fill_viridis(limits = c(0, 1)) +
  coord_equal() +
  theme_bw() +
  labs(x = "Predicted", y = "Ground", 
        title = expression("(A)"~log[10](LV)~"group"))

fig2b <- ggplot(WAB_results, aes(x = predicted, y = ground, fill = p)) +
  geom_tile() +
  geom_label(aes(label = str_glue("{n}\n{round(p, 2)}")), fill = "white",
              size = label_size, alpha = 0.75) +
  scale_x_discrete(position = "top", limits = unique(WAB_results$ground)) +
  scale_y_discrete(limits = rev(unique(WAB_results$ground))) +
  viridis::scale_fill_viridis(limits = c(0, 1)) +
  coord_equal() +
  theme_bw() +
  labs(x = "Predicted", y = "Ground", title = "(B) WAB group")


fig2a + fig2b + plot_layout(guides = "collect")

ggsave("plots/fig2-logLV_WAB.png", width = 6.5)

# All MCC 

all_mcc <- tribble(
    ~comparison, ~mcc,
    "hemi_hc", hemi_hconly_mcc,
    "hemi_all", hemi_mcc,
    "patient_rh", patient_rhonly_mcc,
    "patient_all", patient_class_mcc,
    "patient_all_rhonly", patient_class_mcc2,
    "logLV", logLV_mcc,
    "WAB", WAB_mcc,
    "rWAB", rWAB_mcc
  ) %>%
  unnest(mcc)

ggplot(all_mcc, aes(x = comparison)) +
  geom_pointrange(aes(ymin = mcc_ci_lo, y = mcc_mean, ymax = mcc_ci_hi)) +
  geom_hline(yintercept = 0, color = "red") +
  theme_bw()

# ln network

ln_damage <- read_csv("results/loocv/input-lhs_outcome-LNyn_loocv.csv",
                           show_col_types = FALSE) 

ln_damage_acc <- acc_se(
    sum(ln_damage$is_ln_lesioned == ln_damage$is_ln_lesioned_predicted), 
    n = nrow(ln_damage)
  )

ln_damage_mcc <- bootstrap_mcc(ln_damage, ground = "is_ln_lesioned",
                                predicted = "is_ln_lesioned_predicted")

# aphasia

aphasia <- read_csv("results/loocv/input-lhs_outcome-aphasia_loocv.csv",
                           show_col_types = FALSE) 

aphasia_acc <- acc_se(
    sum(aphasia$aphasia == aphasia$aphasia_predicted), 
    n = nrow(aphasia)
  )

aphasia_mcc <- bootstrap_mcc(aphasia, ground = "aphasia",
                                predicted = "aphasia_predicted")
