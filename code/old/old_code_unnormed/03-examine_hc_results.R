library(tidyverse)

library(doParallel)
library(foreach)

select <- dplyr::select

# Parallel setup ====

cores <- detectCores()
cl <- makeCluster(cores[1] - 1) # not to overload your computer
registerDoParallel(cl)

# Boostrapping code ====

setwd("~/MyDrive/Projects/hemisphere_fingerprinting_stroke/code/")
source("../../hemisphere_fingerprinting/code/modeling/analyze_results/bootstrap_mcc.R")

# Get data ====

demographics <- read_csv("../data/demographics.csv", show_col_types = FALSE)

result_files <- list.files("results/", "*.csv", full.names = TRUE,
                           recursive = TRUE) %>%
  str_subset("HC")

results0 <- tibble(f = result_files) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE, .progress = TRUE)
  )

# Predicting patient status
patient_files <- list.files("results/", "*.csv", full.names = TRUE,
                           recursive = TRUE) %>%
  str_subset("outcome-patient")

patients0 <- tibble(f = patient_files) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE, .progress = TRUE)
  )

## regular methods ====

results <- results0 %>%
  filter_out(
    str_detect(f, "loocv")
  ) %>%
  unnest(data) %>%
  select(-`...1`) %>%
  mutate(
    held_out = str_extract(f, "group-[^_]*"),
    method = str_extract(f, "method-[^_]*") %>%
      str_remove("method-")
  ) %>%
  select(-f) %>%
  arrange(sub, method, ground)

results_mcc <- results %>%
  group_by(method, held_out) %>%
  nest() %>%
  ungroup() %>%
  mutate(
    n = map_int(data, nrow),

    mcc = map(data, ~bootstrap_mcc(data = .x, R = 1000, method = 2)),
    accuracy = map_dbl(data, ~sum(.x$ground == .x$predicted) / nrow(.x)),
    accuracy_se = sqrt(accuracy * (1 - accuracy) / n)
  )

results_mcc2 <- results_mcc %>%
  select(-data) %>%
  unnest(mcc)

ggplot(results_mcc2 , aes(x = method, y = mcc_mean, color = method,
  shape = held_out)) +
  geom_pointrange(aes(ymin = mcc_ci_lo, ymax = mcc_ci_hi),
                  position = position_dodge2(width = 0.25)) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, by = 0.25)) +
  theme_bw() +
  labs(x = "Method", y = "MCC")

ggplot(results_mcc2 , aes(x = method, y = accuracy, color = method,
                          shape = held_out)) +
  geom_pointrange(aes(ymin = accuracy - accuracy_se,
                      ymax = accuracy + accuracy_se),
                  position = position_dodge2(width = 0.25)) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, by = 0.25)) +
  theme_bw() +
  labs(x = "Method", y = "Accuracy")

# Patients ====

patients <- patients0 %>%
  unnest(data) %>%
  select(-`...1`) %>%
  mutate(
    held_out = str_extract(f, "group-[^_]*"),
    method = str_extract(f, "method-[^_]*") %>%
      str_remove("method-")
  ) %>%
  select(-f) %>%
  arrange(sub, method, ground) %>%
  left_join(demographics) %>%
  mutate(

    correct_patient = if_else(ground == predicted, 1, 0)

  )

patients_lda <- patients %>%
  filter(
    ground == "LHS",
    method == "lda"
  ) %>%
  mutate(
    log10LV = log10(LesionVolume)
  )

# Are larger strokes more likely to be identified as controls?
lv_log <- glm(correct_patient ~ log10LV, data = patients_lda,
              family = "binomial")

exp(coef(lv_log))

optimize(lv_log, y = 0.5)

predict_x <- seq(0, 3, by = 0.1)
lv_log_predicted <- tibble(
  x = predict_x,
  y_hat = predict(lv_log, newdata = tibble(log10LV = predict_x),
                  type = "response")
)

ggplot(patients_lda, aes(x = log10(LesionVolume), y = correct_patient)) +
  geom_jitter(width = 0, height = 0.05, alpha = 0.5) +
  geom_line(
    data = lv_log_predicted,
    aes(x = x, y = y_hat)
  ) +
  theme_bw() +
  labs(title = "Correct prediction of patient status from RH connectivity")

patients_mcc <- patients %>%
  group_by(method) %>%
  nest() %>%
  ungroup() %>%
  mutate(
    n = map_int(data, nrow),

    mcc = map(data, ~bootstrap_mcc(data = .x, R = 1000, method = 2)),
    accuracy = map_dbl(data, ~sum(.x$ground == .x$predicted) / nrow(.x)),
    accuracy_se = sqrt(accuracy * (1 - accuracy) / n)
  ) %>%
  unnest(mcc)

ggplot(patients_mcc , aes(x = method, y = mcc_mean, color = method)) +
  geom_pointrange(aes(ymin = mcc_ci_lo, ymax = mcc_ci_hi),
                  position = position_dodge2(width = 0.25)) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, by = 0.25)) +
  theme_bw() +
  labs(x = "Method", y = "MCC")

ggplot(patients_mcc , aes(x = method, y = accuracy, color = method)) +
  geom_pointrange(aes(ymin = accuracy - accuracy_se,
                      ymax = accuracy + accuracy_se),
                  position = position_dodge2(width = 0.25)) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, by = 0.25)) +
  theme_bw() +
  labs(x = "Method", y = "Accuracy")

lesioned_pct <- read_csv("../data/lesioned_pct.csv", show_col_types = FALSE)

# Are larger lesions in certain ROIs more likely to lead to correct detection?
patients_roi <- patients_lda %>%
  select(sub, ground, predicted) %>%
  mutate(
    correct = ground == predicted
  ) %>%
  left_join(lesioned_pct) %>%
  group_by(ROI) %>%
  nest() %>%
  mutate(

    binomial = map(data, ~glm(correct ~ lesioned_pct, data = .x,
                              family = "binomial")),

    coef = map_dbl(binomial, ~coef(.x)[2]),
    p = map_dbl(binomial, ~summary(.x)$coefficients[2, 4]),

  ) %>%
  arrange(p)

# Hold this out bc FDR needs the values in the proper vector
patients_roi$pcor <- p.adjust(patients_roi$p, method = "fdr")

# LOOCV ====

results_loocv <- results0 %>%
  filter(
    str_detect(f, "loocv")
  ) %>%
  unnest(data) %>%
  select(-`...1`) %>%
  mutate(
    outcome = str_extract(f, "outcome-[^_]*"),
    held_out = str_extract(f, "group-[^_]*"),
    method = str_extract(f, "method-[^_]*") %>%
      str_remove("method-")
  ) %>%
  select(-f) %>%
  arrange(sub, method, ground)

results_loocv_mcc <- results %>%
  group_by(method) %>%
  nest() %>%
  ungroup() %>%
  mutate(
    n = map_int(data, nrow),

    mcc = map(data, ~bootstrap_mcc(data = .x, R = 1000, method = 2)),
    accuracy = map_dbl(data, ~sum(.x$ground == .x$predicted) / nrow(.x)),
    accuracy_se = sqrt(accuracy * (1 - accuracy) / n),

    held_out = "loocv"
  )

results_loocv_mcc2 <- results_loocv_mcc %>%
  select(-data) %>%
  unnest(mcc)

ggplot(results_loocv_mcc2,
        aes(x = method, y = mcc_mean, color = method)) +
  geom_pointrange(aes(ymin = mcc_ci_lo, ymax = mcc_ci_hi),
                  position = position_dodge2(width = 0.25)) +
  scale_y_continuous(limits = c(0, 1.1), breaks = seq(0, 1, by = 0.25)) +
  theme_bw() +
  labs(x = "Method", y = "MCC")

# Combine ====

all_mcc2 <- bind_rows(results_mcc2, results_loocv_mcc2)

ggplot(all_mcc2,
        aes(x = method, y = mcc_mean, color = held_out, shape = held_out)) +
  geom_pointrange(aes(ymin = mcc_ci_lo, ymax = mcc_ci_hi),
                  position = position_dodge2(width = 0.25),
                  size = 1) +
  geom_hline(yintercept = 0, color = "red") +
  scale_y_continuous(limits = c(0, NA), breaks = seq(0, 1, by = 0.25)) +
  scale_shape_manual(values = c(15, 16, 17, 18, 20)) +
  theme_bw() +
  labs(x = "Method", y = "MCC", title = "MCC")

ggplot(all_mcc2,
        aes(x = method, y = accuracy, color = held_out, shape = held_out)) +
  geom_pointrange(aes(ymin = accuracy - accuracy_se,
                      ymax = accuracy + accuracy_se),
                  position = position_dodge2(width = 0.25),
                  size = 1) +
  geom_hline(yintercept = 0.5, color = "red") +
  scale_y_continuous(limits = c(0.5, NA), breaks = seq(0, 1, by = 0.1)) +
  scale_shape_manual(values = c(15, 16, 17, 18, 20)) +
  theme_bw() +
  labs(x = "Method", y = "Accuracy", title = "Accuracy")

# Everybody ====

everybody_files <- list.files("results/", "*.csv", full.names = TRUE,
                           recursive = TRUE) %>%
  str_subset("everybody")

everybody_results0 <- tibble(f = result_files) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE)
  )

everybody_results <- everybody_results0 %>%
  unnest(data) %>%
  select(-`...1`) %>%
  mutate(
    method = str_extract(f, "method-[^_]*") %>%
      str_remove("method-")
  ) %>%
  select(-f) %>%
  arrange(sub, method, ground) %>%
  separate_wider_delim(c(ground, predicted), delim = "_",
                       names = c("patient", "hemi"), names_sep = "_",
                       too_few = "align_start",
                       cols_remove = FALSE) %>%
  rename(
    ground = ground_ground,
    predicted = predicted_predicted
  ) %>%
  filter_out(
    ground %in% c("L", "R")
  )

everybody_mcc <- everybody_results %>%
  group_by(method) %>%
  nest() %>%
  mutate(

    n = map_int(data, nrow),

    mcc_overall = map(data, ~bootstrap_mcc(data = .x, R = 1000, method = 2)),
    accuracy_overall = map_dbl(data, ~sum(.x$ground == .x$predicted) / n),

    mcc_hemi = map(data, ~bootstrap_mcc(data = .x, R = 1000,
                                        ground = "ground_hemi",
                                        predicted = "predicted_hemi",
                                        method = 2)),
    accuracy_hemi = map_dbl(data,
                            ~sum(.x$ground_hemi == .x$predicted_hemi) / n),


    mcc_patient = map(data, ~bootstrap_mcc(data = .x, R = 1000,
                                        ground = "ground_patient",
                                        predicted = "predicted_patient",
                                        method = 2)),
    accuracy_patient = map_dbl(data,
                               ~sum(.x$ground_patient == .x$predicted_patient) / n),

    # accuracy__overall_se = sqrt(accuracy * (1 - accuracy) / n),
    across(matches("^accuracy_"), ~.x + sqrt(.x * (1 - .x) / n),
           .names = "{.col}_ci_hi"),

    across(matches("^accuracy_[^_]*$"), ~.x - sqrt(.x * (1 - .x) / n),
           .names = "{.col}_ci_lo")

  ) %>%
  unnest(starts_with("mcc"), names_sep = "_")

everybody_mcc2 <- everybody_mcc %>%
  select(-data, -n) %>%
  pivot_longer(-method) %>%
  separate_wider_delim(name, delim = "_", names = c("metric", "ground", "name"),
                       too_many = "merge", too_few = "align_start") %>%
  mutate(
    name = str_remove(name, "mcc_") %>%
      replace_na("mean")
  ) %>%
  pivot_wider()

red_lines <- tribble(
  ~metric, ~ground, ~random,
  "mcc", "hemi", 0,
  "mcc", "overall", 0,
  "mcc", "patient", 0,
  "accuracy", "hemi", 0.5,
  "accuracy", "overall", 1 / 3,
  "accuracy", "patient", 0.5,
)

ggplot(everybody_mcc2, aes(x = method, color = ground)) +
  geom_pointrange(aes(ymin = ci_lo, y = mean, ymax = ci_hi),
                  position = position_dodge(width = 0.25)) +
  geom_hline(
    data = red_lines,
    aes(yintercept = random),
    color = "red"
  ) +
  facet_grid(cols = vars(metric), rows = vars(ground)) +
  theme_bw()

ggplot(filter(everybody_mcc2, metric == "mcc"), aes(x = method, color = ground)) +
  geom_pointrange(aes(ymin = ci_lo, y = mean, ymax = ci_hi),
                  position = position_dodge(width = 0.25)) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_bw()

