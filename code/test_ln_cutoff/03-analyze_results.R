library(tidyverse)
library(viridis)

setwd("~/MyDrive/Projects/hemisphere_fingerprinting_stroke/code/test_ln_cutoff/")

source("../acc_se.R")
source("../../../hemisphere_fingerprinting/code/bootstrap_mcc.R")

ln_sizes_cm <- c(L = 185912, R = 183872) * 0.001

# Parallel setup ====

library(parallel)
library(doParallel)
library(foreach)

cores <- detectCores()
cl <- makeCluster(cores[1] - 1) 
registerDoParallel(cl)

####

result_files <- list.files("results/", full.names = TRUE)

data0 <- tibble(f = result_files) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE)
  )

data <- data0 %>%
  mutate(
    cutoff = str_extract(f, "lesion-0..."),
  ) %>%
  select(-f) %>%
  mutate(

    n = map_int(data, nrow),
    unlesioned = map_int(data, ~sum(.x$is_lhln_lesioned == "lesioned")),
    pct_unlesioned = unlesioned / n,
    
    acc = map(data, 
                ~acc_se(sum(.x$is_lhln_lesioned == .x$is_lhln_lesioned_predicted),
                        n = nrow(.x))
            ),

    mcc = map(data, bootstrap_mcc, ground = "is_lhln_lesioned", 
                predicted = "is_lhln_lesioned_predicted", .progress = TRUE)

  )

acc <- data %>%
  select(-data, -mcc) %>%
  unnest(acc) %>%
  mutate(
    cutoff = as.numeric(str_remove(cutoff, "lesion-")),
    name = rep(c("acc", "acc_se", "ci1", "ci2"), length.out = nrow(.))
  ) %>%
  filter(
    name != "acc_se"
  ) %>%
  pivot_wider(values_from = acc)

mcc <- data %>%
  select(-data, -acc) %>%
  unnest(mcc) %>%
  mutate(
    cutoff = as.numeric(str_remove(cutoff, "lesion-")),
  )

ggplot(acc, aes(x = cutoff, y = acc)) +
  geom_pointrange(aes(ymin = ci1, ymax = ci2)) +
  geom_line() +
  geom_label(aes(label = round(pct_unlesioned, 2))) +
  geom_hline(yintercept = 0.5, color = "red") +
  scale_y_continuous(
    limits = c(NA, 1),                
  ) +
  scale_x_continuous(
    sec.axis = sec_axis(~ .x * ln_sizes_cm["L"], name = expression(cm^3))
  ) +
  theme_bw() +
  labs(x = "% of LN", y = "Accuracy (95% CI)")

ggplot(mcc, aes(x = cutoff, y = mcc_mean)) +
  geom_pointrange(aes(ymin = mcc_ci_lo, ymax = mcc_ci_hi)) +
  geom_line() +
  geom_label(aes(label = round(pct_unlesioned, 2))) +
  geom_hline(yintercept = 0, color = "red") +
  scale_y_continuous(limits = c(NA, 1)) +
  scale_x_continuous(
    sec.axis = sec_axis(~ .x * ln_sizes_cm["L"], name = expression(cm^3))
  ) +
  theme_bw() +
  labs(x = "% of LN", y = "MCC (95% CI)")

both <- left_join(acc, mcc) %>%
  rename(
    acc_mean = acc, 
    acc_lo = ci1,
    acc_hi = ci2,
    mcc_lo = mcc_ci_lo,
    mcc_hi = mcc_ci_hi
  ) %>%
  pivot_longer(names_to = c("metric", "ci"), names_sep = "_",
                matches("acc|mcc")) %>%
  pivot_wider(names_from = ci)


ggplot(both, aes(x = cutoff, y = mean, color = metric)) +
  geom_hline(yintercept = c(0, 0.5), color = c("blue", "red")) +
  geom_pointrange(aes(ymin = lo, ymax = hi)) +
  geom_line() +
  geom_label(aes(label = round(pct_unlesioned, 2))) +
  scale_y_continuous(limits = c(NA, 1)) +
  scale_x_continuous(
    sec.axis = sec_axis(~ .x * ln_sizes_cm["L"], name = expression(cm^3))
  ) +
  theme_bw() +
  labs(x = "% of LN", y = "Metric (95% CI)")

# confusion

data2 <- data %>%
  mutate(
    confusion = map(data, 
            ~group_by(.x, is_lhln_lesioned, is_lhln_lesioned_predicted) %>%
                summarize(n2 = n())
        )
  ) %>%
  select(cutoff, n, confusion) %>%
  unnest(confusion) %>%
  mutate(
    pct = n2 / n
  )
 
ggplot(data2, 
    aes(y = is_lhln_lesioned, x = is_lhln_lesioned_predicted, fill = pct)) +
  geom_tile() +
  geom_text(aes(label = str_glue("{round(100 * pct)}%"))) +
  scale_x_discrete(position = "top") +
  scale_fill_viridis(limits = c(0, 1)) +
  facet_wrap(vars(cutoff)) +
  theme_bw() +
  labs(x = "Predicted", y = "Ground")
