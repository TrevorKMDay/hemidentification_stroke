library(tidyverse)

setwd("~/MyDrive/Projects/hemisphere_fingerprinting_stroke/code")

holdout_result_files <- list.files("predicted/holdouts/", pattern = "*heldout*",
                                    full.names = TRUE)

holdout_results_0 <- tibble(f = holdout_result_files) %>%
  mutate(
    
    data = map(f, read_csv, show_col_types = FALSE),
    heldout = str_extract(f, "heldout-...") %>%
      str_replace("heldout", "sub"),
    outcome = str_extract(f, "outcome-[^_]*") %>%
      str_remove("outcome-")

  ) 

holdout_results <- holdout_results_0 %>%
  select(heldout, data, outcome) %>%
  unnest(data) %>%
  select(heldout, sub, PatientStatus, hemi, outcome, date, mri_chronicity_mo, 
          WAB_value, LD1) %>%
  group_by(heldout, sub, PatientStatus) %>%
  mutate(
    days = interval(min(date), date) / days(1),
    sub_heldout = heldout == sub
  )

multiple_visits <- holdout_results %>%
  group_by(heldout, sub, PatientStatus) %>%
  filter(
    length(unique(date)) > 1
  ) %>%
  pull(sub) %>%
  unique()

holdout_results_onlylongitudinal <- filter(holdout_results, sub %in% multiple_visits)

ggplot(filter(holdout_results_onlylongitudinal, outcome == "hemi"), 
        aes(x = days, y = LD1)) +
  geom_point(aes(group = sub, color = sub_heldout), alpha = 0.5) +
  geom_line(aes(group = interaction(sub, hemi), color = sub_heldout), 
            alpha = 0.5) +
  scale_color_manual(values = c("black", "red")) +
  facet_wrap(vars(heldout)) +
  theme_bw() +
  labs(x = "Days after first scan", y = "Hemisphere axis score (LD1)",
        title = "All participants, pred=hemi")

ggplot(filter(holdout_results_onlylongitudinal, outcome == "class"), 
        aes(x = days, y = LD1)) +
  geom_point(aes(group = sub, color = sub_heldout), alpha = 0.5) +
  geom_line(aes(group = interaction(sub, hemi), color = sub_heldout), 
            alpha = 0.5) +
  scale_color_manual(values = c("black", "red")) +
  facet_wrap(vars(heldout)) +
  theme_bw() +
  labs(x = "Days after first scan", y = "Hemisphere axis score (LD1)",
        title = "All participants, pred=class")

ggplot(filter(holdout_results_onlylongitudinal, outcome == "class"), 
        aes(x = days, y = LD1,  color = sub_heldout)) +
  geom_point(aes(group = sub), alpha = 0.5) +
  geom_line(aes(group = interaction(heldout, hemi)), alpha = 0.5) +
  scale_color_manual(values = c("black", "red")) +
  facet_wrap(vars(sub), scales = "free_x") +
  theme_bw() +
  theme(legend.position = "bottom")+ 
  labs(x = "Days after first scan", y = "Hemisphere axis score (LD1)",
        title = "All participants per holdout group, pred=class")

ggplot(filter(holdout_results_onlylongitudinal, outcome == "class"), 
        aes(x = days, y = LD1, color = sub)) +
  geom_point(alpha = 0.5) +
  geom_line(aes(group = interaction(sub, hemi)), alpha = 0.5) +
  geom_smooth(aes(group = hemi), method = "lm") +
  facet_wrap(vars(heldout)) +
  theme_bw() +
  labs(x = "Days after first scan", y = "Hemisphere axis score (LD1)",
        title = "All participants per holdout group, pred=class")

# Variance 

variance <- holdout_results_onlylongitudinal %>%
  filter(
    heldout == sub
  ) %>%
  select(-heldout) %>%
  group_by(sub) %>%
  mutate(
    max_interval = interval(min(date), max(date)) / days(1)
  ) %>%
  group_by(outcome, sub, PatientStatus, hemi, max_interval) %>%
  summarize(
    ld1_sd = sd(LD1),
  )

mean_hc_variance <- variance %>%
  filter(
    PatientStatus == "HC"
  ) %>%
  group_by(outcome, PatientStatus, hemi) %>%
  summarize(
    mean_sd = mean(ld1_sd),
    sd_sd = sd(ld1_sd)
  ) %>%
  filter(
    hemi == "R"
  )

library(ggrepel)

ggplot(variance) +
  geom_point(aes(x = max_interval, y = ld1_sd, color = PatientStatus), 
              size = 3) +
  geom_label_repel(aes(label = str_glue("{sub} {hemi}"), x = max_interval, 
                    y = ld1_sd, fill = PatientStatus), 
                   size = 3) +
  geom_line(aes(x = max_interval, y = ld1_sd,  group = sub, color = PatientStatus)) +
  geom_hline(
    data = mean_hc_variance,
    aes(yintercept = mean_sd - sd_sd)
  ) +
  geom_hline(
    data = mean_hc_variance,
    aes(yintercept = mean_sd + sd_sd)
  ) +
  facet_wrap(vars(outcome)) +
  theme_bw()

# All rows 

allrows_hemi <- read_csv("predicted/input-all_outcome-hemi_predictedon-allrows_results.csv")

allrows_class <- read_csv("predicted/input-all_outcome-class_predictedon-allrows_results.csv") %>%
  arrange(sub, date)

allrows_class %>%
  group_by(sub, PatientStatus) %>%
  summarize(
    n = length(unique(date))
  ) %>%
  arrange(desc(n)) %>%
  filter(
    PatientStatus == "HC"
  )

ggplot(allrows_hemi, aes(x = LD1, y = sub, color = PatientStatus)) +
  geom_point()

ggplot(allrows_class, aes(x = LD1, y = LD2)) +
  geom_path(aes(group = interaction(sub, hemi))) +
  geom_point(aes(shape = interaction(PatientStatus, hemi), 
                 color = as.factor(scan_no), fill = PatientStatus), 
             size = 2, shape = 21) +
  scale_x_continuous(limits = c(NA, 3)) +
  scale_y_continuous(limits = c(NA, 3)) +
  theme_bw()

long_hc <- allrows_class %>%
  filter(
    PatientStatus == "HC",
    sub %in% multiple_visits
  )

# Variance 2 

variance2 <- allrows_class %>%
  group_by(sub) %>%
  mutate(
    max_interval = interval(min(date), max(date)) / days(1)
  ) %>%
  group_by(sub, PatientStatus, hemi, max_interval) %>%
  summarize(
    ld1_sd = sd(LD1),
  ) %>%
  filter_out(
    max_interval == 0
  )

ggplot(variance2) +
  geom_line(aes(x = max_interval, y = ld1_sd, group = sub)) +
  geom_point(aes(x = max_interval, y = ld1_sd, 
                  color = interaction(PatientStatus, hemi)), 
              size = 3) +
  geom_label_repel(aes(label = str_glue("{sub} {hemi}"), x = max_interval, 
                    y = ld1_sd, fill = interaction(PatientStatus, hemi)), 
                   size = 3) +
  theme_bw()

