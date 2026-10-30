library(tidyverse)

# Overall results ====

setwd("~/MyDrive/Projects/hemisphere_fingerprinting_stroke/code/dim_reduce/")

results0 <- read_csv("group-HC_method-LOOCV_results.csv") %>%
  select(-`...1`)

results <- results0 %>%
  group_by(k) %>%
  nest() %>%
  mutate(
    n = map_int(data, nrow),
    total_correct = map_int(data, ~sum(.x$hemi == .x$predicted)),
    p = total_correct / n,
    se = sqrt((p * (1 - p)) / n)
  )

ggplot(results, aes(x = k, y = p)) +
  geom_pointrange(aes(ymin = p - se, ymax = p + se)) +
  geom_line()

# Selector scores ====

selectorfiles <- list.files(pattern = ".*scores.csv")

sel0 <- tibble(f = selectorfiles) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE),
  )

sel <- sel0 %>%
  unnest(data) %>%
  group_by(feature) %>%
  summarize(
    mean = mean(score),
    sd = sd(score)
  ) %>%
  arrange(desc(mean)) %>%
  mutate(
    rank = row_number()
  )

ggplot(sel, aes(x = rank, y = mean)) +
  geom_point()
