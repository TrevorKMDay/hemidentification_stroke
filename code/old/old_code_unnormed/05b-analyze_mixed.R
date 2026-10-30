library(tidyverse)

setwd("~/MyDrive/Projects/hemisphere_fingerprinting_stroke/code/")

mixed_files <- list.files("predicted/", "Lpct", full.names = TRUE)

data0 <- tibble(f = mixed_files) %>%
  mutate(
    data = map(f, read_csv, show_col_types = FALSE)
  )

data <- data0 %>%
  unnest(data) %>%
  mutate(
    model = str_extract(f, "model-[^_]*"),
    Lpct = str_extract(f, "Lpct-0..") %>%
      str_remove("Lpct-") %>%
      as.numeric()
  ) %>%
  select(model, Lpct, sub, predicted, LD1, LD2)

data_hemi <- data %>%
  filter(
    model == "model-hemi"
  ) %>%
  mutate(
    Lpct_jtd = jitter(Lpct)
  )

ggplot(data_hemi, aes(x = Lpct_jtd, y = LD1)) +
  geom_point(aes(color = predicted)) +
  geom_line(aes(group = sub), alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw()

data_hemi_summary <- data_hemi %>%
  summarize(
    .by = Lpct,
    n = n(),
    predR = sum(predicted == "R")
  ) %>%
  mutate(
    pct_predR = predR / n
  )

ggplot(data_hemi_summary, aes(x = Lpct, y = pct_predR)) +
  geom_abline(slope = -1, intercept = 1, color = "red") +
  geom_line() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(limits = c(0, 1)) +
  theme_bw()

data_class <- data %>%
  filter(
    model == "model-class"
  ) %>%
  mutate(
    Lpct_jtd = jitter(Lpct)
  ) %>%
  pivot_longer(c(LD1, LD2))

ggplot(data_class, aes(x = Lpct_jtd, y = value)) +
  geom_point(aes(color = predicted)) +
  geom_line(aes(group = sub), alpha = 0.5) +
  scale_x_continuous(limits = c(0, 1)) +
  facet_wrap(vars(name), scales = "free") +
  theme_bw()


data_hemi %>%
  filter(
    Lpct == 0.1
  ) %>%
  View()

data_hemi %>%
  filter(
    sub == "sub-HYX"
  )
