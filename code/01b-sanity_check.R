library(tidyverse)

all_values <- read_rds("all_values_temp.rds")

# Compare patients/HC as a sanity check ====

by_status <- all_values %>%
  group_by(PatientStatus, ROI1, ROI2) %>%
  summarize(
    mean_r = mean(value, na.rm = TRUE)
  )

axis_order <- rois_andrew %>%
  filter_out(
    # Remove non-cortical ROIs
    Label %in% c(rois_cerebellum, rois_subcortical)
  ) %>%
  arrange(Hemisphere, Lobe, Gyrus, Label)

ggplot(by_status, aes(x = ROI1, y = ROI2, fill = mean_r)) +
  geom_tile() +
  scale_x_discrete(limits = axis_order$Label) +
  scale_y_discrete(limits = axis_order$Label) +
  coord_equal() +
  facet_wrap(vars(PatientStatus)) +
  theme_bw()

# Difference between groups
difference <- by_status %>%
  ungroup() %>%
  pivot_wider(names_from = PatientStatus, values_from = mean_r) %>%
  mutate(
    LHSmHC = LHS - HC
  )

ggplot(difference, aes(x = ROI1, y = ROI2, fill = LHSmHC)) +
  geom_tile() +
  scale_x_discrete(limits = axis_order$Label) +
  scale_y_discrete(limits = axis_order$Label) +
  scale_fill_gradient2() +
  theme_bw()
