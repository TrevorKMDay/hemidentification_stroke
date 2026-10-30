library(tidyverse)
library(reticulate)
library(infer)

library(ggseg)
library(ggsegBrainnetome)

setwd("~/MyDrive/Projects/hemisphere_fingerprinting_stroke/code")

# Get the true model =====

true_model <- py_load_object("results/input-all_outcome-class_model.pickle")
true_features <- true_model[[2]]$feature_names_in_
true_scalings <- as_tibble(true_model[[2]]$scalings_)

# Create nice tible of true results
scalings <- true_scalings %>%
  add_column(feature = true_features, .before = 1) %>%
  select_all(~str_replace(., "V", "LD"))

# What p-threshold are we using for initial feature selection?
p_thresh <- 0.05

if (!file.exists("p_simulation.rds")) {

  # Get the bootstrapped scalings
  n <- 10000
  bs0 <- read_csv(str_glue("bs_results/group-all_outcome-class_bs-{n}.zip"),
                  show_col_types=FALSE, num_threads=10)

  bs <- bs0 %>%
    pivot_longer(starts_with("LD"), names_to = "LD") %>%
    group_by(feature, LD) %>%
    nest() 

  ld1 <- bs %>%
    filter(
      LD == "LD1"
    ) %>%
    left_join(
      select(scalings, feature, LD1),
      join_by(feature)
    ) %>%
    mutate(
      
      p = map2_dbl(data, LD1,
                    ~get_p_value(tibble(scaling1 = .x$value), .y,
                                direction = "two_sided")[[1]],
                  .progress = TRUE),

    )

  table(ld1$p < p_thresh)

  ld2 <- bs %>%
    filter(
      LD == "LD2"
    ) %>%
    left_join(
      select(scalings, feature, LD2),
      by = join_by(feature)
    ) %>%
    mutate(
      
      p = map2_dbl(data, LD2,
                    ~get_p_value(tibble(scaling1 = .x$value), .y,
                                direction = "two_sided")[[1]],
                  .progress = TRUE),

    )

  table(ld2$p < p_thresh)
  table(ld1$p < p_thresh, ld2$p < p_thresh)

  save_me <- bind_rows(ld1, ld2) %>%
    select(-data)

  # Don't save out the null distributions to save space
  write_rds(save_me, "p_simulation.rds")
  
} else {

  x <- read_rds("p_simulation.rds")

  ld1 <- filter(x, LD == "LD1")
  ld2 <- filter(x, LD == "LD2")

}

# By region =====

ld1_by_region <- ld1 %>%
  ungroup() %>%
  select(feature, p) %>%
  mutate(
    sig = p < p_thresh
  ) %>%
  separate_wider_delim(feature, delim = "_", names = c(NA, "ROI1", "ROI2")) %>%
  pivot_longer(c(ROI1, ROI2), values_to = "ROI") %>%
  group_by(ROI) %>%
  summarize(
    sig = sum(sig)
  ) %>%
  arrange(desc(sig))


ld2_by_region <- ld2 %>%
  ungroup() %>%
  select(feature, p) %>%
  mutate(
    sig = p < p_thresh
  ) %>%
  separate_wider_delim(feature, delim = "_", names = c(NA, "ROI1", "ROI2")) %>%
  pivot_longer(c(ROI1, ROI2), values_to = "ROI") %>%
  group_by(ROI) %>%
  summarize(
    sig = sum(sig)
  ) %>%
  arrange(desc(sig))

# My data removes special characters from the region names, but ggseg still
# has them. Create a table to replace the cleaned names with their actual names.
fix_names <- brainnetome()$core %>%
  select(region) %>%
  mutate(
    region = str_remove(region, "_[LR]$"),
    clean_name = str_remove_all(region, "[/+.]")
  )

by_region <- left_join(ld1_by_region, ld2_by_region, by = join_by(ROI),
                        suffix = c(".ld1", ".ld2")) %>%
  mutate(

    # Recode values to match Brainnetome atlas
    region = recode_values(ROI, from = fix_names$clean_name, 
                           to = fix_names$region) %>%
      paste0(., "_L"),

    # Don't plot 0 values
    across(starts_with("sig"), ~replace_values(.x, 0 ~ NA))

  )
    

cor(by_region$sig.ld1, by_region$sig.ld2, use = "c")

ggplot(by_region, aes(x = sig.ld1, y = sig.ld2)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_bw()

by_region_wide <- by_region %>%
  select(-ROI) %>%
  pivot_longer(starts_with("sig"), names_to = "LD", names_pattern = "sig.(.*)")

breaks <- c(1, 5, 10, 15, 20)

ggplot(filter(by_region_wide, LD == "ld1")) +
  geom_brain(
    atlas = brainnetome(),
    hemi = "left",
    aes(fill = sqrt(value))
  ) +
  viridis::scale_fill_viridis(limits = c(NA, NA), 
                              breaks = sqrt(breaks),
                              labels = breaks) +
  facet_grid(rows = vars(LD)) +
  theme_bw() +
  labs(title = "Regions with >0 significant connections (LD1)",
        fill = "N", 
        caption = "Fills sqrt(val)") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

ggplot(filter(by_region_wide, LD == "ld2")) +
  geom_brain(
    atlas = brainnetome(),
    hemi = "left",
    aes(fill = sqrt(value))
  ) +
  viridis::scale_fill_viridis(limits = c(NA, NA), 
                              breaks = sqrt(breaks),
                              labels = breaks) +
  facet_grid(rows = vars(LD)) +
  theme_bw() +
  labs(title = "Regions with >0 significant connections (LD2)",
        fill = "N", 
        caption = "Fills sqrt(val)") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

# Bootstrap connections 

# Number and probability of sig connection for binom
n_sig <- c("ld1" = sum(ld1$p < p_thresh), "ld2" = sum(ld2$p < p_thresh))
p_sig <- n_sig / nrow(ld1)

# how to do a binomial test:
# binom.test(21, nrow(by_region), p = p_sig["ld1"], alternative = "g")

by_region2 <- by_region %>%
  mutate(

    across(starts_with("sig."), ~replace_na(., 0)),

    sig.ld1.binom = map(sig.ld1, 
                        ~binom.test(x = .x, n = nrow(by_region), 
                                    p = p_sig["ld1"], alternative = "g")),
    sig.ld2.binom = map(sig.ld2, 
                        ~binom.test(x = .x, n = nrow(by_region), 
                                    p = p_sig["ld2"], alternative = "g")),

    sig.ld1.p = map_dbl(sig.ld1.binom, ~ .x$p.value),
    sig.ld2.p = map_dbl(sig.ld2.binom, ~ .x$p.value),

    sig.ld1_thresh = if_else(sig.ld1.p < .1, sig.ld1, NA),
    sig.ld2_thresh = if_else(sig.ld2.p < .1, sig.ld2, NA),

  )

# Merge together to put in one plot
by_region3 <- by_region2 %>%
  select(region, ends_with("thresh")) %>%
  pivot_longer(ends_with("thresh"), names_to = c("LD"), 
                names_pattern = "sig.(.*)_thresh")

ggplot(filter(by_region3, LD == "ld1")) +
  geom_brain(
    atlas = brainnetome(),
    hemi = "left",
    aes(fill = sqrt(value))
  ) +
  viridis::scale_fill_viridis(limits = c(NA, NA), 
                              breaks = sqrt(breaks),
                              labels = breaks) +
  facet_grid(rows = vars(LD)) +
  theme_bw() +
  labs(title = "Regions with more connections than expected by chance (LD1)",
        fill = "N") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())

ggplot(filter(by_region3, LD == "ld2")) +
  geom_brain(
    atlas = brainnetome(),
    hemi = "left",
    aes(fill = sqrt(value))
  ) +
  viridis::scale_fill_viridis(limits = c(NA, NA), 
                              breaks = sqrt(breaks),
                              labels = breaks) +
  facet_grid(rows = vars(LD)) +
  theme_bw() +
  labs(title = "Regions with more connections than expected by chance (LD2)",
        fill = "N") +
  theme(axis.text = element_blank(), axis.ticks = element_blank())


by_region2 %>%
  filter(
    !is.na(sig.ld1_thresh),
    !is.na(sig.ld2_thresh)
  )
  