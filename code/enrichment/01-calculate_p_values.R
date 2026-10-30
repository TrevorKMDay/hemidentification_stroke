library(tidyverse)
library(infer)

# Setup ====

setwd("~/MyDrive/Projects/hemisphere_fingerprinting_stroke/code/enrichment/")
select <- dplyr::select

source("00-enrichment_functions.R")

## Read in bootstrapping ====

n <- 5000
bs_rds <- str_glue("bs_n-{n}.rds")

if (!file.exists(bs_rds)) {

  # Load the data from the ZIP file, which is not optimally compressed and 
  # move it to an R compressed file, which is more compressed and faster for
  # R to read.

  bs_file <- str_glue("bs_results/data-hclhs.null_shuffle-class_n-{n}_bsscalings.zip")

  bs <- read_csv(bs_file, show_col_types = FALSE) %>%
    select(feature, bs, starts_with("LD"))

  qs2::qs_save(bs, bs_rds, compress_level = 9)
    
} else {

  bs <- qs2::qs_read(bs_rds, nthreads = 10)

}

bs_rois <- bs %>%
  select(feature) %>%
  distinct() %>%
  separate_wider_delim(feature, "_", names = c(NA, "ROI1", "ROI2")) 

sort(unique(c(bs_rois$ROI1, bs_rois$ROI2)))

# Results ====

true_model_obj <- reticulate::py_load_object("../results/input-all_outcome-class_model.pickle")

features <- true_model_obj[[2]]$feature_names_in_

if (!(length(unique(features)) == length(unique(bs$feature)))) {

  message(str_glue("There is a mismatch between the number of features in the", 
                   "null model and the true model!"))

}

true_model <- as_tibble(true_model_obj[[2]]$scalings_) %>%
  add_column(feature = features, .before = 1) %>%
  select_all(~str_replace(., "V", "scaling"))

## Hemisphere axis ====

p_values <- bs %>%
  group_by(feature) %>%
  nest() %>%
  left_join(true_model, join_by(feature)) %>%
  mutate(

    p_hemi = map2_dbl(data, scaling1,
                      ~get_p_value(tibble(scaling1 = .x$LD1), .y,
                                    direction = "two_sided")[[1]],
                      .progress = TRUE),

    p_stroke = map2_dbl(data, scaling2,
                        ~get_p_value(tibble(scaling2 = .x$LD2), .y,
                                      direction = "two_sided")[[1]],
                        .progress = TRUE),

  )

# How many had simulated p-values of 0:
table(p_values$p_hemi == 0, p_values$p_stroke == 0)

p_values2 <- p_values %>%
  mutate(

    # Replace 0 values with the lowest possible p-values, otherwise actually
    # do the transform
    across(starts_with("p_"), 
            ~recode_values(.x, 0 ~ log10(3 / n), default = log10(.x)),
            .names = "log_{.col}")

  ) %>%
  select(-data)

write_rds(p_values2, "hemi.stroke_ps.rds")

p_values3 <- p_values2 %>%
  select(feature, starts_with("log_")) %>%
  pivot_longer(-feature)

ggplot(p_values3, aes(x = -value)) +
  geom_histogram(boundary = 0) +
  facet_grid(rows = vars(name)) +
  theme_bw()

  