library(tidyverse)
library(reticulate)

all <- read_csv("input/group-HC_parc-brainnetome_loocv_rZ.csv",
                show_col_types = FALSE)

all_long <- all %>%
  select(-class) %>%
  pivot_longer(-c(sub, hemi, PatientStatus, group)) %>%
  pivot_wider(names_from = hemi)

mix_LR <- function(data, ratio) {

  x <- data %>%
    mutate(
      avg = map2_dbl(L, R, ~weighted.mean(c(.x, .y), ratio))
    )

  y <- x %>%
    select(-L, -R) %>%
    pivot_wider(values_from = avg) %>%
    mutate(
      # Model fitting code expects some value
      hemi = "avg",
      class = "LHS_avg"
    )

  return(y)

}

L_amt <- seq(0.1, 0.9, by = 0.1)

for (L in L_amt) {

  avg <- mix_LR(all_long, c(L, 1 - L))

  write_csv(avg, str_glue("input/mixed/group-HC_parc-brainnetome_",
                          "Lpct-{L}_rZ.csv"))

  py_save_object(avg, str_glue("input/mixed/group-HC_parc-brainnetome_",
                               "Lpct-{L}_rz.pickle"))

  message(L)

}
