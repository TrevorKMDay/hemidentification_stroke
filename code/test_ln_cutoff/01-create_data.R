library(tidyverse)

setwd("~/MyDrive/Projects/hemisphere_fingerprinting_stroke/code/test_ln_cutoff/")

# These values are cu mm
ln_sizes_cm <- c(L = 185912, R = 183872) * 0.001

max_size_cm <- 100 # DeMarco et al. (2022)
ln_maxpct <- 100 / ln_sizes_cm

# Load data

data1 <- read_csv("../input/group-all_parc-brainnetome_normed.csv", 
                    show_col_types = FALSE)

data <- data1 %>%
    select(sub, date, hemi, PatientStatus, status, ends_with("pctlesioned"),
            starts_with("cx_")) %>%
    filter(
        PatientStatus == "LHS"
    )

# Process

# Weird floating-point errors, just round
cutoffs <- round(seq(0, ln_maxpct["L"], by = 0.05), 2)
dir.create("inputs", showWarnings = FALSE)

for (cutoff in cutoffs) {

    tdata <- data %>%
        mutate(
            is_lhln_lesioned = if_else(lh_ln_pctlesioned > cutoff, "lesioned",
                                        "not_lesioned"),
            .after = PatientStatus
        )

    write_csv(tdata, 
        str_glue("inputs/input-lhs_lesion-{sprintf('%.2f', cutoff)}.csv"))
    
    reticulate::py_save_object(tdata,
        str_glue("inputs/input-lhs_lesion-{sprintf('%.2f', cutoff)}.pickle"))
    
}
