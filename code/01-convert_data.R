library(tidyverse)
library(reticulate)
library(groupdata2)
library(googlesheets4)

setwd("~/MyDrive/Projects/hemisphere_fingerprinting_stroke/code/")

# Setup ====

# Keep these in a file to not put on GH
rename_subs <- read_csv("relabel_subs.csv", col_names = c("from", "to"),
                        show_col_types = FALSE)

## ROI names

rois <- readxl::read_xlsx("../data/language_network_ids.xlsx")
rois_andrew <- readxl::read_excel("../data/BN_Atlas_274_combined.xlsx")

rois_cerebellum <- rois_andrew$Label[rois_andrew$Lobe == "Cerebellum"]
rois_subcortical <- rois_andrew$Label[rois_andrew$Lobe == "Subcortical"]

rois_cortical <- rois_andrew %>%
  filter_out(
    Label %in% rois_cerebellum |
    Label %in% rois_subcortical
  )

## Get lesions

roi_labels <- str_remove_all(rois_andrew$Label, "[/.+]")

lesion0 <- read_csv("datasets/lesion_pct_by_division.csv", 
                    show_col_types = FALSE)

lesion_wide <- lesion0 %>%
  select(sub, date, LR, entire_corticalhemilesion_mm3) %>%
  distinct() %>%
  pivot_wider(id_cols = c(sub, date), names_from = LR, 
              values_from = entire_corticalhemilesion_mm3) %>%
  select(-V)

is_ln_lesioned <- lesion0 %>%
  select(sub, date, LR, division, starts_with("total")) %>%
  filter(
    division %in% c("Language", "Semantic")
  ) %>%
  group_by(sub, date, LR) %>%
  summarize(
    total_lesion = sum(total_lesion_mm3),
    total_volume = sum(total_volume_mm3)
  ) %>%
  mutate(
    pct = total_lesion / total_volume 
  ) 

ggplot(is_ln_lesioned, aes(pct)) +
  geom_histogram(boundary = 0)

# Only consider the LN lesioned if it's at least x% lesioned
lesion_threshold <- 0.05

is_ln_lesioned2 <- is_ln_lesioned %>%
  pivot_wider(id_cols = c(sub, date), names_from = LR, values_from = pct) %>%
  mutate(
    is_ln_lesioned = case_when(
      L > lesion_threshold & R > lesion_threshold ~ "both",
      L > lesion_threshold & R <= lesion_threshold ~ "LH",
      L <= lesion_threshold & R > lesion_threshold ~ "RH",
      L <= lesion_threshold & R <= lesion_threshold ~ "none"
    )
  ) %>%
  rename(
    lh_ln_pctlesioned = L,
    rh_ln_pctlesioned = R
  )

table(is_ln_lesioned2$is_ln_lesioned)

## Get project labels

projects <- read_csv("../data/projectnames.csv", show_col_types = FALSE) %>%
  mutate(
    sub = str_glue("sub-{SubjectCode}")
  ) %>%
  select(sub, ScanDate, projectname) %>%
  distinct() %>%
  filter_out(
    is.na(projectname)
  )

table(projects$projectname, useNA = "a")

## Build demographic table ====

demo1 <- read_csv("datasets/demographics.csv", show_col_types = FALSE)

no_sex_subs <- demo1$sub[is.na(demo1$Sex)]

## List of available data ====

data_files <- list.files("../data/individual_matrices_andrew/",
                            pattern = "sub-.*.csv", full.names = TRUE)

# Save this out 
scans1 <- tibble(f = data_files) %>%
  mutate(
    sub = str_extract(f, "sub-[^_]*"),
    date = str_remove(str_extract(f, "date-[0-9]*"), "date-") %>%
      as.Date(format = "%Y%m%d"),
    PatientStatus = case_when(
      str_detect(f, "patient-0") ~ "HC",
      str_detect(f, "patient-1") ~ "LHS"
    )
  ) %>%
  left_join(projects, by = join_by(sub, date == ScanDate))

# Exclusions list
gs4_auth("td758@georgetown.edu")
exclusions <- read_sheet("13cEBidZSOCNMWOI4JxfaS2mIxahjjLciv7xKKT3hH-A") %>%
  mutate(
    sub = paste0("sub-", sub)
  )

scans2 <- scans1 %>%
  left_join(exclusions, by = join_by(sub)) 

# Anyone labeled "exclude" with any note
excluded_subs_doc <- na.omit(scans2$sub[str_detect(scans2$status, "exclude")])

# Previously were excluding Anna's ppts, now we're gonna use them
annas_projs <- c("Vispat", "Greenwald KL2")

excluded_visits_anna <- scans2 %>%
  group_by(sub, PatientStatus) %>%
  summarize(
    all_annas = sum(projectname %in% annas_projs) == n()
  ) %>%
  filter(all_annas) %>%
  arrange(PatientStatus)

# These participants are excluded - keep track of thme
excluded_visits_anna %>%
  filter(
    PatientStatus == "LHS"
  ) %>%
  select(sub) %>%
  write_csv("lhs_only_in_Annas_projects.csv")

# Exclusions
excluded_subs <- sort(
    c(
      # Manual exclusions
      excluded_subs_doc, 
      # Participants who only ever did Anna's projects - assume they are RHS
      # (even though the label says 'LHS')
      excluded_visits_anna$sub[excluded_visits_anna$PatientStatus == "LHS"]
    )
  )

scans <- scans2 %>%
  filter_out(
    sub %in% excluded_subs
  ) %>%
  left_join(lesion_wide, by = join_by(sub, date)) %>%
  mutate(
    status = replace_na(status, "train_and_test"),
    sub = replace_values(sub, from = rename_subs$from, to = rename_subs$to)
  ) 

stroke_lr_ratio <- scans %>%
  mutate(
    stroke_ratio = (L - R) / (L + R)
  )

scans %>%
  group_by(PatientStatus, R > 0) %>%
  summarize(
    n = n()
  )

table(scans$projectname, useNA = "a")
length(unique(scans$sub))

# Number of ppts with each # of visits
table(table(scans$sub))

write_csv(scans, "participants_with_mri.csv")

demo <- left_join(demo1, exclusions) %>%
  mutate(
    status = replace_na(status, "train_and_test")
  )

scans_demo <- scans %>%
  left_join(demo, 
            by = join_by(sub, status, PatientStatus, date == ScanDate),
            relationship = "one-to-one") %>%
  mutate(

    PatientStatus = case_when(
      !is.na(PatientStatus) ~ PatientStatus,
      str_detect(f, "patient-1") ~ "LHS",
      str_detect(f, "patient-0") ~ "HC",
    )

  )

no_proj <- scans_demo %>%
  filter(
    is.na(projectname)
  )

no_status <- scans_demo %>%
  filter(
    is.na(PatientStatus)
  )

scans_final <- scans_demo %>%
  filter_out(
    is.na(projectname)
  ) %>%
  select(-f) %>%
  group_by(sub) %>%
  mutate(
    scan_no = row_number(), 
    best_scan = scan_no == max(scan_no),
    .after = date,
    status = if_else(best_scan, status, "test_long"),
    group = if_else(is.na(group), PatientStatus, group)
  )

table(scans_final$PatientStatus, scans_final$status, useNA = "i")

## ROIs =====

read_in_matrix <- function(f, cols) {

  mat <- read_csv(f, col_names = cols, show_col_types = FALSE) %>%
    mutate(
        across(everything(), ~replace_values(.x, Inf ~ NA))
    ) %>%
    add_column(ROI1 = cols, .before = 1)

  return(mat)

}

reformat <- function(mat) {

  mat2 <- mat %>%
    pivot_longer(-ROI1)

}

# Lesion locations ====

# TO DO: FIX THIS, IT HAS BEEN CALCULATED IN ANOTHER SCRIPT AND SHOULD RELY ON 
# THAT

# Load data ====

count_NaNs <- function(x) {

  xn <- select(x, where(is.numeric)) %>%
    as.matrix()

  nans <- sum(is.nan(xn))

  return(nans)

}

# TO DO: NO NEED TO DROP NO-PROJECT SUBS 

excluded_subs_regex <- str_remove_all(paste0(c(excluded_subs), collapse = "|"), 
                                      "sub-")

# Remove the excluded participants
data_files_to_read <- str_subset(data_files, excluded_subs_regex, 
                                  negate = TRUE) 
                              
# Read in connectivity data
data0 <- tibble(f = data_files_to_read) %>%
  mutate(
    cormat = map(f, read_in_matrix, cols = rois_andrew$Label, 
                  .progress = TRUE)
  )

data <- data0 %>%
  mutate(
    sub = str_extract(f, "sub-..."),
    date = str_extract(f, "date-[0-9]*") %>%
      str_remove("date-") %>%
      ymd(),
    long = map(cormat, ~pivot_longer(., -ROI1, names_to = "ROI2"))
  ) %>%
  select(sub, date, cormat, long)

nans1 <- data %>%
  mutate(
    nans = map_int(cormat, count_NaNs, .progress = TRUE)
  ) %>%
  arrange(desc(nans), sub)

sessions_with_nans <- filter(nans1, nans > 0) %>%
  select(sub, date)

# Extract all values, remove self-to-self
all_values1 <- data %>%
  anti_join(
    # Remove sessions with NaNs
    sessions_with_nans, 
    join_by(sub, date)
  ) %>%
  select(-cormat) %>%
  mutate(
    # fix labeling error
    sub = replace_values(sub, from = rename_subs$from, to = rename_subs$to)
  ) %>%
  left_join(
    scans_final,
    by = join_by(sub, date)
  ) %>%
  unnest(long) %>%
  filter_out(
    when_any(
      # Remove r=1
      ROI1 == ROI2,
      # Remove non-cortical ROIs
      ROI1 %in% c(rois_cerebellum, rois_subcortical),
      ROI2 %in% c(rois_cerebellum, rois_subcortical),
      sub %in% no_sex_subs
    )
  ) 

all_values1 %>%
  select(sub, PatientStatus, date) %>%
  distinct() %>%
  summarize(
    .by = c(sub, PatientStatus),
    n = n()
  ) %>%
  filter(
    n > 1
  )

table(all_values1$status, useNA = "ifany")

missing_status <- all_values1 %>%
  filter(
    is.na(status)
  ) %>% 
  select(sub, date) %>%
  distinct() %>%
  pull(sub) %>%
  unique()

missing_status %in% demo$sub

# 'all_values' is the good, clean data, _testonly is everyone else
all_values <- filter(all_values1, status == "train_and_test")
all_values_testonly <- filter(all_values1, status != "train_and_test")

# Just to check
missing_values <- all_values %>%
  filter(
    is.na(value)
  )

write_rds(all_values, "all_values_temp.rds")

# Reformat for training

reformatted <- all_values %>%
  ungroup() %>%
  mutate(
    # Remove non-alphanumeric characters
    across(starts_with("ROI"), ~str_remove_all(., "[/.+]")),
    aphasia = if_else(WAB_value < 93.8, "aphasia", "no_aphasia")
  ) %>%
  separate_wider_delim(ROI1, delim = "_", names = c("ROI1", "H1")) %>%
  separate_wider_delim(ROI2, delim = "_", names = c("ROI2", "H2")) 

# Check for weird characters
rois <- unique(c(reformatted$ROI1, reformatted$ROI2))

# Remove redundant columns, i.e. V2_V1 is redundant with V1_V2
lower_triangle <- expand_grid(roi1 = rois, roi2 = rois) %>%
  mutate(
    cx = str_glue("cx_{roi1}_{roi2}")
  ) %>%
  pivot_wider(names_from = roi2, values_from = cx) %>%
  column_to_rownames("roi1") %>%
  as.matrix()

# Remove the redundant values
lower_triangle[upper.tri(lower_triangle, diag = TRUE)] <- NA

lower_triangle2 <- lower_triangle %>%
  as_tibble(rownames = "roi1") %>%
  pivot_longer(-roi1, names_to = "roi2", values_drop_na = TRUE)

connections_to_keep <- lower_triangle2$value

# Check for non-alphanumeric characters - when this is empty, there are none
# left
str_subset(rois, "[^A-Za-z0-9]")

# What are the missing values? Are there any missing within-RH?
# table(reformatted$H1, is.nan(reformatted$value))

# missing_RH <- reformatted %>%
#   filter(
#     H1 == "R" & H2 == "R",
#     is.nan(value) | is.na(value)
#   )

# This is sub, date, demo cols ..., cx_ cols ...
reformatted2 <- reformatted %>%
  mutate(
    connection = str_glue("cx_{ROI1}_{ROI2}"),
  ) %>%
  filter(
    H1 == H2,
    connection %in% connections_to_keep
  ) %>%
  select(-ROI1, -ROI2, -H2) %>%
  pivot_wider(names_from = connection) %>%
  rename(
    hemi = H1
  ) %>%
  mutate(
    class = str_glue("{PatientStatus}_{hemi}H"),
  ) %>%
  arrange(sub, hemi)

qs2::qs_save(reformatted2, "all_data_reformatted.rds")

# Are there any NaNs in the data check
missing_values2 <- reformatted2 %>%
  select(-matches("EHI|WAB")) %>%
  group_by(sub, hemi, Sex, PatientStatus, age_at_stroke, LesionVolume) %>%
  nest() %>%
  mutate(
    nans = map_int(data, count_NaNs)
  ) %>%
  filter(
    nans > 0
  )

# This should be 0
sum(missing_values2$nans > 0)

# Create folds ====

set.seed(20007)

# How many lefties are there?
table(reformatted2$PatientStatus, reformatted2$EHI_value > 0) / 2

reformatted2 %>%
  select(sub, PatientStatus, EHI_value) %>%
  distinct() %>%
  summarize(
    .by = PatientStatus,
    n = sum(EHI_value <= 0)
  )

# Normed ====

missing_demographics <- all_values %>%
  ungroup() %>%
  select(sub, Sex, PatientStatus, EHI_age, LesionVolume, WAB_value) %>%
  distinct() %>%
  filter(
    when_any(
      is.na(EHI_age), 
      is.na(Sex),
      PatientStatus == "LHS" & (is.na(LesionVolume) | is.na(WAB_value))
    )
  )

left_join(missing_demographics, projects, by = join_by(sub))

# Create new balanced folds
folds2 <- all_values %>%
  filter_out(
    sub %in% missing_demographics$sub
  ) %>%
  ungroup() %>%
  select(sub, Sex, PatientStatus, EHI_age, LesionVolume, WAB_value) %>%
  distinct() %>%
  mutate(
    sub = factor(sub),
    sex_group = str_glue("{Sex}_{PatientStatus}")
  ) %>%
  fold(k = 3, id_col = "sub", cat_col = "sex_group", num_col = "EHI_age") %>%
  rename(
    fold = .folds
  ) %>%
  ungroup() %>%
  mutate(
    sub = as.character(sub),
    logLV = log10(LesionVolume),
    logLV_g = cut(logLV, 
                  breaks = quantile(logLV, na.rm = TRUE),
                  include.lowest = TRUE),
    WAB_value_g = cut(WAB_value, 
                      include.lowest = TRUE,
                      breaks = quantile(WAB_value, na.rm = TRUE))
  )

table(folds2$PatientStatus, folds2$logLV_g, useNA = "a")
table(folds2$PatientStatus, folds2$WAB_value_g, useNA = "a")

folds2 %>%
  filter(
    PatientStatus == "LHS",
    is.na(WAB_value_g)
  )

residual_WAB <- lm(WAB_value ~ LesionVolume, 
                   data = filter(folds2, PatientStatus == "LHS"))

residual_WAB_broom <- broom::augment(residual_WAB) %>%
  left_join(
    filter(folds2, PatientStatus == "LHS"), .
  ) %>%
  select(sub, LesionVolume, .resid, .std.resid) %>%
  rename(
    residual_wab = .resid,
    residual_wab_std = .std.resid
  )

ggplot(residual_WAB_broom, aes(x = LesionVolume, y = residual_wab_std)) +
  geom_point() +
  theme_bw()

folds3 <- left_join(folds2, residual_WAB_broom) %>%
  mutate(
    residual_wab_group = case_when(
      # Break into groups - lower values of WAB are worse, so negative 
      # residuals represent individuals who are doing worse than their lesion
      # size would predict, and vice versa. Set +/- 1 SD to be "as expected"
      # +/-1 SD was too wide.
      abs(residual_wab_std) < 0.5 ~ "rzd_wab_avg",
      residual_wab_std <= -0.5 ~ "rzd_wab_worse",
      residual_wab_std >= 0.5 ~ "rzd_wab_better",
    )
  )

folds3 %>%
  filter(
    PatientStatus == "LHS",
    is.na(WAB_value_g)
  )

all_corr_values <- all_values %>%
  mutate(
    across(starts_with("ROI"), ~str_remove_all(., "[/.+]")),
    aphasia = if_else(WAB_value < 93.8, "aphasia", "no_aphasia")
  ) %>%
  separate_wider_delim(ROI1, delim = "_", names = c("ROI1", "hemi")) %>%
  separate_wider_delim(ROI2, delim = "_", names = c("ROI2", "hemi2")) %>%
  filter(
    hemi == hemi2
  ) %>%
  mutate(
    connection = str_glue("cx_{ROI1}_{ROI2}"),
  ) %>%
  left_join(is_ln_lesioned2)

scalings <- all_corr_values %>%
  group_by(ROI1, ROI2, connection) %>%
  nest() %>%
  mutate(
    connection = str_glue("cx_{ROI1}_{ROI2}"),
    m_value = map_dbl(data, ~mean(.x$value)),
    sd_value = map_dbl(data, ~sd(.x$value)),
  ) %>%
  select(ROI1, ROI2, connection, m_value, sd_value) 

write_rds(scalings, "MSD_to_scale.rds")

normed_wide <- all_corr_values %>%
  left_join(scalings, by = join_by(ROI1, ROI2, connection)) %>%
  filter(
    connection %in% connections_to_keep 
  ) %>%
  filter_out(
    when_any(
      PatientStatus == "LHS" & hemi == "L",
      sub %in% missing_demographics$sub
    )
  ) %>%
  mutate(
    # Scale each value according to the calculated distribution
    value = (value - m_value) / sd_value
  ) %>%
  ungroup() %>%
  select(-hemi2, -m_value, -sd_value,  -starts_with("ROI")) %>%
  pivot_wider(names_from = connection) %>%
  mutate(
    class = str_glue("{PatientStatus}_{hemi}H"),
  ) %>%
  left_join(
    select(folds3, -logLV),
    # by = join_by(sub, PatientStatus)
  ) %>%
  select(
    # Move demographic columns to the left
    all_of(str_subset(colnames(.), "^cx_", negate = TRUE)), 
    everything()
  )

normed_wide$sub[is.na(normed_wide$logLV_g) & normed_wide$PatientStatus == "LHS"]
table(normed_wide$PatientStatus, normed_wide$logLV_g, useNA = "a")

write_csv(normed_wide, "input/group-all_parc-brainnetome_normed.csv")
py_save_object(normed_wide, "input/group-all_parc-brainnetome_normed.pickle")

# Transform the held-out subs

all_corr_values_testonly <- all_values_testonly %>%
  mutate(
    across(starts_with("ROI"), ~str_remove_all(., "[/.+]"))
  ) %>%
  separate_wider_delim(ROI1, delim = "_", names = c("ROI1", "hemi")) %>%
  separate_wider_delim(ROI2, delim = "_", names = c("ROI2", "hemi2")) %>%
  filter(
    hemi == hemi2
  ) %>%
  mutate(
    connection = str_glue("cx_{ROI1}_{ROI2}"),
  ) 

normed_wide_testonly <- all_corr_values_testonly %>%
  left_join(scalings, by = join_by(ROI1, ROI2, connection)) %>%
  filter(
    connection %in% connections_to_keep 
  ) %>%
  filter_out(
    when_any(
      PatientStatus == "LHS" & hemi == "L",
      sub %in% missing_demographics$sub
    )
  ) %>%
  mutate(
    # Scale each value according to the calculated distribution
    value = (value - m_value) / sd_value
  ) %>%
  ungroup() %>%
  select(-hemi2, -starts_with("ROI"), -m_value, -sd_value) %>%
  pivot_wider(names_from = connection, values_from = value) %>%
  mutate(
    class = str_glue("{PatientStatus}_{hemi}H"),
  ) %>%
  select(
    # Move demographic columns to the left
    all_of(str_subset(colnames(.), "^cx_", negate = TRUE)), 
    everything()
  )

write_csv(normed_wide_testonly, 
          "input/group-testonly_parc-brainnetome_normed.csv")

py_save_object(normed_wide_testonly, 
               "input/group-testonly_parc-brainnetome_normed.pickle")

# All rows

all_rows <- bind_rows(normed_wide, normed_wide_testonly)

write_csv(all_rows, 
          "input/group-allrows_parc-brainnetome_normed.csv")

py_save_object(all_rows, 
               "input/group-allrows_parc-brainnetome_normed.pickle")

all_rows %>%
  summarize(
    .by = c(PatientStatus, sub),
    n = length(unique(date))
  ) %>%
  filter(
    n > 1
  ) %>%
  arrange(desc(n))

all_rows %>%
  select(sub, PatientStatus, status) %>%
  distinct() %>%
  group_by(PatientStatus) %>%
  summarize(
    n = n()
  )

# Demographics table

demo_cols <- normed_wide %>%
  select(-starts_with("cx_"), -hemi, -class, -best_scan) %>%
  distinct()

demo_fmt <- demo_cols %>%
  group_by(group) %>%
  summarize(

    n = as.character(n()),

    sex = str_glue("{sum(Sex == 'M')} ({round(100 * sum(Sex == 'M') / as.integer(n))})% M"),
    age = str_glue("{round(mean(mri_age), 1)} (SD: {round(sd(mri_age), 1)})"),
    age_range = paste(round(range(mri_age), 1), collapse = " -- "),

    EHI = str_glue("{round(mean(EHI_value), 1)} (SD: {round(sd(EHI_value), 1)})"),
    EHI_range = paste(round(range(EHI_value)), collapse = " -- "),
    n_lefties = str_glue("{sum(EHI_value <= 0)} ({round(100 * sum(EHI_value <= 0) / as.integer(n))})"),

    age_str = str_glue("{round(mean(age_at_stroke), 1)} (SD: {round(sd(age_at_stroke), 1)})"),
    age_str_range = paste(round(range(age_at_stroke), 1), collapse = " -- "),

    chron_at_mri = str_glue("{round(mean(mri_chronicity_mo), 1)} (SD: {round(sd(mri_chronicity_mo), 1)})"),
    chron_range = paste(round(range(mri_chronicity_mo), 1), collapse = " -- "),
    
    WAB = str_glue("{round(mean(WAB_value), 1)} (SD: {round(sd(WAB_value), 1)})"),
    WAB_range = paste(range(WAB_value), collapse = " -- "),
    WAB_n_aphasia = str_glue("{sum(WAB_value <= 93.8)} ({round(100 * sum(WAB_value <= 93.8) / as.integer(n))})"),

    LV = str_glue("{round(mean(LesionVolume), 1)} (SD: {round(sd(LesionVolume), 1)})"),
    LV_range = paste(range(LesionVolume), collapse = " -- "),


  ) %>%
  pivot_longer(-group) %>%
  pivot_wider(names_from = group)


demo_cols_testonly <- normed_wide_testonly %>%
  select(-starts_with("cx_"), -hemi, -class, -best_scan) %>%
  distinct() %>%
  mutate(
    group = case_when(
      status == "test_long" & PatientStatus == "HC" ~ "HC_long",
      status == "test_long" & PatientStatus == "LHS" ~ "LHS_long",
      status == "test_only" & PatientStatus == "LHS" ~ "LplusRHS",
    )
  )

demo_fmt_testonly <- normed_wide_testonly %>%
  filter(
    # Don't double-count HCs
    hemi == "R"
  ) %>%
  group_by(sub) %>%
  mutate(
    across(c(Sex, LesionVolume), ~if_else(row_number() > 1, NA, .x))
  ) %>%
  select(-starts_with("cx_"), -hemi) %>%
  ungroup() %>%
  group_by(group) %>%
  summarize(

    n_ses = as.character(n()),
    n_indiv = as.character(length(unique(sub))),

    sex = str_glue("{sum(Sex == 'M', na.rm = TRUE)} ({round(100 * sum(Sex == 'M', na.rm = TRUE) / as.integer(n_indiv))})% M"),

    age = str_glue("{round(mean(mri_age), 1)} (SD: {round(sd(mri_age), 1)})"),
    age_range = paste(round(range(mri_age), 1), collapse = " -- "),

    # EHI = str_glue("{round(mean(EHI_value), 1)} (SD: {round(sd(EHI_value), 1)})"),
    # EHI_range = paste(round(range(EHI_value)), collapse = " -- "),
    # n_lefties = str_glue("{sum(EHI_value <= 0)} ({round(100 * sum(EHI_value <= 0) / as.integer(n))})"),

    age_str = str_glue("{round(mean(age_at_stroke), 1)} (SD: {round(sd(age_at_stroke), 1)})"),
    age_str_range = paste(round(range(age_at_stroke), 1), collapse = " -- "),

    chron_at_mri = str_glue("{round(mean(mri_chronicity_mo), 1)} (SD: {round(sd(mri_chronicity_mo), 1)})"),
    chron_range = paste(round(range(mri_chronicity_mo), 1), collapse = " -- "),
    
    WAB = str_glue("{round(mean(WAB_value), 1)} (SD: {round(sd(WAB_value), 1)})"),
    WAB_range = paste(range(WAB_value), collapse = " -- "),
    # WAB_n_aphasia = str_glue("{sum(WAB_value <= 93.8)} ({round(100 * sum(WAB_value <= 93.8) / as.integer(n))})"),

    LV = str_glue("{round(mean(LesionVolume), 1)} (SD: {round(sd(LesionVolume), 1)})"),
    LV_range = paste(range(LesionVolume), collapse = " -- "),

  ) %>%
  pivot_longer(-group) %>%
  pivot_wider(names_from = group)

demo_cols_testonly %>%
  filter(
    PatientStatus == "LHS",
    length(unique(date)) > 1
  ) %>%
  group_by(sub) %>%
  summarize(
    aph = sd(WAB_value <= 93.8),
    lv = sd(LesionVolume)
  )

demo_fmt_testonly2 <- normed_wide_testonly %>%
  filter(
    # Don't double-count HCs
    hemi == "R"
  ) %>%
  group_by(sub) %>%
  mutate(
    across(c(Sex, LesionVolume), ~if_else(row_number() > 1, NA, .x))
  ) %>%
  select(-starts_with("cx_"), -hemi) %>%
  group_by(sub, group) %>%
  summarize(
    n = n(),
    meanEHI = mean(EHI_value),
    n_aphasic = sum(WAB_value <= 93.8),
    LesionVolume = mean(LesionVolume, na.rm = TRUE)
  ) %>%
  mutate(
    pct_aphasic = n_aphasic / n
  )

demo_fmt_testonly2 %>%
  group_by(group) %>%
  summarize(

    n_ses = as.character(n()),
    n_indiv = as.character(length(unique(sub))),

    EHI = str_glue("{round(mean(meanEHI), 1)} (SD: {round(sd(meanEHI), 1)})"),
    EHI_range = paste(round(range(meanEHI)), collapse = " -- "),
    n_lefties = str_glue("{sum(meanEHI <= 0)} ({round(100 * sum(meanEHI <= 0) / as.integer(n_indiv))})"),

    LV = str_glue("{round(mean(LesionVolume), 1)} (SD: {round(sd(LesionVolume), 1)})"),
    LV_range = paste(range(LesionVolume), collapse = " -- "),

  ) %>%
  pivot_longer(-group) %>%
  pivot_wider(names_from = group)

