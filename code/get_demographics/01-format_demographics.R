library(tidyverse)
library(reticulate)
library(groupdata2)
library(readxl)
library(googlesheets4)

setwd("~/MyDrive/Projects/hemisphere_fingerprinting_stroke/code/get_demographics/")

# Setup ====

exclusions <- read_sheet("13cEBidZSOCNMWOI4JxfaS2mIxahjjLciv7xKKT3hH-A") %>%
  mutate(
    sub = paste0("sub-", sub)
  )

data_files <- list.files("../../data/individual_matrices_andrew/",
                            pattern = "sub-.*.csv", full.names = TRUE)

data_files <- read_csv("../participants_with_mri.csv", show_col_types = FALSE)

# SHM was accidentally entered as 'CAT'
participants_with_imaging <- sort(c(unique(data_files$sub), "sub-SHM"))

write_csv(tibble(sub = str_remove(participants_with_imaging, "sub-")),
          "subs_demo.csv")

reduce_values <- function(vec, avg=FALSE, min=FALSE) {

  uniq <- unique(na.omit(vec))

  if (avg) 
    uniq <- mean(uniq)
  
  if (min)
    uniq <- min(uniq)

  if (length(uniq) == 0)
    return(NA)
  else
    return(uniq)

}

# Load data ====

na_vals <- c("-99", "Unknown", "N/A")

demo0_1 <- read_xls("data_exporter/demo_allalltime_260814.xls", na = na_vals)
demo0_2 <- read_xls("data_exporter/demo_allalltime_addtl_261006.xls", 
                    na = na_vals)
demo0_3 <- read_xls("data_exporter/demo_allalltime_addtl_261014.xls", 
                    na = na_vals)
demo0_4 <- read_xls("data_exporter/demo_allalltime_addtl_261016.xls", 
                    na = na_vals)
ehi_wab0 <- read_xls("data_exporter/demo_EHIandWAB_260814.xls", na = na_vals)

# Extra rows of data from elsewhere
manual_patch <- read_csv("manual_patch.csv", col_types = "ccccccDncn", 
                          show_col_types = FALSE)

demo0 <- bind_rows(demo0_1, demo0_2, demo0_3, demo0_4, ehi_wab0) %>%
  mutate(

    sub = str_glue("sub-{SubjectCode}"), .before = SubjectCode,

    PatientStatus = if_else(PatientStatus == 1, "LHS", "HC"),
    
    RawTestItemNumber = case_when(
        TestNameLong == "Demographics" & RawTestItemNumber == "1.1" ~ "DOS",
        TestNameLong == "Demographics" & RawTestItemNumber == "1.2" ~ "DOB",
        TestNameLong == "Demographics" & RawTestItemNumber == "1.5" ~ "Sex",
      ),

    # Fix a bug
    PatientStatus = case_when(
      sub == "sub-PSN" ~ "LHS",
      .default = PatientStatus
    ),

    # A new participant label was erroneously created for SHM
    # SHM and CAT are the same person
    sub = replace_values(sub, from = "CAT", to = "SHM")

  ) %>%
  bind_rows(
    # These rows are filled in post-renaming for human legibility
    manual_patch
  ) %>% 
  arrange(sub) %>%
  fill(PatientStatus, .by = sub, .direction = "downup") 

# Missing data
data_files %>%
  filter(
    !(sub %in% demo0$sub)
  ) %>%
  select(-f) %>%
  write_csv("current_missing_demo.csv")

dob <- demo0 %>%
  select(sub, RawTestItemNumber, RawTestResponseValue) %>%
  filter(
    RawTestItemNumber %in% c("DOB", "DOS"),
    sub %in% participants_with_imaging
  ) %>%
  group_by(sub, RawTestItemNumber) %>%
  mutate(
    RawTestResponseValue = RawTestResponseValue %>%
      replace_values(
        c("Unknown", "Unscored") ~ NA
      )
  ) %>%
  reframe(
    # Take date of first stroke; and earliest birthday when two values
    RawTestResponseValue = reduce_values(RawTestResponseValue, min = TRUE)
  ) %>%
  pivot_wider(names_from = RawTestItemNumber, 
              values_from = RawTestResponseValue) %>%
  mutate(
    age_at_stroke = interval(DOB, DOS) / years(1)
  )

sex_etc <- demo0 %>%
  select(sub, Sex, PatientStatus, RawTestItemNumber, RawTestResponseValue) %>%
  filter(
    sub %in% participants_with_imaging,
    RawTestItemNumber == "Sex"
  ) %>%
  mutate(
    # Fill in missing Sex values from the response row where needed
    Sex = case_when(
      tolower(RawTestResponseValue) == "female" ~ "F",
      tolower(RawTestResponseValue) == "male" ~ "M",
      .default = Sex
    )
  ) %>%
  distinct() %>%
  group_by(sub, PatientStatus) %>%
  reframe(
    Sex = reduce_values(Sex)
  ) %>%
  left_join(dob, by = join_by(sub))

# Check for multiple values
sex_etc %>%
  group_by(sub) %>%
  filter(
    n() > 1
  )

sex_etc %>%
  filter(
    is.na(Sex) | (PatientStatus == "LHS" & is.na(age_at_stroke))
  )

# Find best-matching EHI and WAB

scan_dates <- read_csv("../../data/scan_dates.csv", show_col_types = FALSE) %>%
  group_by(SubjectCode) %>%
  mutate(
    sub = str_glue("sub-{SubjectCode}"), .before = SubjectCode,
    sub = replace_values(sub, "sub-CAT" ~ "sub-SHM")
  ) %>%
  arrange(SubjectCode, ScanDate) %>%
  ungroup() %>%
  left_join(
    distinct(select(demo0, sub, PatientStatus)),
    by = join_by(sub),
    relationship = "many-to-one"
  ) %>%
  select(sub, PatientStatus, ScanDate)

write_csv(scan_dates, "scan_dates.csv")

basic_demo <- left_join(sex_etc, scan_dates, 
                        by = join_by(sub, PatientStatus)) %>%
  mutate(
    mri_age = interval(DOB, ScanDate) / years(1),
    mri_chronicity_mo = interval(DOS, ScanDate) / months(1)
  )

# BX ====

ehi_wab <- demo0 %>%
  filter(
    when_any(
      SummaryScoreName == "LateralityIndex",
      TestNameLong == "WAB Aphasia Quotient",
    ),
    sub %in% participants_with_imaging
  ) %>%
  mutate(

    name = if_else(TestNameLong == "WAB Aphasia Quotient", "WAB", 
                    SummaryScoreName),

    value = if_else(TestNameLong == "WAB Aphasia Quotient", 
                    as.numeric(RawTestResponseValue), SummaryScoreValue),

    # Another manual repaire
    value = case_when(
      sub == "sub-RPC" & TestDate == ymd("2019-10-23") & name == "WAB" ~ 99.8,
      .default = value
    )

  ) %>%
  filter_out(
    # Skip Anna's post-morbid EHIs (for LHS, HCs are fine)
    when_all(
      PatientStatus == "LHS",
      ProjectName %in% c("Unknown", "Vispat"),
      SummaryScoreName == "LateralityIndex"
    )
  ) %>%
  select(sub, TestDate, name, value, ProjectName) %>%
  distinct() %>%
  left_join(scan_dates, by = join_by(sub), relationship = "many-to-many") %>%
  mutate(
    test_scan_diff = interval(ScanDate, TestDate) / days(1)
  ) %>%
  group_by(sub, TestDate) %>%
  fill(ProjectName, .direction = "downup") %>%
  distinct()

ehi_wab_closest <- ehi_wab %>%
  group_by(sub, name, ScanDate) %>%
  filter(
    abs(test_scan_diff) == min(abs(test_scan_diff))
  ) %>%
  ungroup() %>%
  select(sub, ScanDate, TestDate, name, value) %>%
  mutate(
    name = replace_values(name, "LateralityIndex" ~ "EHI"),
  ) %>%
  rename(
    date = TestDate,
  ) %>%
  pivot_wider(id_cols = c(sub, ScanDate), 
              values_from = c(value, date),
              names_glue = "{name}_{.value}",
              names_vary = "slowest") %>%
  left_join(dob, by = join_by(sub)) %>%
  mutate(
    EHI_age = interval(DOB, EHI_date) / years(1),
    WAB_age = interval(DOB, WAB_date) / years(1),
  )

# Lesion volume 

lv <- read_xlsx("../../data/classifier_output1.xlsx") %>%
  mutate(
    PatientStatus = case_when(
      PatientStatus == 1 ~ "LHS",
      PatientStatus == 0 ~ "HC",
      PatientStatus == -99 ~ NA
    )
  ) %>%
  mutate(
    sub = str_glue("sub-{SubjectCode}"), .before = SubjectCode,
  ) %>%
  bind_rows(demo0) %>%
  filter(
    PatientStatus == "LHS",
    !is.na(LesionVolume)
  ) %>%
  select(sub, PatientStatus, LesionVolume) %>%
  arrange(sub) %>%
  distinct() 

# Final demographics ====

demo <- left_join(basic_demo, ehi_wab_closest, 
                  join_by(sub, DOB, DOS, age_at_stroke, ScanDate)) %>%
  left_join(
    select(lv, -PatientStatus),
    by = join_by(sub)
  ) %>%
  select(-DOB, -DOS, -ends_with("Date"), ScanDate) %>%
  mutate(
    # Manually fill in missing values
    Sex = if_else(sub %in% c("sub-ICY", "sub-STP"), "M", Sex),
    logLV = log10(LesionVolume),
    group = if_else(sub %in% exclusions$sub, "LHS_plusRH", PatientStatus)
  )

missing_all <- demo %>%
  filter(
    when_any(
      is.na(Sex),
      is.na(mri_age),
      is.na(EHI_value)
    )
  )

missing_lhs <- demo %>%
  filter(
    PatientStatus == "LHS",
    when_any(
      is.na(age_at_stroke),
      is.na(mri_chronicity_mo),
      is.na(WAB_value),
      is.na(LesionVolume)
    )
  )

subs_missing_demo <- participants_with_imaging[!(participants_with_imaging %in% demo$sub)]
subs_missing_demo

write_csv(demo, "../datasets/demographics.csv")

projects <- read_csv("../../data/projectnames.csv", show_col_types = FALSE) %>%
  filter(
    SubjectCode %in% str_remove(subs_missing_demo, "sub-") |
      SubjectCode %in% missing_lhs$sub
  )

# Format for publication

demo_nolong <- demo %>%
  group_by(sub) %>%
  filter(
    ScanDate == max(ScanDate)
  )

demo_fmt <- demo_nolong %>%
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
 
clipr::write_clip(demo_fmt)

# Check against exclusions

crl_notes <- read_xlsx("crl_pptnotes_for_studyincl.xlsx", skip = 2)

crl_notes_this_sample <- crl_notes %>%
  filter(
    str_glue("sub-{`Subject ID`}") %in% demo$sub
  )

write_csv(crl_notes_this_sample, "trevor_sample_notes.csv")

# demo 2 

demo_longonly <- demo %>%
  group_by(sub) %>%
  filter(
    n() > 1,
    ScanDate < max(ScanDate)
  ) %>%
  select(-EHI_age, -WAB_age) %>%
  mutate(

    EHI_value = mean(EHI_value),

    any_aphasia = any(WAB_value < 93.8),

    # Remove all but one value per ppt so averages don't count some multiple 
    # times
    across(c(Sex, age_at_stroke, EHI_value, WAB_value, LesionVolume, any_aphasia),
            ~if_else(row_number() > 1, NA, .x)),


  )

demo_fmt_2 <- demo_longonly %>%
  group_by(group) %>%
  summarize(

    n = as.character(length(unique(sub))),
    sessions = as.character(n()),

    sex = str_glue("{sum(Sex == 'M', na.rm = TRUE)} ({round(100 * sum(Sex == 'M', na.rm = TRUE) / as.integer(n))})% M"),
    age = str_glue("{round(mean(mri_age, na.rm = TRUE), 1)} (SD: {round(sd(mri_age, na.rm = TRUE), 1)})"),
    age_range = paste(round(range(mri_age, na.rm = TRUE), 1), collapse = " -- "),

    EHI = str_glue("{round(mean(EHI_value, na.rm = TRUE), 1)} (SD: {round(sd(EHI_value, na.rm = TRUE), 1)})"),
    EHI_range = paste(round(range(EHI_value, na.rm = TRUE)), collapse = " -- "),
    n_lefties = str_glue("{sum(EHI_value <= 0, na.rm = TRUE)} ({round(100 * sum(EHI_value <= 0, na.rm = TRUE) / as.integer(n))})"),

    age_str = str_glue("{round(mean(age_at_stroke, na.rm = TRUE), 1)} (SD: {round(sd(age_at_stroke, na.rm = TRUE), 1)})"),
    age_str_range = paste(round(range(age_at_stroke, na.rm = TRUE), 1), collapse = " -- "),

    chron_at_mri = str_glue("{round(mean(mri_chronicity_mo, na.rm = TRUE), 1)} (SD: {round(sd(mri_chronicity_mo, na.rm = TRUE), 1)})"),
    chron_range = paste(round(range(mri_chronicity_mo, na.rm = TRUE), 1), collapse = " -- "),
    
    WAB = str_glue("{round(mean(WAB_value, na.rm = TRUE), 1)} (SD: {round(sd(WAB_value, na.rm = TRUE), 1)})"),
    WAB_range = paste(range(WAB_value, na.rm = TRUE), collapse = " -- "),
    WAB_n_aphasia = str_glue("{sum(any_aphasia, na.rm = TRUE)} ({round(100 * sum(any_aphasia, na.rm = TRUE) / as.integer(n))})"),

    LV = str_glue("{round(mean(LesionVolume, na.rm = TRUE), 1)} (SD: {round(sd(LesionVolume, na.rm = TRUE), 1)})"),
    LV_range = paste(range(LesionVolume, na.rm = TRUE), collapse = " -- "),


  ) %>%
  pivot_longer(-group) %>%
  pivot_wider(names_from = group)

demo_fmt_2

clipr::write_clip(demo_fmt_2)

aphasia_count <- demo_longonly2 %>%
  group_by(sub) %>%
  filter_out(
    PatientStatus == "HC"
  ) %>%
  summarize(
    n = n(),
    aphasia = sum(WAB_value < 93.8)
  ) %>%
  mutate(
    p = aphasia / n 
  ) %>%
  arrange(desc(p))
