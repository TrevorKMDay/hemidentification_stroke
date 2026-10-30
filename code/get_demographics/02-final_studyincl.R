library(tidyverse)
library(readxl)
library(writexl)

setwd("~/MyDrive/Projects/hemisphere_fingerprinting_stroke/code/get_demographics/")

# Get notes ===

notes_file <- "crl_pptnotes_for_studyincl.xlsx"

patients <- read_xlsx(notes_file, sheet = "Patients", skip = 2)
controls <- read_xlsx(notes_file, sheet = "Controls", skip = 2)
case_studies <- read_xlsx(notes_file, sheet = "Case Studies", skip = 2)
patients2 <- read_xlsx(notes_file, sheet = "Pre-BUILD Patients", skip = 2)
controls2 <- read_xlsx(notes_file, sheet = "Pre-BUILD Controls", skip = 2)

all_notes <- bind_rows(patients, controls, case_studies, patients2, 
                        controls2) %>%
  rename(
    sub = `Subject ID`
  )

# Get dataset ====

my_data_0 <- read_csv("../input/group-allrows_parc-brainnetome_normed.csv",
                      show_col_types = FALSE) 

my_data <- my_data_0 %>%
  select(sub, PatientStatus, status) %>%
  distinct() %>%
  mutate(
    sub = str_remove(sub, "sub-")
  )

my_data_checked <- my_data %>%
  mutate(
    in_spreadsheet = sub %in% all_notes$sub
  ) %>%
  left_join(all_notes) 

writexl::write_xlsx(my_data_checked, "trevor_inclusion.xlsx")

table(my_data_checked$PatientStatus, my_data_checked$in_spreadsheet)
