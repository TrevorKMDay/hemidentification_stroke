library(tidyverse)
library(reticulate)
library(groupdata2)
library(neurobase)
library(corrr)

setwd("~/MyDrive/Projects/hemisphere_fingerprinting_stroke/code/")

# Setup ====

## Participants =====

data_files <- list.files("../data/individual_matrices_andrew/",
                            pattern = "sub-.*.csv", full.names = TRUE)

lhs_subs_with_data <- data_files %>%
  str_subset("patient-1") %>%
  str_extract("sub-...")

## Projects 

projects <- read_csv("../data/projectnames.csv", show_col_types = FALSE) %>%
  mutate(
    sub = str_glue("sub-{SubjectCode}")
  ) %>%
  rename(date = ScanDate) %>%
  select(sub, date, projectname) %>%
  distinct() %>%
  filter_out(
    is.na(projectname)
  )

## ROIs

rois <- readxl::read_xlsx("../data/language_network_ids.xlsx")

rois_andrew <- readxl::read_excel("../data/BN_Atlas_274_combined.xlsx")

rois_cerebellum <- rois_andrew$Label[rois_andrew$Lobe == "Cerebellum"]
rois_subcortical <- rois_andrew$Label[rois_andrew$Lobe == "Subcortical"]

rois_cortical <- rois_andrew %>%
  filter_out(
    Label %in% rois_cerebellum |
    Label %in% rois_subcortical
  )

all_rois <- str_remove_all(rois_andrew$Label, "[/.+]")

# Get each ROI volume 

bn <- readnii(str_glue("~/MyDrive/Projects/brainnetome/my_version/",
                       "BN_Atlas_274_combined_2mm.nii.gz"))

bn_vox <- tibble(ID = 1:max(bn)) %>%
  mutate(
    roi_n_vox = map_int(ID, ~sum(bn == .x), .progress = TRUE),
    roi_vol_mm3 = roi_n_vox * 8,
  ) %>%
  left_join(
    select(rois_andrew, ID, Label)
  ) %>%
  select(ID, Label, roi_n_vox, roi_vol_mm3) %>%
  mutate(
    Label = str_remove_all(Label, "[/.+]")
  )

# Which network 

roi_assignments0 <- readxl::read_xlsx("../data/language_network_ids.xlsx") 

roi_assignments1 <- roi_assignments0 %>%
  select(ID, Label, Lobe, three_net_division) %>%
  mutate(
    Label = str_remove_all(Label, "[/.+]")
  )

roi_assignments <- tibble(Label = all_rois) %>%
  left_join(roi_assignments1, by = join_by(Label)) %>%
  mutate(

    three_net_division = case_when(
      Label %in% rois_cerebellum ~ "Cerebellum",
      Label %in% rois_subcortical ~ "Subcortical",
      .default = three_net_division
    ),

    language_division = case_when(
      three_net_division %in% c("Language", "Semantic") & Lobe == "Frontal" ~ 
        "langnet_front",
      # The one parietal area is a39rv, which is angular gyrus, so merge 
      # with parietal
      three_net_division %in% c("Language", "Semantic") & 
        Lobe %in% c("Temporal", "Parietal") ~ "langnet_temp",
      three_net_division %in% c("Language", "Semantic") & Lobe == "Limbic" ~ 
        "langnet_limbic",
      TRUE ~ three_net_division
    ),

    ID = row_number()

  )

# roi_assignments %>%
#   filter(
#     three_net_division == "Language",
#     Lobe == "Parietal"
#   )

# Load amount of lesion percentage

roi_labels <- str_remove_all(rois_andrew$Label, "[/.+]")

lesion0 <- read_csv("../data/lesions.txt", show_col_types = FALSE)
colnames(lesion0) <- c("sub", "date", roi_labels)

lhs_subs_with_data %in% str_glue("sub-{lesion0$sub}")

# Previously, removed non-cortical ROIs
lesion <- lesion0 %>%
  filter(
    sub %in% str_remove(lhs_subs_with_data, "sub-")
  ) 

lesion_long <- lesion %>%
  pivot_longer(-c(sub, date), names_to = "ROI", values_to = "lesioned_pct") %>%
  filter(
    # Only LH lesions
    # str_detect(ROI, "_L$")
  ) %>%
  mutate(
    sub = str_glue("sub-{sub}"),
  ) %>%
  left_join(
    select(bn_vox, ID, Label, roi_vol_mm3), 
    by = join_by(ROI == Label)
  ) %>%
  left_join(
    roi_assignments,
    by = join_by(ROI == Label, ID)
  ) %>%
  separate_wider_delim(ROI, delim = "_", names = c("ROI", "LR")) %>%
  mutate(
    total_lesion = lesioned_pct * roi_vol_mm3
  ) 

lesioned_ROIs <- unique(lesion_long$ROI[lesion_long$lesioned_pct > 0])

mean_by_roi <- lesion_long %>%
  group_by(ROI, LR) %>%
  summarize(
    n = n(),
    mean_pct = mean(lesioned_pct)
  )

mean_by_roi_wide <- mean_by_roi %>%
  select(-n) %>%
  pivot_wider(names_from = LR, values_from = mean_pct)

lesion_grouped1 <- lesion_long %>%
  group_by(sub, date, LR, three_net_division) %>%
  summarize(
    total_lesion_mm3 = round(sum(total_lesion)),
    total_volume_mm3 = sum(roi_vol_mm3)
  ) %>%
  mutate(
    total_division_pct = total_lesion_mm3 / total_volume_mm3,
    log_total_lesion = log10(total_lesion_mm3), 
    log_total_lesion_with0s = if_else(is.infinite(total_lesion_mm3), 0,
                                        total_lesion_mm3)
  )

lesion_grouped2 <- lesion_long %>%
  group_by(sub, date, LR, language_division) %>%
  mutate(
    total_lesion = lesioned_pct * roi_vol_mm3
  ) %>%
  summarize(
    total_lesion_mm3 = round(sum(total_lesion)),
    total_volume_mm3 = sum(roi_vol_mm3)
  ) %>%
  mutate(
    total_division_pct = total_lesion_mm3 / total_volume_mm3,
    log_total_lesion = log10(total_lesion_mm3), 
    log_total_lesion_with0s = if_else(is.infinite(total_lesion_mm3), 0,
                                        total_lesion_mm3)
  )

lesion_grouped <- bind_rows(
    rename(lesion_grouped1, division = three_net_division),
    rename(lesion_grouped2, division = language_division)
  ) %>%
  distinct() %>%
  arrange(sub) %>%
  group_by(sub, LR) %>%
  mutate(
    entire_corticalhemilesion_mm3 = sum(
        total_lesion_mm3[division %in% c("RoB", "Language", "Semantic")]
      ),
    log_total_lesion = replace_values(log_total_lesion, -Inf ~ NA)
  ) 
 
ggplot(lesion_grouped, aes(x = division, y = total_lesion_mm3, fill = LR)) +
  geom_boxplot() +
  theme_bw()

ggplot(lesion_grouped, aes(x = division, y = total_division_pct, fill = LR)) +
  geom_boxplot() +
  theme_bw()

ggplot(lesion_grouped, aes(x = division, y = log_total_lesion, fill = LR)) +
  geom_boxplot() +
  theme_bw()

ggplot(lesion_grouped, aes(x = division, y = sub, fill = total_division_pct)) +
  geom_tile() +
  scale_fill_gradient2() +
  theme_bw()

write_csv(lesion_grouped, "datasets/lesion_pct_by_division.csv")

lesion_grouped_wide <- lesion_grouped %>%
  select(sub, date, division, LR, entire_corticalhemilesion_mm3, 
          total_lesion_mm3) %>%
  pivot_wider(names_from = c(division), values_from = total_lesion_mm3) %>%
  left_join(projects, by = join_by(sub, date), relationship = "many-to-one")

RH_lesions <- lesion_grouped_wide %>%
  filter(
    LR == "R",
    entire_corticalhemilesion_mm3 > 0
  ) %>%
  arrange(desc(entire_corticalhemilesion_mm3)) %>%
  select(-date) %>%
  distinct()

lesion_grouped_total_byhemi <- lesion_grouped_wide %>%
  pivot_wider(id_cols = c(sub, date, projectname), names_from = LR, 
              values_from = entire_corticalhemilesion_mm3)

ggplot(lesion_grouped_total_byhemi, aes(x = L, y = R)) +
  geom_point() +
  scale_x_sqrt() +
  scale_y_sqrt() +
  theme_bw()

RH_lesions2 <- lesion_grouped_total_byhemi %>%
  filter(
    R > 0
  ) %>%
  arrange(desc(R))
