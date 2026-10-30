library(tidyverse)

setwd("~/MyDrive/Projects/hemisphere_fingerprinting_stroke/code/enrichment/")

# Brainnetome ====

rois <- readxl::read_xlsx("../../data/language_network_ids.xlsx")

rois_andrew <- readxl::read_excel("../../data/BN_Atlas_274_combined.xlsx")

rois_cerebellum <- rois_andrew$Label[rois_andrew$Lobe == "Cerebellum"]
rois_subcortical <- rois_andrew$Label[rois_andrew$Lobe == "Subcortical"]

rois_cortical <- rois_andrew %>%
  filter_out(
    Label %in% rois_cerebellum | Label %in% rois_subcortical
  )

# Glasser ====

networks0 <- read_csv("~/MyDrive/Projects/brainnetome/my_version/brainnetome_to_glassernetworks.csv", 
                      show_col_types = FALSE)

nets <- networks0 %>%
  filter(
    label %in% rois_cortical$LabelNoHemisphere
  ) %>%
  mutate(

    label = str_remove_all(label, "[/.+]"),

    net_short = recode_values(
      net,
      "Auditory" ~ "AudN",
      "CinguloOpercular" ~ "CON",
      "Default" ~ "DMN",
      "DorsalAttention" ~ "DAN",
      "Frontoparietal" ~ "FPN",
      "Language" ~ "LN",
      "OrbitoAffective" ~ "OAN",
      "PosteriorMultimodal" ~ "pMN",
      "Somatomotor" ~ "SMN",
      "VentralMultimodal" ~ "vMN",
      "Visual" ~ "VisN",
      "AudLang" ~ "AudLN",
      default = net
    )

  ) %>%
  rename(
    roi = label
  ) %>%
  select(roi, net, net_short)

nets_ordered <- nets %>%
  filter_out(
    net_short %in% c("OAN", "vMN")
  ) %>%
  arrange(net_short, roi) %>%
  mutate(
    order = row_number(),
  ) %>%
  group_by(net) %>%
  mutate(
    order_in_net = row_number()
  )

write_rds(nets_ordered, "nets-glasser_ordered.rds")

# For plotting
nets_squares <- nets_ordered %>%
  filter(
    order_in_net == 1 | order_in_net == max(order_in_net)
  ) %>%
  mutate(
    startend = if_else(order_in_net == 1, "start", "end")
  ) %>%
  pivot_wider(id_cols = c(net, net_short), names_from = startend,
              values_from = roi) %>%
  mutate(
    # Replace 1 with 0 so it doesn't overlap the first ROI
    start_order = match(start, nets_ordered$roi),
    end_order = match(end, nets_ordered$roi)
  )

write_rds(nets_squares, "nets-glasser_squares.rds")

# Yeo ====

yeo_networks0 <- readxl::read_xlsx("../../data/language_network_ids.xlsx") %>%
  select(ID, Label, starts_with("Yeo"))

yeo_networksL <- yeo_networks0 %>%
  filter(
    str_detect(Label, "_L$")
  ) %>%
  rename(
    label = Label,
    net = Yeo_RSN
  ) %>%
  select(-Yeo_RSN_ADJ) %>%
  mutate(
    label = str_remove(label, "_L")
  )

networks0$label[!(networks0$label %in% networks0$label)]

yeo_nets <- yeo_networksL %>%
  filter(
    label %in% str_remove(rois_cortical$Label, "_[LR]$")
  ) %>%
  mutate(

    label = str_remove_all(label, "[/.+]"),

    net_short = recode_values(
      net,
      "Auditory" ~ "AudN",
      "CinguloOpercular" ~ "CON",
      "Default" ~ "DMN",
      "DorsalAttention" ~ "DAN",
      "Frontoparietal" ~ "FPN",
      "Language" ~ "LN",
      "OrbitoAffective" ~ "OAN",
      "PosteriorMultimodal" ~ "pMN",
      "Somatomotor" ~ "SMN",
      "VentralAttention" ~ "VAN",
      "VentralMultimodal" ~ "vMN",
      "Visual" ~ "VisN",
      "AudLang" ~ "AudLN",
      default = net
    )

  ) %>%
  rename(
    roi = label
  ) %>%
  select(roi, net, net_short) 

yeo_nets %>%
  filter(
    net %in% c("Subcortical", "Limbic")
  )

yeo_nets_ordered <- yeo_nets %>%
  filter_out(
    net_short %in% c("Subcortical")
  ) %>%
  arrange(net_short, roi) %>%
  mutate(
    order = row_number(),
  ) %>%
  group_by(net) %>%
  mutate(
    order_in_net = row_number()
  )

write_rds(yeo_nets_ordered, "nets-yeo_ordered.rds")

# For plotting
yeo_nets_squares <- yeo_nets_ordered %>%
  filter(
    order_in_net == 1 | order_in_net == max(order_in_net)
  ) %>%
  mutate(
    startend = if_else(order_in_net == 1, "start", "end")
  ) %>%
  pivot_wider(id_cols = c(net, net_short), names_from = startend,
              values_from = roi) %>%
  mutate(
    # Replace 1 with 0 so it doesn't overlap the first ROI
    start_order = match(start, yeo_nets_ordered$roi),
    end_order = match(end, yeo_nets_ordered$roi)
  )

write_rds(yeo_nets_squares, "nets-yeo_squares.rds")

# Plot 

library(ggseg)
library(ggsegBrainnetome)

clean_names <- rois_andrew %>%
  select(LabelNoHemisphere) %>%
  distinct() %>%
  mutate(
    label_clean = str_remove_all(LabelNoHemisphere, "[/.+]"),
    region = str_glue("{LabelNoHemisphere}_L")
  ) %>%
  left_join(
    brainnetome()$core
  )

yeo_nets_plotme <- yeo_nets %>%
  left_join(
    clean_names,
    by = join_by(roi == label_clean)
  ) %>%
  mutate(
    hemi = "left",
  )

ggplot(yeo_nets_plotme) + 
  geom_brain(
    atlas = brainnetome(), 
    position = position_brain(hemi ~ view), 
    hemi = "left",
    aes(fill = net)
  ) +
  scale_alpha_continuous(limits = c(0.0, 1)) +
  theme_bw() +
  theme(legend.position = "bottom")

# Yeo v2 ====

yeo_networksL2 <- yeo_networks0 %>%
  filter(
    str_detect(Label, "_L$")
  ) %>%
  rename(
    label = Label,
    net = Yeo_RSN_ADJ
  ) %>%
  select(-Yeo_RSN) %>%
  mutate(
    label = str_remove(label, "_L")
  )

yeo_nets2 <- yeo_networksL2 %>%
  filter(
    label %in% str_remove(rois_cortical$Label, "_[LR]$")
  ) %>%
  mutate(

    label = str_remove_all(label, "[/.+]"),

    net_short = recode_values(
      net,
      "Auditory" ~ "AudN",
      "CinguloOpercular" ~ "CON",
      "Default" ~ "DMN",
      "DorsalAttention" ~ "DAN",
      "Frontoparietal" ~ "FPN",
      "Language" ~ "LN",
      "OrbitoAffective" ~ "OAN",
      "PosteriorMultimodal" ~ "pMN",
      "Somatomotor" ~ "SMN",
      "VentralAttention" ~ "VAN",
      "VentralMultimodal" ~ "vMN",
      "Visual" ~ "VisN",
      "AudLang" ~ "AudLN",
      default = net
    )

  ) %>%
  rename(
    roi = label
  ) %>%
  select(roi, net, net_short) 

yeo_nets2 %>%
  filter(
    net %in% c("Subcortical", "Limbic")
  )

yeo_nets2_ordered <- yeo_nets2 %>%
  filter_out(
    net_short %in% c("Subcortical")
  ) %>%
  arrange(net_short, roi) %>%
  mutate(
    order = row_number(),
  ) %>%
  group_by(net) %>%
  mutate(
    order_in_net = row_number()
  )

write_rds(yeo_nets2_ordered, "nets-yeo2_ordered.rds")

# For plotting
yeo_nets2_squares <- yeo_nets2_ordered %>%
  filter(
    order_in_net == 1 | order_in_net == max(order_in_net)
  ) %>%
  mutate(
    startend = if_else(order_in_net == 1, "start", "end")
  ) %>%
  pivot_wider(id_cols = c(net, net_short), names_from = startend,
              values_from = roi) %>%
  mutate(
    # Replace 1 with 0 so it doesn't overlap the first ROI
    start_order = match(start, yeo_nets_ordered$roi),
    end_order = match(end, yeo_nets_ordered$roi)
  )

write_rds(yeo_nets2_squares, "nets-yeo2_squares.rds")

# Plot 

library(ggseg)
library(ggsegBrainnetome)

yeo_nets2_plotme <- yeo_nets2 %>%
  left_join(
    clean_names,
    by = join_by(roi == label_clean)
  ) %>%
  mutate(
    hemi = "left",
  )

ggplot(yeo_nets2_plotme) + 
  geom_brain(
    atlas = brainnetome(), 
    position = position_brain(hemi ~ view), 
    hemi = "left",
    aes(fill = net)
  ) +
  scale_alpha_continuous(limits = c(0.0, 1)) +
  theme_bw() +
  theme(legend.position = "bottom")
