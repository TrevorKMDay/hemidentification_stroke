library(tidyverse)
library(patchwork)
library(viridis)

# Setup ====

setwd("~/MyDrive/Projects/hemisphere_fingerprinting_stroke/code/enrichment/")

# Load necessary functions
source("00-enrichment_functions.R")

nets_squares <- read_rds("nets-yeo2_squares.rds") %>%
  mutate(
    start_order = start_order - 0.5,
    end_order = end_order + 0.5
  )

# This comes in with OAN (n=1) and vMN (n=3) omitted
nets_ordered <- read_rds("nets-yeo2_ordered.rds")


brainnetome_nets <- readxl::read_xlsx(str_glue("~/MyDrive/Projects/",
                                               "hemisphere_fingerprinting_stroke/", 
                                               "data/language_network_ids.xlsx")) %>%
  filter(
    str_detect(Label, "_L")
  ) %>%
  mutate(
    label = str_remove_all(Label, "[/.+]|_L"),
  ) %>%
  rename(
    net = Yeo_RSN
  ) %>%
  select(ID, label, net) %>%
  arrange(net, label) %>%
  filter(
    label %in% nets_ordered$roi
  ) %>%
  mutate(
    # Fix order
    order = row_number()
  ) 

# Load the previously calculated p-values
# This takes awhile, so only re-run script 01 when needed to update this
test_1class_ps <- read_rds("hemi.stroke_ps.rds") %>%
  mutate(
    feature = str_remove(feature, "cx_")
  ) %>%
  separate_wider_delim(feature, delim = "_", names = c("ROI1", "ROI2"),
                        cols_remove = FALSE) %>%
  filter(
    ROI1 %in% nets_ordered$roi,
    ROI2 %in% nets_ordered$roi
  )

## Features p < .001 ====

# Get the most important features with p values < .001 and < .01 for info/
# exploratory purposes

high_value_p001 <- test_1class_ps %>%
  filter(
    p_hemi < .001
  ) 

high_value_nodes_p001 <- tibble(rois = c(high_value_p001$ROI1,
                                         high_value_p001$ROI2)) %>%
  group_by(rois) %>%
  summarize(
    n = n()
  ) %>%
  arrange(desc(n)) %>%
  left_join(
    select(brainnetome_nets, label, net), 
    join_by(rois == label)
  )

length(unique(high_value_nodes_p001$rois))

## Features p < .01 ====

high_value_p01 <- test_1class_ps %>%
  filter(
    p_hemi <= .01
  ) 

n_hv_nodes_p01 <- length(unique(high_value_nodes_p001$rois))

high_value_nodes_p01 <- tibble(rois = c(high_value_p01$ROI1,
                                        high_value_p01$ROI2)) %>%
  group_by(rois) %>%
  summarize(
    n = n()
  ) %>%
  mutate(
    p = round(n / n_hv_nodes_p01, 3)
  ) %>%
  arrange(desc(n))

# Visualize the distriution of p values
ggplot(test_1class_ps, aes(x = -1 * log_p_hemi)) +
  geom_histogram(binwidth = 0.1, boundary = 0) +
  geom_vline(xintercept = -1 * log(c(0.05, 0.01))) +
  theme_bw()

lt_p05 <- sum(test_1class_ps$p_hemi < .05)
lt_p01 <- sum(test_1class_ps$p_hemi < .01)
lt_p001 <- sum(test_1class_ps$p_hemi < .001)

# Visualizations ====

# raw p-value ====

g_input <- test_1class_ps %>%
  select(feature, scaling1, p_hemi) %>%
  separate_wider_delim(feature, "_", names = c("ROI1", "ROI2"), 
                       cols_remove = FALSE) %>%
  mutate(
    across(c(ROI1, ROI2), ~factor(.x, levels = brainnetome_nets$label))
  ) %>%
  arrange(ROI1, ROI2)

g_input %>%
  filter(
    is.na(ROI1) | is.na(ROI2)
  )

g1 <- graph_to_tibble(g_input, "scaling1") %>%
  mutate(
    ROI1_order = map_int(ROI1, ~ brainnetome_nets$order[brainnetome_nets$label == .x]),
    ROI2_order = map_int(ROI2, ~ brainnetome_nets$order[brainnetome_nets$label == .x])
  )

# pepper plot ====

p_input <- g_input %>%
  mutate(
    sig = p_hemi < .01
  )

pepper <- graph_to_tibble(p_input, "sig") %>%
  mutate(
    ROI1_order = map_int(ROI1, ~ brainnetome_nets$order[brainnetome_nets$label == .x]),
    ROI2_order = map_int(ROI2, ~ brainnetome_nets$order[brainnetome_nets$label == .x])
  ) %>%
  rename(
    sig = value
  ) %>%
  mutate(
    sig = as.logical(sig)
  )

# Quickly validate pepper matrix
ggplot(pepper, aes(x = ROI1_order, y = ROI2_order, fill = sig)) +
  geom_tile()

brainnetome_nets$order[!(brainnetome_nets$order %in% pepper$ROI1_order)]
brainnetome_nets$label[!(brainnetome_nets$order %in% pepper$ROI1_order)]

## Colored plot ====

labels <- c("", nets_squares$net_short)

ld1_plot <- ggplot(NULL) +
  geom_tile(
    data = g1,
    aes(x = ROI1_order, y = ROI2_order, fill = value),
    linewidth = 0.2
  ) +
  scale_x_continuous(labels = labels,
                     breaks = c(1, nets_squares$end_order),
                     minor_breaks = NULL) +
  scale_y_continuous(labels = labels,
                     breaks = c(1, nets_squares$end_order),
                     minor_breaks = NULL,
                     trans = "reverse") +
  scale_fill_gradient2(na.value = NA, limits = c(NA, NA)) +
  geom_rect(
    data = nets_squares,
    aes(xmin = start_order, xmax = end_order,
        ymin = start_order, ymax = end_order,
        color = net_short),
    fill = NA
  ) +
  coord_equal() +
  theme_bw() +
  theme(axis.ticks = NULL, axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = "ROI 1", y = "ROI 2", color = "Network", fill = "Scaling",
       title = "(A)") +
  guides(color = guide_legend(ncol = 2))

## Pepper plot =====

pepper_plot <- ggplot(NULL) +
  geom_tile(
    data = pepper,
    aes(x = ROI1_order, y = ROI2_order, fill = sig, alpha = sig),
    linewidth = 0.2
  ) +
  scale_x_continuous(labels = labels,
                     breaks = c(1, nets_squares$end_order),
                     minor_breaks = NULL) +
  scale_y_continuous(labels = labels,
                     breaks = c(1, nets_squares$end_order),
                     minor_breaks = NULL,
                     trans = "reverse") +
  scale_fill_manual(values = c("white", "black"), guide = "none",
                    na.value = "white") +
  geom_rect(
    data = nets_squares,
    aes(xmin = start_order, xmax = end_order,
        ymin = start_order, ymax = end_order,
        color = net_short),
    fill = NA
  ) +
  coord_equal() +
  theme_bw() +
  theme(axis.ticks = NULL, axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "ROI 1", y = NULL, color = "Network", fill = "Sig.",
       alpha = "Sig.",
       title = "(B)")+
  guides(color = guide_legend(ncol = 2), alpha = "none")

ld1_plot + pepper_plot + plot_layout(axis_titles = "collect")

# Bootstrap models ====

true_nets <- pepper %>%
  filter(
    ROI1 %in% nets_ordered$roi,
    ROI2 %in% nets_ordered$roi
  ) %>%
  summarize_tibble_by_net(., nets_ordered) %>%
  mutate(
    prop = replace(prop, prop == 0, NA)
  ) %>%
  remove_tri()

bootstrapped <- pepper %>%
  filter(
    ROI1 %in% nets_ordered$roi,
    ROI2 %in% nets_ordered$roi
  ) %>% 
  bootstrap_pmatrix(., nets_ordered)

bootstrapped1 <- remove_tri(bootstrapped) %>%
  mutate(

    sig = replace(sig, sig == 0, NA),

    prop = replace(prop, prop == 0, NA),

    sig_label = case_when(
      pvalue < .001 ~ "***",
      pvalue < .01 ~ "**",
      pvalue < .05 ~ "*",
      TRUE ~ ""
    ),

    pvalue = replace(pvalue, pvalue == 0, 0.0003)

  )

bootstrapped1$p_fdr <- p.adjust(bootstrapped1$pvalue, method = "fdr")

bootstrapped1 %>%
  filter(
    pvalue < .05
  )

net_plot <- ggplot(bootstrapped1, aes(x = net1, y = net2, fill = 100 * prop)) +
  geom_tile() +
  geom_text(aes(label = sig_label)) +
  viridis::scale_fill_viridis(limits = c(0, NA)) +
  scale_y_discrete(limits = rev) +
  coord_equal() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, title = "(C)", fill = "%")

bootstrapped1 %>%
  filter(
    p_fdr < .05
  )

net_plot
