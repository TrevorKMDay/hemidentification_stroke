library(tidyverse)

setwd("~/MyDrive/Projects/hemisphere_fingerprinting_stroke/code/alternative_algorithms")

results <- read_csv("tsne_results_n2.csv", show_col_types = FALSE)

ggplot(results, aes(x = tsne_0, y = tsne_1)) +
  geom_point(aes(color = class)) +
  geom_line(aes(group = sub)) +
  theme_bw()

results_n3 <- read_csv("tsne_results_n3.csv", show_col_types = FALSE)

ggplot(results_n3, aes(x = tsne_0, y = tsne_1)) +
  geom_point(aes(color = class)) +
  geom_line(aes(group = sub)) +
  theme_bw()

ggplot(results_n3, aes(x = tsne_0, y = tsne_2)) +
  geom_point(aes(color = class)) +
  geom_line(aes(group = sub)) +
  theme_bw()

ggplot(results_n3, aes(x = tsne_1, y = tsne_2)) +
  geom_point(aes(color = class)) +
  geom_line(aes(group = sub)) +
  theme_bw()
