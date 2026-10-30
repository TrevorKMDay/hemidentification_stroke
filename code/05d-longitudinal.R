library(tidyverse)

ld1ld2 <- read_csv("datasets/ld1ld2_results.csv", show_col_types = FALSE)

longitudinal_ppts <- unique(ld1ld2$sub[ld1ld2$status == "test_long"])

ld1ld2$longitudinal <- ld1ld2$sub %in% longitudinal_ppts

# This includes the training scan, and the rest
results_longonly <- filter(ld1ld2, longitudinal)

ggplot(results_longonly, aes(x = days_since_first_scan / 30, y = LD1)) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(sub, hemi))) +
  geom_smooth(method = "lm", aes(linetype = hemi)) +
  theme_bw() +
  labs(x = "Months since first scan", y = "LD1") +
  facet_wrap(vars(PatientStatus)) 

ggsave("plots/longitudinal_ld1.png", width = 6.5, height = 4)

lm_randeff_time <- lmer(abs(LD1) ~ days_since_first_scan + (1 | sub), 
                        data = results_longonly)

summary(lm_randeff_time)

lm_randeff_time_LHS <- lmer(abs(LD1) ~ days_since_first_scan + (1 | sub), 
                            data = filter(results_longonly, class == "LHS_RH"))

summary(lm_randeff_time_LHS)

ggplot(results_longonly, aes(x = days_since_first_scan / 30, y = LD2, color = hemi)) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(sub, hemi))) +
  geom_smooth(method = "lm", aes(linetype = hemi)) +
  theme_bw() +
  labs(x = "Months since first scan", y = "LD1", color = "hemi") +
  facet_wrap(vars(PatientStatus))

ggsave("plots/longitudinal_ld2.png", width = 6.5, height = 4)

lm_randeff_time2 <- lmer(LD2 ~ days_since_first_scan + (1 | sub), 
                         data = results_longonly)

summary(lm_randeff_time2)

lm_randeff_time2_LHS <- lmer(LD2 ~ days_since_first_scan + (1 | sub), 
                             data = filter(results_longonly, class == "LHS_RH"))

summary(lm_randeff_time2_LHS)
