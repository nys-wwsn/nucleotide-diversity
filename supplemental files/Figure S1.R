###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-25

# Supplemental Figure S1

# ---------------------------------------
# Packages
# ---------------------------------------

library(dplyr)
library(ggplot2)
library(ggpubr)
library(showtext)
library(lme4)
library(lmerTest)
library(performance)
library(lubridate)
library(gridExtra)
library(cowplot)
library(flextable)
library(performance)
library(rcompanion)
library(lme4)
library(lmerTest)
library(nlme)
library(MASS)
library(glmmTMB)

# --------------------------------------
# PLOT THEMES
# --------------------------------------
source("seq diversity - plot themes.R")

# ---------------------------------------
# FUNCTIONS
# ---------------------------------------
source("seq diversity - functions.R")

# ---------------------------------------
# Load data
# ---------------------------------------

# combined data files for each geography
load(file = "data/combined_data.Rdata")

# max and min weeks for the study
min_week <- "2023-01-01"
max_week <- "2025-04-20"

# mean depth per sample - data prior to July 2024
depth_1 <- read.csv("data/other data/depth_stats_20250424.txt") %>%
  tidyr::separate(
    sample, into = c("date", "cdc_id"), sep = 8
  ) %>%
  mutate(date = lubridate::ymd(date),
         mean_depth = mean) %>%
  select(date, cdc_id, mean_depth)

# mean depth post july 2024
depth_2 <- read.csv("data/other data/genomwide_depth_2025-07-02.csv") %>%
  mutate(date = ymd(date)) %>%
  select(date, cdc_id, mean_depth) %>%
  tidyr::separate(
    cdc_id, into = c("cdc_id", "drop"), sep = 12
  )

# bind together, then merge to the dat_sewershed object
depth <- bind_rows(depth_1, depth_2) %>%
  group_by(cdc_id, week = floor_date(date, unit = "weeks", week_start = 7)) %>%
  summarize(mean_depth = mean(mean_depth, na.rm =  TRUE)
  ) %>%
  ungroup() %>%
  rename(facility_id = cdc_id)

dat_sewershed <- left_join(dat_sewershed, depth, by = c("facility_id", "week"))

# ---------------------------------------
# Supplemental Figure S1 - depth and ct
# ---------------------------------------

# Depth with conc and Depth with ct - 4 figure panel with correlations and 
# time series.

depth_ct <- 
  ggplot(data = dat_sewershed,
         aes(x = mean_depth,
             y = mean_ct))+
  geom_point(alpha = 0.5, color = "tomato")+
  stat_cor(method = "spearman")+
  theme_dth_1+
  labs(x = "Mean depth",
       y = "Mean Ct",
       title = "Ct ~ depth")
depth_ct

# ct over time

ct_time <- 
  ggplot(data = dat_sewershed,
         aes(x = week,
             y = mean_ct))+
  geom_point(alpha = 0.5, color = "tomato")+
  theme_dth_1+
  labs(x = "",
       y = "Mean Ct",
       title = "Ct over time")
ct_time

# depth by concentration
depth_conc <- 
  ggplot(data = dat_sewershed,
         aes(x = mean_depth,
             y = mean_sars2_conc))+
  geom_point(alpha = 0.5, color = "orange")+
  stat_cor(method = "spearman")+
  theme_dth_1+
  labs(x = "Mean depth",
       y = "Conc, copies/mL",
       title = "Concentration ~ depth")
depth_conc

# ct over time
conc_time <- 
  ggplot(data = dat_sewershed,
         aes(x = week,
             y = mean_sars2_conc))+
  geom_point(alpha = 0.5, color = "orange")+
  theme_dth_1+
  labs(x = "",
       y = expression(paste("Conc, copies/mL",
       )),
       title = "Concentration over time")
conc_time

# coverage by depth
coverage_plot <- 
  ggplot(data = dat_sewershed,
         aes(x = mean_depth,
             y = coverage
         ))+
  geom_point(alpha = 0.5, color = "dodgerblue")+
  stat_cor(method = "spearman")+
  theme_dth_1 +
  labs(x = "Mean depth",
       y = "Mean coverage (%)",
       title = "Coverage ~ depth")

# coverage over time
coverage_time <- 
  ggplot(data = dat_sewershed,
         aes(x = week,
             y = coverage))+
  geom_point(alpha = 0.5, color = "dodgerblue")+
  theme_dth_1+
  labs(x = "",
       y = "Mean coverage (%)",
       title = "Coverage over time")
coverage_time

# panel figure
plot_grid(depth_ct, depth_conc,coverage_plot,
          ct_time, conc_time,
          coverage_time,
          nrow = 2,
          labels = c("A", "B", "C", "D", "E", "F"))

# save
png("Supplemental files/Figure S1.png",
    units = "in",
    width = 12, height = 6,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
plot_grid(depth_ct, depth_conc,coverage_plot,
          ct_time, conc_time,
          coverage_time,
          nrow = 2,
          labels = c("A", "B", "C", "D", "E", "F"))
showtext_end()
dev.off()