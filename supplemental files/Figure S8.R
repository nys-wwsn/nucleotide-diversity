###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-26

# Supplemental Figure S13

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
library(lmtest)

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
  dplyr::select(date, cdc_id, mean_depth)

# mean depth post july 2024
depth_2 <- read.csv("data/other data/genomwide_depth_2025-07-02.csv") %>%
  mutate(date = ymd(date)) %>%
  dplyr::select(date, cdc_id, mean_depth) %>%
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

# ---------------------------------------------
# Figure S13
# ---------------------------------------------
# Depth from freyja variants (alteration detected and depth recorded)
# compared to depth from the total sample
depth_plot_7 <- 
  ggplot(data = dat_sewershed ,
         aes(x = depth,
             y = mean_depth)
  )+
  geom_point(alpha = 0.5, color = "seagreen")+
  stat_cor(method = "spearman",
           digits = 4)+
  theme_dth_1+
  labs(x = "Depth - alterations only",
       y = "Depth - complete sample")

depth_plot_7

# save
png("Supplemental files/Figure S13.png",
    units = "in",
    width = 5, height = 5,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
depth_plot_7
showtext_end()
dev.off()
