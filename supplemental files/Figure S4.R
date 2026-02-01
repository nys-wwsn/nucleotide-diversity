###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-25

# Supplemental Figure S2

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

# --------------------------------------------------------------------
# Figure S2 - county level ntd concentration and case correlation
# --------------------------------------------------------------------

# ntd pi over time by county
ntd_pi_plot <-
  ggplot(data = dat_county)+
  geom_point(aes(x = week,
                 y= ntd_pi_county_3w),
             alpha = 0.6,
             color = "grey60")+
  theme_dth_1+
  labs(
    title = expression(paste( pi[ww])),
    x = "",
    y = expression(paste( pi[ww]))
  )

ntd_pi_plot

# ntd h over time by county
ntd_h_plot <-
  ggplot(data = dat_county)+
  geom_point(aes(x = week,
                 y= ntd_h_county_3w),
             alpha = 0.6,
             color = "grey60")+
  theme_dth_1+
  labs(
    title = expression(paste("H"[ww])),
    x = "",
    y = expression(paste("H"[ww]))
  )

# variant count over time by county
var_count_plot <- 
  ggplot(data = dat_county)+
  geom_point(aes(x = week,
                 y= n_variants_no_thresh_3w_mean,
                 fill = "grey60"),
             alpha = 0.6,
             color = "grey60")+
  theme_dth_1+
  labs(
    title = "Freyja variant counts",
    x = "",
    y = "Variant counts"
  )+
  scale_fill_manual(values = c("grey60" = "grey60"),
                    name = "",
                    labels = c("County weighted average"))+
  theme(legend.position = "bottom",
        legend.background = element_blank())
var_count_plot

# case incidence over time
case_plot <- 
  ggplot(data = dat_state)+
  geom_bar(aes(x = week,
               y = case_incidence),
           position = "dodge",
           stat = "identity")+
  theme_dth_1+
  labs(
    title = "Case incidence (per 100,000)",
    x = "",
    y = "Cases"
  )
mylegend <- g_legend(var_count_plot)


# put the plots in a panel
plot_grid(ntd_pi_plot, 
          ntd_h_plot, 
          var_count_plot+ 
            theme(legend.position = "none"),
          case_plot, 
          mylegend,
          nrow = 5,
          align = "v",
          axis = "l",
          labels = c("A", "B", "C", "D"),
          rel_heights = c(3,3,3,3,1))

# save
png("Supplemental files/Figure S2.png",
    units = "in",
    width = 8.5, height = 11,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
plot_grid(ntd_pi_plot, 
          ntd_h_plot, 
          var_count_plot+ 
            theme(legend.position = "none"),
          case_plot, 
          mylegend,
          nrow = 5,
          align = "v",
          axis = "l",
          labels = c("A", "B", "C", "D"),
          rel_heights = c(3,3,3,3,1))
showtext_end()
dev.off()
