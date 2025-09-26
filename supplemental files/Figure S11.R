###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-26

# Supplemental Figure S11

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

# ------------------------------------
# Figure S11 - random subsampling
# ------------------------------------


# Randomly subsample each sample to equal read depth and reassess correlation

# relationship prior to subsample
ggplot(data = dat_sewershed %>%
         filter(region != "New York City"))+
  geom_point(aes(x = week, y = ntd_pi, color = facility_id))+
  theme_dth_1+
  theme(legend.position = "none")+
  geom_smooth(aes(x = week, y = ntd_pi),
              method = "loess",
              span = 0.1,
              color = "darkblue",
              lwd = 1.5)+
  labs(x = "",
       y = "S1 NTD Pi wastewater",
       title = "All samples")+
  scale_color_manual(values = MetBrewer::met.brewer(name = "Austria",
                                                    n = 204))


# subsample for depth reading between 100 and 500
set.seed(23)

# Depth subsample - 3 week moving average for state values
# statewide average - generate example for the figure legend

cases <- dat_state %>%
  dplyr::select(week, case_incidence)

# sample 3
sample_3_3w <- dat_sewershed  %>%
  filter(depth >= 500 & depth <= 1000) %>%
  sample_n(1000) %>%
  arrange(week) %>%
  group_by(week)%>%
  summarise(
    # statewide pi
    genomewide_pi_state_3w = weighted.mean(x = genomewide_pi_ma3, 
                                           w = population_served, 
                                           na.rm = TRUE),
    ntd_pi_state_3w = weighted.mean(x = ntd_pi_ma3, 
                                    w = population_served, 
                                    na.rm = TRUE)
  )%>%
  ungroup() %>%
  left_join(cases, by = c("week"))

# generate legend
p_time <- 
  ggplot(data = sample_3_3w)+
  geom_bar(data = dat_state %>%
             filter(group == "State"),
           aes(x = week,
               y = case_incidence,
               fill = "grey60"),
           position = "dodge",
           stat = "identity")+
  geom_line(aes(x = week,
                y = ntd_pi_state_3w*10000,
                color = "darkblue"),
            lwd = 1.5)+
  theme_dth_1+
  scale_color_manual(values = c("darkblue"),
                     labels = c("Pi"),
                     name = "")+
  scale_fill_manual(values = c("grey60"),
                    labels = "Cases",
                    name = "")+
  theme(legend.background = element_blank())+
  theme(legend.position = "bottom")

# extract legend function
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# extract legend
mylegend<-g_legend(p_time)



sample_plots_1 <- time_cor_plot_function(
  dataframe = dat_sewershed,
  depth_start = 100,
  depth_end = 400,
  title = "Depth read of sample 100 to 400\n1000 randomly selected samples"
)

sample_plots_2 <- time_cor_plot_function(
  dataframe = dat_sewershed,
  depth_start = 400,
  depth_end = 800,
  title = "Depth read of sample 400 to 800\n1000 randomly selected samples")

sample_plots_3 <- time_cor_plot_function(
  dataframe = dat_sewershed,
  depth_start = 800,
  depth_end = 1500,
  title = "Depth read of sample > 800\n1000 randomly selected samples")

sample_plots_original <- time_cor_plot_function(
  dataframe = dat_sewershed,
  depth_start = 0,
  depth_end = 45000,
  title = "Depth read of sample (all reads)\n1000 randomly selected samples")

ggarrange(ggarrange(sample_plots_original,
                    sample_plots_1,
                    sample_plots_2,
                    sample_plots_3,
                    nrow = 4),mylegend, nrow=2,heights=c(10, 1))


png("Supplemental files/Figure S11.png",
    units = "in",
    width = 8, height = 12,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)

ggarrange(ggarrange(sample_plots_original,
                    sample_plots_1,
                    sample_plots_2,
                    sample_plots_3,
                    nrow = 4),mylegend, nrow=2,heights=c(10, 1))
showtext_end()
dev.off()