###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-26

# Supplemental Figure S9

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

# -----------------------
# Figure S9 - depth and infections
# -----------------------

p1 <- 
  ggplot(data = dat_state)+
  geom_bar(aes(x = week,
               y = case_incidence,
               fill = "grey60"),
           position = "dodge",
           stat = "identity")+
  geom_line(aes(x = week,
                y = depth_state_3w/10,
                color = "black"),
            lwd = 1)+
  theme_dth_1+
  scale_y_continuous(
    "COVID-19 cases per 100k", 
    sec.axis = sec_axis(~ .*10, name = "Mean depth")
  )+
  scale_fill_manual(values = c("grey60" = "grey60"),
                    name = "",
                    label = c("Cases"))+
  scale_color_manual(values = c("black" = "black"),
                     name = "",
                     label = "Depth")+
  theme(legend.position = "bottom",
        legend.background = element_blank())+
  labs(x = "")

# correlation plot
p2 <- 
  ggplot(data = dat_state,
         aes(x = depth_state_3w,
             y = case_incidence))+
  geom_point(size = 2, color = "black")+
  stat_cor()+
  theme_dth_1+
  labs(x = "Mean depth",
       y = "COVID-19 cases per 100k")

# save
png("Figures/Figure S9.png",
    units = "in",
    width = 8.5, height = 6,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
plot_grid(
  plot_grid(p1 + theme(legend.position = "none"), p2, nrow = 1,
            align = "h",
            axis = "bt",
            rel_widths = c(2,1.5),
            labels = c("A", "B")),
  g_legend(p1),
  nrow = 2,
  rel_heights = c(2,0.5))
showtext_end()
dev.off()