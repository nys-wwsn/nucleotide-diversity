###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2026-02-03
# Last updated 2026-02-03

# Supplemental Figure S4

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

# --------------------------------------------------------------------
# Figure S4 - county level ntd concentration and case correlation
# --------------------------------------------------------------------

# ntd pi over time by county
ntd_pi_plot <-
  ggplot(data = dat_sewershed)+
  geom_point(aes(x = week,
                 y= ntd_pi_ma3),
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
  ggplot(data = dat_sewershed)+
  geom_point(aes(x = week,
                 y= ntd_h_ma3),
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
  ggplot(data = dat_sewershed)+
  geom_point(aes(x = week,
                 y= n_lineages_no_thresh_3w,
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
png("Supplemental files/_Figure S4.png",
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
