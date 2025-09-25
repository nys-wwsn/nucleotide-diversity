###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-25

# Main analysis - Figure 3

# --------------------------------------
# PACKAGES
# --------------------------------------

library(docstring)
library(showtext)
library(ggplot2)
library(sf)
library(zoo)
library(stringr)
library(lubridate)
library(flextable)
library(tidyr)
library(dplyr)
library(ggpubr)
library(gridExtra)
library(confintr)
library(MetBrewer)
library(broom)
library(cowplot)
library(performance)
library(rcompanion)

# Load functions
source("seq diversity - functions.R")

# Load plot themes
source("seq diversity - plot themes.R")

# -----------------------------------------------------------------------------

# --------------------------------------
# DATA LOAD
# --------------------------------------

# combined data files for each geography
load(file = "data/combined_data.Rdata")

# max and min weeks for the study
min_week <- "2023-01-01"
max_week <- "2025-04-20"

# -----------------------------------------------------------------------------


# ------------------------------------------------
# Figure 3 - correlation with increased aggregation
# ------------------------------------------------
# case incidence over time
case_plot <- 
  ggplot(data = dat_state)+
  geom_bar(aes(x = week,
               y = case_incidence,
               fill = "grey60"),
           position = "dodge",
           stat = "identity")+
  geom_line(
    aes(x = week,
        y = n_variants_no_thresh_3w_mean,
        color = "#e76254"),
    lwd = 1
  )+
  geom_line(
    aes(x = week,
        y = ntd_h_state_3w*10000,
        color = "#72bcd5"),
    lwd = 1
  )+
  geom_line(
    aes(x = week,
        y = ntd_pi_state_3w*10000,
        color = "#376795"),
    lwd = 1
  )+
  theme_dth_1+
  labs(
    title = "Case incidence (per 100,000)",
    x = "",
    y = "Cases"
  )+
  theme(legend.position = "bottom",
        legend.background = element_blank())+
  scale_color_manual(
    values = c("#e76254" ="#e76254",
               "#72bcd5" = "#72bcd5",
               "#376795" = "#376795"),
    labels = c(
      expression("S1 NTD "~ pi[ww]),
      expression("S1 NTD H"[ww]),
      "Freyja variant count"),
    
    name = ""
  )+
  scale_fill_manual(values = c("grey60" = "grey60"),
                    name = "",
                    label = "Cases")+
  scale_y_continuous(
    "COVID-19 cases per 100k", 
    sec.axis = sec_axis(~ ., name = expression(paste("Variant count/",
                                                     ~Pi[ww], "/", "H"[ww], ))
    )
  )
case_plot

# hosp incidence
hosp_plot <- 
  ggplot(data = dat_state)+
  geom_bar(aes(x = week,
               y = hosp_incidence,
               fill = "grey70"),
           position = "dodge",
           stat = "identity")+
  geom_bar(aes(fill = "grey60",
               x = week,
               y = 0),
           position = "dodge",
           stat = "identity")+
  geom_line(
    aes(x = week,
        y = n_variants_no_thresh_3w_mean/8,
        color = "#e76254"),
    lwd = 1
  )+
  geom_line(
    aes(x = week,
        y = ntd_h_state_3w*1000,
        color = "#72bcd5"),
    lwd = 1
  )+
  geom_line(
    aes(x = week,
        y = ntd_pi_state_3w*1000,
        color = "#376795"),
    lwd = 1
  )+
  theme_dth_1+
  labs(
    title = "Hospitalizations (per 100,000)",
    x = ""
  )+
  theme(legend.position = "bottom",
        legend.background = element_blank())+
  scale_color_manual(
    values = c("#e76254" ="#e76254",
               "#72bcd5" = "#72bcd5",
               "#376795" = "#376795"),
    labels = c(
      expression("S1 NTD "~ pi[ww]),
      expression("S1 NTD H"[ww]),
      "Freyja variant count"),
    
    name = ""
  )+
  scale_fill_manual(values = c("grey60" = "grey60"),
                    name = "",
                    label = "Cases/Hospitalizations")+
  scale_y_continuous(
    "COVID-19 hospitalizations per 100k", 
    sec.axis = sec_axis(~ .*8, name = expression(paste("Variant count/",
                                                       ~Pi[ww], "/", "H"[ww], ))
    )
  )
hosp_plot

mylegend <- g_legend(hosp_plot)

plot_grid(case_plot + theme(legend.position = "none"),
          hosp_plot + theme(legend.position = "none"),
          mylegend,
          nrow = 3,
          align = "v",
          axis = "l",
          labels = c("A", "B"),
          rel_heights = c(3,3,1))

# save
png("Figures/Figure 3.png",
    units = "in",
    width = 6, height = 8,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
plot_grid(case_plot + theme(legend.position = "none"),
          hosp_plot + theme(legend.position = "none"),
          mylegend,
          nrow = 3,
          align = "v",
          axis = "l",
          labels = c("A", "B"),
          rel_heights = c(3,3,1))
showtext_end()
dev.off()