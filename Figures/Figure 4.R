###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-25

# Main analysis - Figure 5

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

# --------------------------------------------------------------------------
# Figure 5 - lag and lead times for all measures
# --------------------------------------------------------------------------

# case
lag_case <- lag_lead_function(dataframe = dat_state,
                              columns = c(# pi
                                "ntd_pi_state_3w",
                                # H
                                "ntd_h_state_3w",
                                # var count
                                "n_variants_no_thresh_3w_mean",
                                # conc
                                "mean_sars2_conc_state_3w"),
                              lags = seq(1:10))

lag_case_plot <- ggplot(data = lag_case ,
                        aes(x = lag, y = spearman_cor, color = name))+
  geom_point(shape = 1, size = 2)+
  geom_point(shape = 17,
             data = lag_case %>%
               group_by(name) %>%
               filter(spearman_cor == max(spearman_cor)),
             size = 3
  )+
  geom_line()+
  ylim(0,1)+
  theme_dth_1+
  theme(legend.position = "bottom",
        panel.grid.major.x = element_line(color = "gray80"),
        panel.grid.minor.y = element_blank())+
  scale_color_manual(values = c(
    "ntd_pi_state_3w"=pal[6],
    "ntd_h_state_3w"=pal[5],
    "n_variants_no_thresh_3w_mean"=pal[1],
    "mean_sars2_conc_state_3w"=pal[2]
  ),
  name = "",
  labels = c(
    "ntd_pi_state_3w"=expression("S1 NTD "~ pi[ww]),
    "ntd_h_state_3w"=expression( "S1 NTD H"[ww]),
    "n_variants_no_thresh_3w_mean"="Freyja Variant counts",
    "mean_sars2_conc_state_3w"="Concentration"
  )
  )+
  scale_x_continuous(breaks = seq(-5, 5, by = 1),
                     limits = c(-5,5))+
  labs(x = "Lead time (weeks)",
       y = "spearman correlation\ncoefficient",
       title= "Case incidence ~ Wastewater measure for\nstatewide data")+
  annotate("text",
           x = -3, y = 0.09,
           label = "signal lags incidence")+
  annotate("text",
           x = 3, y = 0.09,
           label = "signal leads incidence")+
  guides(color=guide_legend(nrow=2,byrow=TRUE))
lag_case_plot

# hosp

lag_hosp <- lag_lead_function_hosp(dataframe = dat_state,
                                   columns = c(# pi
                                     "ntd_pi_state_3w",
                                     # H
                                     "ntd_h_state_3w",
                                     # var count
                                     "n_variants_no_thresh_3w_mean",
                                     # conc
                                     "mean_sars2_conc_state_3w"),
                                   lags = seq(1:10))

lag_hosp_plot <- ggplot(data = lag_hosp ,
                        aes(x = lag, y = spearman_cor, color = name))+
  geom_point(shape = 1, size = 2)+
  geom_point(shape = 17,
             data = lag_hosp %>%
               group_by(name) %>%
               filter(spearman_cor == max(spearman_cor)),
             size = 3
  )+
  geom_line()+
  ylim(0,1)+
  theme_dth_1+
  theme(legend.position = "bottom",
        panel.grid.major.x = element_line(color = "gray80"),
        panel.grid.minor.y = element_blank())+
  scale_color_manual(values = c(
    "ntd_pi_state_3w"=pal[6],
    "ntd_h_state_3w"=pal[5],
    "n_variants_no_thresh_3w_mean"=pal[1],
    "mean_sars2_conc_state_3w"=pal[2]
  ),
  name = "",
  labels = c(
    "ntd_pi_state_3w"=expression("S1 NTD "~ pi[ww]),
    "ntd_h_state_3w"=expression( "S1 NTD H"[ww]),
    "n_variants_no_thresh_3w_mean"="Freyja Variant counts",
    "mean_sars2_conc_state_3w"="Concentration"
  )
  )+
  scale_x_continuous(breaks = seq(-5, 5, by = 1),
                     limits = c(-5,5))+
  labs(x = "Lead time (weeks)",
       y = "spearman correlation\ncoefficient",
       title= "Hospitalization incidence ~ Wastewater measure for\nstatewide data")+
  annotate("text",
           x = -3, y = 0.09,
           label = "signal lags incidence")+
  annotate("text",
           x = 3, y = 0.09,
           label = "signal leads incidence")


# extract legend function
mylegend <- g_legend(lag_case_plot)


# save
png("Figures/Figure 5.png",
    units = "in",
    width = 6.5, height = 10,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
plot_grid(
  plot_grid(lag_case_plot+ theme(legend.position = "none"),
            lag_hosp_plot + theme(legend.position = "none"),
            labels = c("A", "B"), nrow = 2, ncol = 1),
  plot_grid(NULL, mylegend, NULL, nrow = 1, rel_widths = c(1, 0.5, 1)),
  nrow = 2,
  rel_heights = c(5,1)
)
dev.off()


