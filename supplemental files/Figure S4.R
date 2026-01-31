###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2026-01-27
# Last updated 2026-01-27

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
library(grid)

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


# -----------------
# Figure
# -----------------

# regional boxplot
region_boxplot <- spatial_boxplot_function(dataframe = dat_region,
                                           aggregation_name = "region",
                                           outcome_var = "case_incidence",
                                           aggregation_label = "Regional",
                                           pi_name = "ntd_pi_region_3w",
                                           h_name = "ntd_h_region_3w",
                                           conc_name = "mean_sars2_conc_region_3w",
                                           freyja_name = "n_variants_no_thresh_3w_mean")
region_boxplot


# county
county_boxplot <- spatial_boxplot_function(dataframe = dat_county,
                                           aggregation_name = "county",
                                           outcome_var = "case_incidence",
                                           aggregation_label = "County",
                                           pi_name = "ntd_pi_county_3w",
                                           h_name = "ntd_h_county_3w",
                                           conc_name = "mean_sars2_conc_county_3w",
                                           freyja_name = "n_variants_no_thresh_3w_mean")
county_boxplot

# sewershed
county_cases <- dat_county %>%
  select(county, case_incidence, hosp_incidence, week) %>%
  distinct()
dat_sewershed <- left_join(dat_sewershed,
                           county_cases,
                           by = c("county", "week"))

sewershed_boxplot <- spatial_boxplot_function(
  dataframe = dat_sewershed,
  aggregation_name = "facility_id",
  outcome_var = "case_incidence",
  aggregation_label = "Sewershed",
  pi_name = "ntd_pi_ma3",
  h_name = "ntd_h_ma3",
  conc_name = "mean_sars2_conc_ma3",
  freyja_name = "n_variants_no_thresh_3w"
)

sewershed_boxplot

# get legend
legend <- g_legend(sewershed_boxplot)

case_row <- plot_grid(
  sewershed_boxplot + theme(legend.position = "none"),
  county_boxplot+theme(legend.position = "none"),
  region_boxplot+theme(legend.position = "none"),
  labels = c("A", "B", "C"),
  rel_widths = c(1, 1,1),
  nrow = 1
)

# title
title <- ggdraw() + draw_label("Correlation with case incidence", 
                               fontface='bold')

case_row <- plot_grid(
  title,
  case_row,
  nrow = 2,
  rel_heights = c(0.1,1)
)

# hosp row

# regional boxplot
region_boxplot <- spatial_boxplot_function(dataframe = dat_region,
                                           aggregation_name = "region",
                                           outcome_var = "hosp_incidence",
                                           aggregation_label = "Regional",
                                           pi_name = "ntd_pi_region_3w",
                                           h_name = "ntd_h_region_3w",
                                           conc_name = "mean_sars2_conc_region_3w",
                                           freyja_name = "n_variants_no_thresh_3w_mean")
region_boxplot


# county
county_boxplot <- spatial_boxplot_function(dataframe = dat_county,
                                           aggregation_name = "county",
                                           outcome_var = "hosp_incidence",
                                           aggregation_label = "County",
                                           pi_name = "ntd_pi_county_3w",
                                           h_name = "ntd_h_county_3w",
                                           conc_name = "mean_sars2_conc_county_3w",
                                           freyja_name = "n_variants_no_thresh_3w_mean")
county_boxplot

# sewershed
sewershed_boxplot <- spatial_boxplot_function(
  dataframe = dat_sewershed,
  aggregation_name = "facility_id",
  outcome_var = "hosp_incidence",
  aggregation_label = "Sewershed",
  pi_name = "ntd_pi_ma3",
  h_name = "ntd_h_ma3",
  conc_name = "mean_sars2_conc_ma3",
  freyja_name = "n_variants_no_thresh_3w"
)

sewershed_boxplot

# get legend
legend <- g_legend(sewershed_boxplot)

# title
title <- ggdraw() + draw_label("Correlation with hospital incidence", 
                               fontface='bold')


hosp_row <- plot_grid(
  sewershed_boxplot + theme(legend.position = "none"),
  county_boxplot+theme(legend.position = "none"),
  region_boxplot+theme(legend.position = "none"),
  labels = c("D", "E", "F"),
  rel_widths = c(1, 1,1),
  nrow = 1
)

hosp_row <- plot_grid(
  title,
  hosp_row,
  rel_heights = c(0.1, 1),
  nrow = 2
)


# legend for figure
s_plot <- sewershed_boxplot+
  theme(legend.position = "bottom")
legend <- g_legend(s_plot)

# panel

plot_grid(
  case_row,
  hosp_row,
  legend,
  nrow = 3,
  rel_heights = c(2,2,0.5)
)

# save
png("Figures/Figure 3.png",
    units = "in",
    width = 10, height = 8,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
plot_grid(
  case_row,
  hosp_row,
  legend,
  nrow = 3,
  rel_heights = c(2,2,0.5)
)
showtext_end()
dev.off()

