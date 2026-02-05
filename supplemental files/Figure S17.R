###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2026-01-27
# Last updated 2026-01-27

# Supplemental figure - county level scatterplots

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


# erie plots
dat_erie <- dat_county %>%
  filter(county == "Erie")

erie_pi_plot <- county_cor_plot_function(
  dataframe = dat_county,
  county = "Erie",
  var_name = "ntd_pi_county_3w",
  outcome_name = "case_incidence",
  cor_value = "0.72",
  pal_color = pal[6],
  pop_value = "946,000",
  outcome_text = "Case incidence",
  p_va = "<0.001",
  var_text = expression(pi[ww])
)

# erie h
cor(dat_erie$ntd_h_county_3w,
    dat_erie$case_incidence,
    method = "spearman")

erie_h_plot <- county_cor_plot_function(
  dataframe = dat_county,
  county = "Erie",
  var_name = "ntd_h_county_3w",
  outcome_name = "case_incidence",
  cor_value = "0.76",
  pal_color = pal[5],
  pop_value = "946,000",
  outcome_text = "Case incidence",
  p_va = "<0.001",
  var_text = expression("H"[ww])
)

# erie freyja
cor.test(dat_erie$n_variants_no_thresh_3w_mean,
         dat_erie$case_incidence,
         method = "spearman")

erie_freyja_plot <- county_cor_plot_function(
  dataframe = dat_county,
  county = "Erie",
  var_name = "n_variants_no_thresh_3w_mean",
  outcome_name = "case_incidence",
  cor_value = "0.84",
  pal_color = pal[2],
  pop_value = "946,000",
  outcome_text = "Case incidence",
  p_va = "<0.001",
  var_text = "Freyja variant count"
)

# erie conc plot
cor.test(dat_erie$mean_sars2_conc_county_3w,
         dat_erie$case_incidence,
         method = "spearman")

erie_conc_plot <- county_cor_plot_function(
  dataframe = dat_county,
  county = "Erie",
  var_name = "mean_sars2_conc_county_3w",
  outcome_name = "case_incidence",
  cor_value = "0.62",
  pal_color = pal[1],
  pop_value = "946,000",
  outcome_text = "Case incidence",
  p_va = "<0.001",
  var_text = "Concentration"
)

# first column for the panel
col_1 <- grid.arrange(
  erie_pi_plot,
  erie_h_plot,
  erie_freyja_plot,
  erie_conc_plot,
  nrow = 4,
  top=textGrob("A) Erie county, population 946,000", 
               gp=gpar(fontsize=15))
)

# small population county (Steuben)
# Steuben plots
dat_Steuben <- dat_county %>%
  filter(county == "Steuben")

# Steuben pi
cor.test(dat_Steuben$ntd_pi_county_3w,
         dat_Steuben$case_incidence,
         method = "spearman")

Steuben_pi_plot <- county_cor_plot_function(
  dataframe = dat_county,
  county = "Steuben",
  var_name = "ntd_pi_county_3w",
  outcome_name = "case_incidence",
  cor_value = "0.68",
  pal_color = pal[6],
  pop_value = "93,000",
  outcome_text = "Case incidence",
  p_va = "<0.001",
  var_text = expression(pi[ww])
)

# Steuben h
cor.test(dat_Steuben$ntd_h_county_3w,
         dat_Steuben$case_incidence,
         method = "spearman")

Steuben_h_plot <- county_cor_plot_function(
  dataframe = dat_county,
  county = "Steuben",
  var_name = "ntd_h_county_3w",
  outcome_name = "case_incidence",
  cor_value = "0.66",
  pal_color = pal[5],
  pop_value = "93,000",
  outcome_text = "Case incidence",
  p_va = "<0.001",
  var_text = expression("H"[ww])
)

# Steuben freyja
cor.test(dat_Steuben$n_variants_no_thresh_3w_mean,
         dat_Steuben$case_incidence,
         method = "spearman")

Steuben_freyja_plot <- county_cor_plot_function(
  dataframe = dat_county,
  county = "Steuben",
  var_name = "n_variants_no_thresh_3w_mean",
  outcome_name = "case_incidence",
  cor_value = "0.74",
  pal_color = pal[2],
  pop_value = "93,000",
  outcome_text = "Case incidence",
  p_va = "<0.001",
  var_text = "Freyja variant count"
)

# Steuben conc plot
cor.test(dat_Steuben$mean_sars2_conc_county_3w,
         dat_Steuben$case_incidence,
         method = "spearman")

Steuben_conc_plot <- county_cor_plot_function(
  dataframe = dat_county,
  county = "Steuben",
  var_name = "mean_sars2_conc_county_3w",
  outcome_name = "case_incidence",
  cor_value = "0.60",
  pal_color = pal[1],
  pop_value = "95,000",
  outcome_text = "Case incidence",
  p_va = "<0.001",
  var_text = "Concentration"
)

# first column for the panel
col_2 <- grid.arrange(
  Steuben_pi_plot,
  Steuben_h_plot,
  Steuben_freyja_plot,
  Steuben_conc_plot,
  nrow = 4,
  top=textGrob("B) Steuben county, population 95,000", 
               gp=gpar(fontsize=15))
)

# Onondaga plots
dat_Onondaga <- dat_county %>%
  filter(county == "Onondaga")

# Onondaga pi
cor.test(dat_Onondaga$ntd_pi_county_3w,
         dat_Onondaga$case_incidence,
         method = "spearman")

Onondaga_pi_plot <- county_cor_plot_function(
  dataframe = dat_county,
  county = "Onondaga",
  var_name = "ntd_pi_county_3w",
  outcome_name = "case_incidence",
  cor_value = "0.72",
  pal_color = pal[6],
  pop_value = "468,000",
  outcome_text = "Case incidence",
  p_va = "<0.001",
  var_text = expression(pi[ww])
)

# Onondaga h
cor.test(dat_Onondaga$ntd_h_county_3w,
         dat_Onondaga$case_incidence,
         method = "spearman")

Onondaga_h_plot <- county_cor_plot_function(
  dataframe = dat_county,
  county = "Onondaga",
  var_name = "ntd_h_county_3w",
  outcome_name = "case_incidence",
  cor_value = "0.72",
  pal_color = pal[5],
  pop_value = "468,000",
  outcome_text = "Case incidence",
  p_va = "<0.001",
  var_text = expression("H"[ww])
)

# Onondaga freyja
cor.test(dat_Onondaga$n_variants_no_thresh_3w_mean,
         dat_Onondaga$case_incidence,
         method = "spearman")

Onondaga_freyja_plot <- county_cor_plot_function(
  dataframe = dat_county,
  county = "Onondaga",
  var_name = "n_variants_no_thresh_3w_mean",
  outcome_name = "case_incidence",
  cor_value = "0.60",
  pal_color = pal[2],
  pop_value = "468,000",
  outcome_text = "Case incidence",
  p_va = "<0.001",
  var_text = "Freyja variant count"
)

# Onondaga conc plot
cor.test(dat_Onondaga$mean_sars2_conc_county_3w,
         dat_Onondaga$case_incidence,
         method = "spearman")

Onondaga_conc_plot <- county_cor_plot_function(
  dataframe = dat_county,
  county = "Onondaga",
  var_name = "mean_sars2_conc_county_3w",
  outcome_name = "case_incidence",
  cor_value = "0.67",
  pal_color = pal[1],
  pop_value = "468,000",
  outcome_text = "Case incidence",
  p_va = "<0.001",
  var_text = "Concentration"
)

# first column for the panel
col_3 <- grid.arrange(
  Onondaga_pi_plot,
  Onondaga_h_plot,
  Onondaga_freyja_plot,
  Onondaga_conc_plot,
  nrow = 4,
  top=textGrob("C) Onondaga county, population 468,000", 
               gp=gpar(fontsize=15))
)

# panel
grid.arrange(
  col_1,col_2,col_3, nrow = 1
)

# save
# save
png("Supplemental files/Figure S - unknown.png",
    units = "in",
    width = 12, height = 11,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
grid.arrange(
  col_1,col_2,col_3, nrow = 1
)
showtext_end()
dev.off()

