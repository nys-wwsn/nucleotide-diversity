###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2026-01-27
# Last updated 2026-01-27

# Main analysis - table 1 with comparison between spatial scales and 
# spearman correlation coefficients

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


# sewershed
county_cases <- dat_county %>%
  select(county, case_incidence, hosp_incidence, week) %>%
  distinct()
dat_sewershed <- left_join(dat_sewershed,
                           county_cases,
                           by = c("county", "week"))
# -------------------------
# sewershed  case correlation comparison
# -------------------------
# case
sewershed_case_results <- spatial_cor_compare_function(
  dataframe = dat_sewershed,
  aggregation ="Sewershed",
  outcome_var = "case_incidence",
  pi_var = "ntd_pi_ma3",
  h_var = "ntd_h_ma3",
  conc_var = "mean_sars2_conc_ma3",
  freyja_var = "n_variants_no_thresh_3w"
)

# hosp
sewershed_hosp_results <- spatial_cor_compare_function(
  dataframe = dat_sewershed,
  aggregation ="Sewershed",
  outcome_var = "hosp_incidence",
  pi_var = "ntd_pi_ma3",
  h_var = "ntd_h_ma3",
  conc_var = "mean_sars2_conc_ma3",
  freyja_var = "n_variants_no_thresh_3w"
)
# -------------------------
# county comparison
# -------------------------

# case
county_case_results <- spatial_cor_compare_function(
  dataframe = dat_county,
  aggregation ="County",
  outcome_var = "case_incidence",
  pi_var = "ntd_pi_county_3w",
  h_var = "ntd_h_county_3w",
  conc_var = "mean_sars2_conc_county_3w",
  freyja_var = "n_variants_no_thresh_3w_mean"
)

# hosp
county_hosp_results <- spatial_cor_compare_function(
  dataframe = dat_county,
  aggregation ="County",
  outcome_var = "hosp_incidence",
  pi_var = "ntd_pi_county_3w",
  h_var = "ntd_h_county_3w",
  conc_var = "mean_sars2_conc_county_3w",
  freyja_var = "n_variants_no_thresh_3w_mean"
)

# -------------------------
# regional comparison
# -------------------------

# case
region_case_results <- spatial_cor_compare_function(
  dataframe = dat_region,
  aggregation ="Region",
  outcome_var = "case_incidence",
  pi_var = "ntd_pi_region_3w",
  h_var = "ntd_h_region_3w",
  conc_var = "mean_sars2_conc_region_3w",
  freyja_var = "n_variants_no_thresh_3w_mean"
)

# hosp
region_hosp_results <- spatial_cor_compare_function(
  dataframe = dat_region,
  aggregation ="Region",
  outcome_var = "hosp_incidence",
  pi_var = "ntd_pi_region_3w",
  h_var = "ntd_h_region_3w",
  conc_var = "mean_sars2_conc_region_3w",
  freyja_var = "n_variants_no_thresh_3w_mean"
)

# -------------------------
# state comparison
# -------------------------

# case
state_case_results <- spatial_cor_compare_function(
  dataframe = dat_state,
  aggregation ="State",
  outcome_var = "case_incidence",
  pi_var = "ntd_pi_state_3w",
  h_var = "ntd_h_state_3w",
  conc_var = "mean_sars2_conc_state_3w",
  freyja_var = "n_variants_no_thresh_3w_mean"
)

# hosp
state_hosp_results <- spatial_cor_compare_function(
  dataframe = dat_state,
  aggregation ="State",
  outcome_var = "hosp_incidence",
  pi_var = "ntd_pi_state_3w",
  h_var = "ntd_h_state_3w",
  conc_var = "mean_sars2_conc_state_3w",
  freyja_var = "n_variants_no_thresh_3w_mean"
)


# -------------------------
# combine dfs
# -------------------------

table_df <- bind_rows(
  sewershed_case_results,
  sewershed_hosp_results,
  county_case_results,
  county_hosp_results,
  region_case_results,
  region_hosp_results,
  state_case_results,
  state_hosp_results
)

# save
table_ft <- table_as_flex_function(dataframe = table_df,
                                   title = "Spearman correlations across spatial scales")

# save
# save
save_as_docx(FitFlextableToPage(table_ft), 
             path = paste("Figures/",
                          "Table spearman comparison.docx" ,sep = ""))
