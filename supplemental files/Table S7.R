###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-26

# Supplemental Table S6

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


# -------------------------------------------------------------------------
# Table S6 - Granger causality results
# -------------------------------------------------------------------------

# Granger causality tests for whether time series x predicts time series y
# Both time series are adjusted to be stationary.

# Remove NAs from the data and arrange the time series
s <- dat_state %>%
  filter(!is.na(case_incidence)) %>%
  filter(!is.na(ntd_pi_state_3w)) %>%
  filter(!is.na(mean_sars2_conc_state_3w)) %>%
  filter(!is.na(depth_state_3w)) %>%
  filter(!is.na(n_variants_no_thresh_3w_mean)) %>%
  arrange(week) 

# granger test for s1 ntd pi and cases
g_ntd_pi_cases <- 
  granger_casuality_function(
    dataframe = s,
    y = "case_incidence",
    x = "ntd_pi_state_3w",
    lag = 1,
    group = "S1 NTD Pi",
    outcome = "Case incidence"
  )

# granger test for s1 ntd h and cases
g_ntd_h_cases <- 
  granger_casuality_function(
    dataframe = s,
    y = "case_incidence",
    x = "ntd_h_state_3w",
    lag = 1,
    group = "S1 NTD H",
    outcome = "Case incidence"
  )

# granger test for freyja variant count and  cases
g_var_count_cases <- 
  granger_casuality_function(
    dataframe = s,
    y = "case_incidence",
    x = "n_variants_no_thresh_3w_mean",
    lag = 1,
    group = "Freyja variant count",
    outcome = "Case incidence"
  )

# granger test for s1 ntd pi and hospitalizations
g_ntd_pi_hosp <- 
  granger_casuality_function(
    dataframe = s,
    y = "hosp_incidence",
    x = "ntd_pi_state_3w",
    lag = 1,
    group = "S1 NTD Pi",
    outcome = "Hospitalization incidence"
  )

# granger test for s1 ntd h and hospitalizations
g_ntd_h_hosp <- 
  granger_casuality_function(
    dataframe = s,
    y = "hosp_incidence",
    x = "ntd_h_state_3w",
    lag = 1,
    group = "S1 NTD H",
    outcome = "Hospitalization incidence"
  )

# granger test for freyja variant count and hospitalizations
g_var_count_hosp <- 
  granger_casuality_function(
    dataframe = s,
    y = "hosp_incidence",
    x = "n_variants_no_thresh_3w_mean",
    lag = 1,
    group = "Freyja variant count",
    outcome = "Hospitalization incidence"
  )

# combine into one dataframe
granger_results <- bind_rows(
  g_ntd_pi_cases,
  g_ntd_h_cases,
  g_var_count_cases,
  g_ntd_pi_hosp,
  g_ntd_h_hosp,
  g_var_count_hosp
)

# round results
granger_results <- granger_results %>%
  mutate_if(., is.numeric, round_signifi_function)

# label p values
granger_results <- granger_results %>%
  mutate(`P value` = case_when(
    `P value` > 0.05 ~ as.character(`P value`),
    `P value` < 0.05 & `P value` > 0.01~ paste(as.character(`P value`), "*", 
                                               sep = ""),
    `P value` < 0.01 & `P value` > 0.001 ~ paste(as.character(`P value`), "**", 
                                                 sep = ""),
    `P value` < 0.001 ~ paste("<0.001***")
  ))

# edit column names
colnames(granger_results) <- c("Diversity measure (x)",
                               "Outcome variable (y)",
                               "Model",
                               "F statistic",
                               "P value")

# make it a flextable object
granger_ft <- table_as_flex_function(granger_results,
                                     title = "Table S6: Granger causality results")

granger_ft

# save
save_as_docx(FitFlextableToPage(granger_ft), 
             path = paste("Figures/",
                          "Table S6.docx",sep = ""))
