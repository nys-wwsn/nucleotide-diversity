###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-26

# Supplemental Table S10

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
# -------------------------------------------------------------------
# Table S10 for case + hosp data
# -------------------------------------------------------------------

hosp_summary <- dat_county %>%
  ungroup() %>%
  filter(week >= "2023-01-01" & week <= "2025-04-20") %>%
  summarize(`Number of weeks with data` = length(unique(week)),
            `Mean` = mean(hospitalizations, na.rm = TRUE),
            `Min` = min(hospitalizations, na.rm = TRUE),
            `Median` = median(hospitalizations, na.rm = TRUE),
            `Max` = max(hospitalizations, na.rm = TRUE),
            `SD` = sd(hospitalizations, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(group = "Hospitalizations") %>%
  mutate_if(is.numeric, ~round(., 2))

# incidence
hosp_incidence_summary <- dat_county %>%
  ungroup() %>%
  filter(week >= "2023-01-01" & week <= "2025-04-20") %>%
  summarize(`Number of weeks with data` = length(unique(week)),
            `Mean` = mean(hosp_incidence, na.rm = TRUE),
            `Min` = min(hosp_incidence, na.rm = TRUE),
            `Median` = median(hosp_incidence, na.rm = TRUE),
            `Max` = max(hosp_incidence, na.rm = TRUE),
            `SD` = sd(hosp_incidence, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(group = "Hospital incidence (per 100k)") %>%
  mutate_if(is.numeric, ~round(., 2))

# cases
cases_summary <- dat_county %>%
  ungroup() %>%
  filter(week >= "2023-01-01" & week <= "2025-04-20") %>%
  summarize(`Number of weeks with data` = length(unique(week)),
            `Mean` = mean(cases, na.rm = TRUE),
            `Min` = min(cases, na.rm = TRUE),
            `Median` = median(cases, na.rm = TRUE),
            `Max` = max(cases, na.rm = TRUE),
            `SD` = sd(cases, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(group = "COVID-19 cases") %>%
  mutate_if(is.numeric, ~round(., 2))

# case incidence
case_incidence_summary <- dat_county %>%
  ungroup() %>%
  filter(week >= "2023-01-01" & week <= "2025-04-20") %>%
  summarize(`Number of weeks with data` = length(unique(week)),
            `Mean` = mean(case_incidence, na.rm = TRUE),
            `Min` = min(case_incidence, na.rm = TRUE),
            `Median` = median(case_incidence, na.rm = TRUE),
            `Max` = max(case_incidence, na.rm = TRUE),
            `SD` = sd(case_incidence, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(group = "COVID-19 case incidence (per 100k)") %>%
  mutate_if(is.numeric, ~round(., 2))

# combine
clinical_df_summary <- bind_rows(
  cases_summary,
  case_incidence_summary,
  hosp_summary,
  hosp_incidence_summary
) %>%
  dplyr::select(group, everything())

# make it an ft table
# make it flextable
clinical_ft <- table_as_flex_function(
  dataframe = clinical_df_summary,
  title = "Table S10: Descriptive statistics for COVID-19 clinical data")

# save
save_as_docx(FitFlextableToPage(clinical_ft), 
             path = paste0("Supplemental files/",
                           "Table S10.docx")
)
