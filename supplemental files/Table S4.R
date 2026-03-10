###########
###########
#### GENOME-WIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-25

# Supplemental Table S3

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


# ----------------------------------------------------------------
# Table S3 - GLM results for variance explained by S1 NTD and 
# cases/hospitalizations
# ----------------------------------------------------------------

# state model

# samples per week
samples_state <- dat_sewershed %>%
  group_by(week) %>%
  summarize(seq_samples = sum(samples, na.rm = TRUE),
            conc_samples = sum(conc_samples, na.rm = TRUE)) %>%
  ungroup()

dat_state <- left_join(dat_state, samples_state,
                       by = c("week"))

# remove nas
s <- dat_state %>%
  filter(!is.na(case_incidence)) %>%
  filter(!is.na(ntd_pi_state_3w)) %>%
  filter(!is.na(mean_sars2_conc_state_3w)) %>%
  filter(!is.na(depth_state_3w)) %>%
  filter(!is.na(n_variants_no_thresh_3w_mean))

# gls model
model_time <- gls(case_incidence ~ 
                    + scale(ntd_pi_state_3w)
                  + scale(log(mean_sars2_conc_state_3w))
                  + scale(depth_state_3w)
                  + scale(seq_samples)
                  + scale(conc_samples)
                  + scale(mean_coverage)
                  , data = s,
                  correlation = corAR1(form = ~1|week))

summary(model_time)
performance(model_time)
performance::r2(model_time)

# save output in a table
# table
state_table <- model_summary_function(
  dataframe = s,
  outcome = "case_incidence",
  model = model_time,
  glmer_option = "gls"
)%>%
  mutate(group = "Full model")

#  save table
# edit estimate column
state_table$variable<- c(
  "Intercept",
  "S1 NTD Pi",
  "SARS-CoV-2 concentration",
  "Mean sample depth of read",
  "Number of samples sequenced",
  "Number of samples collected",
  "Mean genome coverage",
  "n",
  "Efron R2",
  "AIC"
)

# round p value
state_table$`p-value` <- round_2(as.numeric(state_table$`p-value`))

# ntd only
# gls model
model_ntd <- gls(case_incidence ~ 
                   + scale(ntd_pi_state_3w)
                 , data = s,
                 correlation = corAR1(form = ~1|week))


# save output in a table
# table
state_table_ntd <- model_summary_function(
  dataframe = s,
  outcome = "case_incidence",
  model = model_ntd,
  glmer_option = "gls"
)%>%
  mutate(group = "S1 NTD model")

#  save table
# edit estimate column
state_table_ntd$variable<- c(
  "Intercept",
  "S1 NTD Pi",
  "n",
  "Efron R2",
  "AIC"
)

# conc only
# gls model
model_conc <- gls(case_incidence ~ 
                    + scale(log(mean_sars2_conc_state_3w))
                  , data = s,
                  correlation = corAR1(form = ~1|week))


# save output in a table
# table
state_table_conc<- model_summary_function(
  dataframe = s,
  outcome = "case_incidence",
  model = model_conc,
  glmer_option = "gls"
)%>%
  mutate(group = "Concentration model")

#  save table
# edit estimate column
state_table_conc$variable<- c(
  "Intercept",
  "SARS-CoV-2 concentration",
  "n",
  "Efron R2",
  "AIC"
)

# make into one table
state_table <- rbind(state_table, state_table_ntd,
                     state_table_conc)
state_table <- state_table %>%
  select(group, everything())

state_table$`p-value` <- round_2(state_table$`p-value`)
# 

# make it an ft table
t <- table_as_flex_function(dataframe = state_table,
                            title = "Tale: Statewide generalized liner model results for S1 NTD")

t <- set_header_labels(t,
                       values = list(
                         group = "Model",
                         variable = "Variable/Metric",
                         Est = "Estimate",
                         `Std. Error` = "Standard Error (SE)",
                         `t value2` = "t value",
                         `Pr(>|t|)` = "P value"
                       ))

t

# save
save_as_docx(FitFlextableToPage(t), 
             path = paste("Supplemental files/",
                          "Table S3.docx",
                          sep = "")
)
