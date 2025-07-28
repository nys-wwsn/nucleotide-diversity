###########
###########
#### GENOMEWIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-07-28

# Main analysis - figures and tables for sequencing diversity paper

# This paper uses several functions and the docstring package to document
# the function usage. You can view help pages for each function by writing in 
# the console ?function_name

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
library(lme4)
library(lmerTest)
library(nlme)

# --------------------------------------
# FUNCTIONS
# --------------------------------------

# make table a flex table
table_as_flex_function <- function(dataframe,
                                   title){
  # docstring information
  #' Change a dataframe to a flextable
  #' 
  #' @description This function changes a dataframe to a flextable object for
  #' creating publication tables. 
  #' 
  #' @param dataframe dataframe. Object to change to flextable.
  #' @param title character string. Title of the table.
  #' 
  #' @return Flextable object.
  
  ft <- flextable(dataframe)
  ft <- theme_vanilla(ft)
  ft <- color(ft, part = "footer", color = "black")
  ft <- set_caption(ft, caption = title)
  ft <- ft %>%
    fontsize(size = 10, part = "all") %>%
    vline(part = "all") %>%
    padding(padding = 0, part = "all") %>%
    align(align = "center") %>%
    align(align = "center", part = "header")
}

# fit to page function
FitFlextableToPage <- function(ft, pgwidth = 6){
  
  # docstring information
  #' Fit flextable output to HTMl word document page.
  #' 
  #' @description This function fits the output from a flextable to the entire
  #' page width of a word document.
  #' 
  #' @param ft flextable object.
  #' @param pgwidth numeric. Page width. Default is 6.
  #' 
  ft_out <- ft %>% autofit()
  
  ft_out <- 
    width(ft_out, width = 
            dim(ft_out)$widths*pgwidth /(flextable_dim(ft_out)$widths))
  return(ft_out)
}

# lag / lead time function for case incidence
lag_lead_function <- function(dataframe,
                              columns,
                              lags){
  
  # docstring information
  #' Estimate the optimal lead or or lag time for time series correlated data 
  #' (case data).
  #' 
  #' @description This function calculates Spearman and Pearson correlations for
  #' the diversity data and case incidence.
  #' 
  #' @param dataframe dataframe. Dataframe containing case incidence and 
  #' nucleotid diveristy measures.
  #' @param columns character vector. Vector of column names in dataframe to 
  #' pass to the correlation calculation step. For example, the column name for 
  #' S1 NTD Pi.
  #' @param lags integer vector. Number of lagged dates to test. For example, to 
  #' test  up to 10 weeks lag, the parameter would be seq(1:10)
  #' 
  #' @details Nucleotide diveristy presents a leading indicator in the time
  #' series plots. This function calculates correlations for leading or lagging
  #' diversity measures with case incidence. As lag increases, a new correlation
  #' is estimated for each lag.

  #' @return Dataframe with Pearson and Spearman correlations for each column
  #' passed in the `columns` parameter.
  
  #' @examples
  #' pi_lag <- lag_lead_function_hosp(dataframe = dat_state,
  #' columns = c("genomewide_pi_state_3w",
  #'             "spike_pi_state_3w",
  #'             "orf_pi_state_3w",
  #'             "ntd_pi_state_3w",
  #'             "rbd_pi_state_3w",
  #'             "cov_mt_2_pi_state_3w"),
  #' lags = seq(1:10))
  #' 
  # lag 0 correlation
  cor_values_0 <- dataframe %>%
    pivot_longer(cols = columns
    ) %>%
    group_by(group, name) %>%
    summarize(spearman_cor = cor(method = "spearman",
                                 value, case_incidence,
                                 use = "na.or.complete"),
              pearson_cor = cor(method = "pearson",
                                value, case_incidence,
                                use = "na.or.complete")
    )
  lag_list <- list()
  # lag 
  for(i in lags){
    cor_values <- dataframe %>%
      pivot_longer(cols = columns
      ) %>%
      group_by(group, name) %>%
      arrange(week) %>%
      mutate(value = lag(value, i, default = NA)) %>%
      ungroup() %>%
      group_by(group, name) %>%
      summarize(spearman_cor = cor(method = "spearman",
                                   value, case_incidence,
                                   use = "na.or.complete"),
                pearson_cor = cor(method = "pearson",
                                  value, case_incidence,
                                  use = "na.or.complete")
      )
    
    # add lag as field in data
    cor_values$lag <- i
    
    # add into list object
    lag_list[[i]] <- cor_values
    
  }
  
  # make into df
  cor_pi_lag <- do.call(rbind, lag_list) 
  
  # try lead time
  leads <- lags
  
  lead_list <- list()
  
  for(i in leads){
    cor_values <- dataframe%>%
      pivot_longer(cols = columns
      ) %>%
      group_by(group, name) %>%
      arrange(week) %>%
      mutate(value = lead(value, i, default = NA)) %>%
      ungroup() %>%
      group_by(group, name) %>%
      summarize(spearman_cor = cor(method = "spearman",
                                   value, case_incidence,
                                   use = "na.or.complete"),
                pearson_cor = cor(method = "pearson",
                                  value, case_incidence,
                                  use = "na.or.complete")
      )
    
    # add lag as field in data
    cor_values$lag <- i * -1
    
    # add into list object
    lead_list[[i]] <- cor_values
    
  }
  
  # make into df
  cor_pi_lead <- do.call(rbind, lead_list) %>%
    bind_rows(cor_values_0) %>%
    mutate(lag = ifelse(is.na(lag), 0, lag)
    ) %>%
    bind_rows(cor_pi_lag)
  
  return(cor_pi_lead)
}

# lag / lead time function for hospital admissions
lag_lead_function_hosp <- function(dataframe,
                                   columns,
                                   lags){
  
  # docstring information
  #' Estimate the optimal lead or or lag time for time series correlated data 
  #' (hospital admissions).
  #' 
  #' @description This function calculates Spearman and Pearson correlations for
  #' the diversity data and hospital incidence.
  #' 
  #' @param dataframe dataframe. Dataframe containing hospital incidence and 
  #' nucleotide diveristy measures.
  #' @param columns character vector. Vector of column names in dataframe to 
  #' pass to the correlation calculation step. For example, the column name for 
  #' S1 NTD Pi.
  #' @param lags integer vector. Number of lagged dates to test. For example, to 
  #' test up to 10 weeks lag, the parameter would be seq(1:10)
  #' 
  #' @details Nucleotide diveristy presents a leading indicator in the time
  #' series plots. This function calculates correlations for leading or lagging
  #' diversity measures with hospital incidence. As lag increases, a new correlation
  #' is estimated for each lag.
  
  #' @return Dataframe with Pearson and Spearman correlations for each column
  #' passed in the `columns` parameter.
  
  #' @examples
  #' pi_lag <- lag_lead_function_hosp(dataframe = dat_state,
  #' columns = c("genomewide_pi_state_3w",
  #'             "spike_pi_state_3w",
  #'             "orf_pi_state_3w",
  #'             "ntd_pi_state_3w",
  #'             "rbd_pi_state_3w",
  #'             "cov_mt_2_pi_state_3w"),
  #' lags = seq(1:10))
  # lag 0 correlation
  cor_values_0 <- dataframe %>%
    pivot_longer(cols = columns
    ) %>%
    group_by(group, name) %>%
    summarize(spearman_cor = cor(method = "spearman",
                                 value, hosp_incidence,
                                 use = "na.or.complete"),
              pearson_cor = cor(method = "pearson",
                                value, hosp_incidence,
                                use = "na.or.complete")
    )
  lag_list <- list()
  # lag 
  for(i in lags){
    cor_values <- dataframe %>%
      pivot_longer(cols = columns
      ) %>%
      group_by(group, name) %>%
      arrange(week) %>%
      mutate(value = lag(value, i, default = NA)) %>%
      ungroup() %>%
      group_by(group, name) %>%
      summarize(spearman_cor = cor(method = "spearman",
                                   value, hosp_incidence,
                                   use = "na.or.complete"),
                pearson_cor = cor(method = "pearson",
                                  value, hosp_incidence,
                                  use = "na.or.complete")
      )
    
    # add lag as field in data
    cor_values$lag <- i
    
    # add into list object
    lag_list[[i]] <- cor_values
    
  }
  
  # make into df
  cor_pi_lag <- do.call(rbind, lag_list) 
  
  # try lead time
  leads <- lags
  
  lead_list <- list()
  
  for(i in leads){
    cor_values <- dataframe%>%
      pivot_longer(cols = columns
      ) %>%
      group_by(group, name) %>%
      arrange(week) %>%
      mutate(value = lead(value, i, default = NA)) %>%
      ungroup() %>%
      group_by(group, name) %>%
      summarize(spearman_cor = cor(method = "spearman",
                                   value, hosp_incidence,
                                   use = "na.or.complete"),
                pearson_cor = cor(method = "pearson",
                                  value, hosp_incidence,
                                  use = "na.or.complete")
      )
    
    # add lag as field in data
    cor_values$lag <- i * -1
    
    # add into list object
    lead_list[[i]] <- cor_values
    
  }
  
  # make into df
  cor_pi_lead <- do.call(rbind, lead_list) %>%
    bind_rows(cor_values_0) %>%
    mutate(lag = ifelse(is.na(lag), 0, lag)
    ) %>%
    bind_rows(cor_pi_lag)
  
  return(cor_pi_lead)
}

# round to two digits
round_2 <- function(x){
  # docstring information
  #' Round numeric data to 3 places
  #' 
  #' @description This function will round values to three places.
  #' 
  #' @param x vector. Numeric vector.
  round(x, 2)
}

# two signifi digits
signif_2 <- function(x){
  # docstring information
  #' Report 3 significant digits
  #' 
  #' @description This function will report 3 significant figures for a numeric
  #' vector
  #' 
  #' @param x vector. Numeric vector.
  signif(x, 2)
}

# round_signifi_function <- function(x){
#   if(x >1){
#     round(x, 3)
#   } else if(x <= 1){
#     signif(x, 3)
#   }
# }

# for results, round values >1, if <1, report 3 significant figures
round_signifi_function <- function(x){
  # docstring information
  #' Round numeric data
  #' 
  #' @description This function will round values greater than 1 to 3 places and
  #' report 3 significant digits for values less than 1.
  #' 
  #' @param x vector. Numeric vector.

  ifelse(x >1,
    round(x, 3),
    signif(x, 3)
    )
}

# model summary function
model_summary_function <- function(dataframe,
                                   outcome,
                                   model,
                                   glmer_option){
  # docstring information
  #' Output model summary information to a table
  #' 
  #' @description This function will take a model object in R and output a
  #' dataframe with the coefficients and model metrics.
  #' 
  #' @param dataframe dataframe. Dataframe containing the data used to fit the 
  #' model.
  #' @param outcome character vector. Outcome variable name from the model.
  #' @param model model object. Model fit to the data.
  #' @param glemer_option character string. Whether the model was fit using the
  #' glmer function.
  #' 
  #' @details This function will take a model fit using the data and output
  #' the coefficients, intercept, standard errors, t values, and p values
  #' for the model. Further, essential model metrics including the number of 
  #' observation, AIC, and Efron R2 are returned.
  
  #' @return Dataframe 

  if(glmer_option == "No"){
    
    dataframe$outcome_name <- dataframe[[outcome]]
    # Effron R2 - compares the actual to the predicted values 
    # against the residuals
    Actual <- dataframe$case_incidence
    Predicted = predict(model, type="response")
    Residuals = residuals(model)
    
    effron_r <- efronRSquared(residual = Residuals, 
                              predicted = Predicted, 
                              statistic = "EfronRSquared")
    effron_r <- as.data.frame(cbind("Efron R2", as.numeric(effron_r)))
    colnames(effron_r) <- c("variable", "Est")
    aic <- AIC(model)
    aic <- as.data.frame(cbind("AIC", aic))
    colnames(aic) <- c("variable", "Est")
    
    # model coefficients, se, and p values
    cc <- coef(summary(model))
    cc <- within(as.data.frame(cc),
                 {   `Std. Error` <- `Std. Error`
                 `t value` <- Estimate/`Std. Error`
                 `Pr(>|t|)` <- 2*pt(-abs(`t value`), df=length(dataframe)-1)
                 })
    
    #add ci
    #calculate se from dataframe
    se <- (cc$`Std. Error`)
    # table of estimates with 95% CI
    (tab <- cbind(Est = cc$Estimate, LL = cc$Estimate - 1.96 * se, 
                  UL = cc$Estimate + 1.96 *
                    se))
    tab2 <- cbind(cc, tab)
    
    # add stars based on the p value 
    tab2$`t value2` <- ifelse(tab2$`Pr(>|t|)` < 0.1, 
                              paste(round(tab2$`t value`, 2), "*", sep = ""), 
                              round(tab2$`t value`, 2))
    tab2$`t value2` <- ifelse(tab2$`Pr(>|t|)` < 0.05, 
                              paste(round(tab2$`t value`, 2), "**", sep = ""), 
                              round(tab2$`t value`, 2))
    tab2$`t value2` <- ifelse(tab2$`Pr(>|t|)` < 0.01, 
                              paste(round(tab2$`t value`, 2), "***", sep = ""), 
                              round(tab2$`t value`, 2))
    # add variable names
    tab2$variable <- row.names(tab2)
    # grab columns and export table
    tab_export <- tab2 %>% dplyr::select(variable, Est, `Std. Error`, `t value2`, 
                                         `Pr(>|t|)`)
    # add number of observations
    n <- nobs(model)
    n <- as.data.frame(cbind("n", n))
    colnames(n) <- c("variable", "Est")
    
    tab_export$Est <- as.character(tab_export$Est)
    tab_export_final <- bind_rows(tab_export, n, effron_r, aic)
    
    tab_export_final$Est <- round_signifi_function(
      as.numeric(tab_export_final$Est)
    )
    tab_export_final$`Std. Error` <- round_signifi_function(
      as.numeric(tab_export_final$`Std. Error`)
    )
    tab_export_final$`Pr(>|t|)` <- round_signifi_function(
      as.numeric(tab_export_final$`Pr(>|t|)`)
    )
    
    # round values
    tab_export_final <- tab_export_final %>%
      mutate_if(., is.numeric, round_signifi_function)
    
    return(tab_export_final)
  } else if(glmer_option == "Yes"){
    
    dataframe$outcome_name <- dataframe[[outcome]]
    
    # Effron R2 - compares the actual to the predicted values 
    # against the residuals
    Actual <- dataframe$case_incidence
    Predicted = predict(model, type="response")
    Residuals = residuals(model)
    
    effron_r <- efronRSquared(residual = Residuals, 
                              predicted = Predicted, 
                              statistic = "EfronRSquared")
    effron_r <- as.data.frame(cbind("Efron R2", as.numeric(effron_r)))
    colnames(effron_r) <- c("variable", "Est")
    aic <- AIC(model)
    aic <- as.data.frame(cbind("AIC", aic))
    colnames(aic) <- c("variable", "Est")
    
    # model coefficients, se, and p values
    cc <- coef(summary(model))
    cc <- within(as.data.frame(cc),
                 {   `Std. Error` <- `Std. Error`
                 `t value` <- Estimate/`Std. Error`
                 `Pr(>|t|)` <- 2*pt(-abs(`t value`), df=length(dataframe)-1)
                 })
    
    #add ci
    #calculate se from dataframe
    se <- (cc$`Std. Error`)
    # table of estimates with 95% CI
    (tab <- cbind(Est = cc$Estimate, LL = cc$Estimate - 1.96 * se, UL = cc$Estimate + 1.96 *
                    se))
    tab2 <- cbind(cc, tab)
    
    # add stars based on the p value 
    tab2$`t value2` <- ifelse(tab2$`Pr(>|t|)` < 0.1, 
                              paste(round(tab2$`t value`, 2), "*", sep = ""), 
                              round(tab2$`t value`, 2))
    tab2$`t value2` <- ifelse(tab2$`Pr(>|t|)` < 0.05, 
                              paste(round(tab2$`t value`, 2), "**", sep = ""), 
                              round(tab2$`t value`, 2))
    tab2$`t value2` <- ifelse(tab2$`Pr(>|t|)` < 0.01, 
                              paste(round(tab2$`t value`, 2), "***", sep = ""), 
                              round(tab2$`t value`, 2))
    
    # add variable names
    tab2$variable <- row.names(tab2)
    # grab columns and export table
    tab_export <- tab2 %>% dplyr::select(variable, Est, `Std. Error`, `t value2`, 
                                         `Pr(>|t|)`)
    # add number of observations
    n <- nobs(model)
    n <- as.data.frame(cbind("n", n))
    colnames(n) <- c("variable", "Est")
    
    tab_export$Est <- as.character(tab_export$Est)
    tab_export_final <- bind_rows(tab_export, n, effron_r, aic)
    
    tab_export_final$Est <- as.numeric(tab_export_final$Est)
    
    # round values
    tab_export_final <- tab_export_final %>%
      mutate_if(., is.numeric, round_signifi_function)
    
    return(tab_export_final)
  } else if(glmer_option == "gls"){
    
    dataframe$outcome_name <- dataframe[[outcome]]
    
    # Effron R2 - compares the actual to the predicted values 
    # against the residuals
    Actual <- dataframe$case_incidence
    Predicted = predict(model, type="response")
    Residuals = residuals(model)
    
    effron_r <- efronRSquared(residual = Residuals, 
                              predicted = Predicted, 
                              statistic = "EfronRSquared")
    effron_r <- as.data.frame(cbind("Efron R2", as.numeric(effron_r)))
    colnames(effron_r) <- c("variable", "Est")
    aic <- AIC(model)
    aic <- as.data.frame(cbind("AIC", aic))
    colnames(aic) <- c("variable", "Est")
    
    # model coefficients, se, and p values
    cs <- as.data.frame(summary(model)$tTable)
    cc <- within(as.data.frame(cs),
                 {   `Std.Error` <- `Std.Error`
                 `t-value` <- Value/`Std.Error`
                 `p-value` <- 2*pt(-abs(`t-value`), df=length(dataframe)-1)
                 })
    
    #add ci
    #calculate se from dataframe
    se <- (cc$`Std.Error`)
    # table of estimates with 95% CI
    (tab <- cbind(Est = cc$Value, LL = cc$Value - 1.96 * se, UL = cc$Value + 1.96 *
                    se))
    tab2 <- cbind(cc, tab)
    
    # add stars based on the p value 
    tab2$`t value2` <- ifelse(tab2$`p-value` < 0.1, 
                              paste(round(tab2$`t-value`, 2), "*", sep = ""), 
                              round(tab2$`t-value`, 2))
    tab2$`t value2` <- ifelse(tab2$`p-value` < 0.05, 
                              paste(round(tab2$`t-value`, 2), "**", sep = ""), 
                              round(tab2$`t-value`, 2))
    tab2$`t value2` <- ifelse(tab2$`p-value` < 0.01, 
                              paste(round(tab2$`t-value`, 2), "***", sep = ""), 
                              round(tab2$`t-value`, 2))
    
    # add variable names
    tab2$variable <- row.names(tab2)
    # grab columns and export table
    tab_export <- tab2 %>% dplyr::select(variable, Est, `Std.Error`, `t value2`, 
                                         `p-value`)
    # add number of observations
    n <- nobs(model)
    n <- as.data.frame(cbind("n", n))
    colnames(n) <- c("variable", "Est")
    
    tab_export$Est <- as.character(tab_export$Est)
    tab_export_final <- bind_rows(tab_export, n, effron_r, aic)
    
    tab_export_final$Est <- as.numeric(tab_export_final$Est)
    
    # round values
    tab_export_final <- tab_export_final %>%
      mutate_if(., is.numeric, round_signifi_function)
    
    return(tab_export_final)
  } 
  
}

# extract legend function
g_legend<-function(a.gplot){
  # docstring information
  #' Extract legend from ggplot object
  #' 
  #' @description This function will pull the legend from a ggplot as a grob
  #' and let the user add the legend to a panel figure.
  #' 
  #' @param a.gplot ggplot object.
  #' 
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# -----------------------------------------------------------------------------

# --------------------------------------
# PLOT THEMES
# --------------------------------------

font_add_google("Roboto Condensed", family = "roboto_c")

# map
theme_dth_maps <- 
  theme_void(base_size = 12)+
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(family = "roboto_c"),
        axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray80"),
        panel.grid.minor.y = element_line(color = "gray80"),
        strip.background = element_blank(),
        legend.title = element_text(hjust = 0.5),
        legend.background = element_rect(color = "gray80"),
        plot.subtitle = element_text(hjust=0.5)
  )

# plot theme
theme_dth_1 <- 
  theme_bw(base_size = 12)+
  theme(plot.title = element_text(hjust = 0.5),
        text=element_text(family = "roboto_c"),
        axis.ticks = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "gray80"),
        panel.grid.minor.y = element_line(color = "gray80"),
        strip.background = element_blank(),
        legend.title = element_text(hjust = 0.5),
        legend.background = element_rect(color = "gray80")
  )

pal <- met.brewer(name = "Austria",
                  n = 6)



# plot for sample data with sequencing locations
plot_bp_labels <-ggplot()+
  annotate("rect",
           xmin = 0, xmax = 30000,
           ymin = -0.0025, ymax = 0,
           alpha = 0.5,
           fill = "green")+
  annotate("text",
           label = "SARS-CoV-2 Genome",
           x = 3000,
           y = -0.0007
  )+
  annotate("rect",
           xmin = 10063, xmax = 11842,
           ymin = -0.0025, ymax = 0,
           fill = pal[4])+
  annotate("text",
           label = "NSPs\n5+6",
           x = 11000,
           y = -0.0008,
           size = 3.5
  )+
  annotate("rect",
           xmin = 20661, xmax = 21549,
           ymin = -0.0025, ymax = 0,
           fill = pal[1])+
  annotate("text",
           label = "CoV\nMT2",
           x = 21150,
           y = -0.0008,
           color = "white",
           size = 3
  )+
  annotate("rect",
           xmin = 21598, xmax = 22474,
           ymin = -0.0025, ymax = 0,
           fill = pal[3])+
  annotate("text",
           label = "NTD",
           x = 22050,
           y = -0.0008,
           color = "white",
           size = 3.5,
           angle = 90
  )+
  annotate("rect",
           xmin = 22516, xmax = 23185,
           ymin = -0.0015, ymax = 0,
           fill = pal[2])+
  annotate("text",
           label = "RBD",
           x = 22900,
           y = -0.0006,
           color = "white",
           size = 3,
           angle = 90
  )+
  annotate("rect",
           xmin = 21563, xmax = 25384,
           ymin = -0.0025, ymax = -0.0015,
           fill = pal[5])+
  annotate("text",
           label = "Spike protein",
           x = 23250,
           y = -0.002,
           color = "black",
           size = 4
  )

# -----------------------------------------------------------------------------

# --------------------------------------
# DATA LOAD
# --------------------------------------

# combined data files for each geography
load(file = "data/combined_data.Rdata")

# max and min weeks for the study
min_week <- "2023-01-01"
max_week <- "2025-04-20"

# load county spatial data
ny_counties <- st_read("data/Counties-Shoreline/Counties_Shoreline.shp") %>%
  mutate(county = NAME)
ny_counties$county[ny_counties$county == "St Lawrence"] <- "St Lawrence"

# -----------------------------------------------------------------------------

# --------------------------------------
# ANALYSIS
# --------------------------------------

ggplot(data = dat_state,
       aes(x = ntd_pi_state_3w,
           y = case_incidence))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  facet_wrap(~group)

ggplot(data = dat_state,
       aes(x = ntd_h_state_3w,
           y = case_incidence))+
  geom_point()+
  ggpubr::stat_cor(method = "spearman")+
  facet_wrap(~group)


cor(dat_state$n_variants_no_thresh_3w_mean,
    dat_state$case_incidence,
    method = "spearman",
    use = "na.or.complete")

# -------------------------------------------------------------------------

# --------------------------------------
# DESCRIPTIVE STATISTICS
# --------------------------------------

# Table for diversity data

# merge lab data
# summarize by pcr lab
diversity_df_summary <- dat_sewershed %>%
  group_by(pcr_lab) %>%
  filter(!is.na(genomewide_pi)) %>%
  summarize(`Number of sampling sites` = length(unique(facility_id)),
            `Number of samples` = sum(samples, na.rm = TRUE),
            `First collection week` = min(week, na.rm = TRUE),
            `Last collection week` = max(week, na.rm= TRUE),
            `Mean Pi` = signif(mean(genomewide_pi, na.rm = TRUE),2),
            `Min Pi` = signif(min(genomewide_pi, na.rm = TRUE), 2),
            `Median Pi` = signif(median(genomewide_pi, na.rm = TRUE), 2),
            `Max Pi` = signif(max(genomewide_pi, na.rm = TRUE), 2),
            `SD Pi` = signif(sd(genomewide_pi, na.rm = TRUE), 2)
            
  ) %>%
  ungroup() %>%
  mutate(pcr_lab = case_when(
    pcr_lab == "NYC" ~ "NYC DOH",
    pcr_lab == "genesee_orleans_health" ~ "GO Health",
    pcr_lab == "quadrant" ~ "Quadrant",
    pcr_lab == "suny_buffalo" ~ "SUNY Buffalo",
    pcr_lab == "suny_stony_brook" ~ "SUNY Stony Brook",
    pcr_lab == "wadsworth" ~ "Wadsworth Center"
  )) %>%
  rename(`PCR Lab` = pcr_lab) 


# save table as flextable with rounded values
pi_ft <- table_as_flex_function(
  dataframe = diversity_df_summary,
  title = "Table: Descriptive statstics for genomewide Pi")

# save
save_as_docx(FitFlextableToPage(pi_ft), 
             path = paste("figures and tables/",
                          "Table - genomewide pi descriptive statistics.docx")
)

# --------------------------------------------------------------------------

# Table for conc data
# summarize by pcr lab
conc_df_summary <- dat_sewershed %>%
  ungroup() %>%
  filter(!is.na(pcr_lab)) %>%
  group_by(pcr_lab) %>%
  summarize(`Number of sampling sites` = length(unique(facility_id)),
            `Number of samples` = sum(samples, na.rm  = TRUE),
            `First collection week` = min(week, na.rm = TRUE),
            `Last collection week` = max(week, na.rm= TRUE),
            `Mean concentration` = signif(mean(mean_sars2_conc, na.rm = TRUE),
                                          2),
            `Min concentration` = signif(min(mean_sars2_conc, na.rm = TRUE),
                                         2),
            `Median concentration` = signif(median(mean_sars2_conc, 
                                                   na.rm = TRUE),
                                            2),
            `Max concentration` = signif(max(mean_sars2_conc, na.rm = TRUE), 2),
            `SD concentration` = signif(sd(mean_sars2_conc, na.rm = TRUE),
                                        2)
            
  ) %>%
  ungroup() %>%
  mutate(pcr_lab = case_when(
    pcr_lab == "NYC" ~ "NYC DOH",
    pcr_lab == "genesee_orleans_health" ~ "GO Health",
    pcr_lab == "quadrant" ~ "Quadrant",
    pcr_lab == "suny_buffalo" ~ "SUNY Buffalo",
    pcr_lab == "suny_stony_brook" ~ "SUNY Stony Brook",
    pcr_lab == "wadsworth" ~ "Wadsworth Center"
  )) %>%
  rename(`PCR Lab` = pcr_lab) %>%
  mutate(`Min concentration` = ifelse(
    `Min concentration` == 1, 0, `Min concentration`
  ))

# make it flextable
conc_ft <- table_as_flex_function(
  dataframe = conc_df_summary,
  title = "Table: Descriptive statistics for concentration data")

# save
save_as_docx(FitFlextableToPage(conc_ft), 
             path = paste0("figures and tables",
                           "/Table - concentration descriptive statistics.docx")
)

# -------------------------------------------------------------------

# Table for case + hosp data

hosp_summary <- dat_county %>%
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
  select(group, everything())

# make it an ft table
# make it flextable
clinical_ft <- table_as_flex_function(
  dataframe = clinical_df_summary,
  title = "Table: Descriptive statistics for COVID-19 clinical data")

# save
save_as_docx(FitFlextableToPage(clinical_ft), 
             path = paste0("figures and tables/",
                           "Table - clinical descriptive statistics.docx")
)

# --------------------------------------------------------------------

# MAP OF THE SAMPLING SITES AND POSITIVE DETECTION RATE
# summarize by pcr lab
detection_rate <- dat_sewershed %>%
  filter(week >= ymd(min_week) & week <= ymd(max_week)) %>%
  filter(!is.na(genomewide_pi)) %>%
  mutate(sars_positive = ifelse(mean_sars2_conc > 1, 1, 0)) %>%
  group_by(facility_id, wwtp_latitude, wwtp_longitude,capacity_mgd
           ) %>%
  summarize(`Number of positive samples` = sum(sars_positive, na.rm = TRUE),
            `Number of samples` = n()
  ) %>%
  ungroup() %>%
  mutate(`Detection rate` = `Number of positive samples`/`Number of samples`) 

# plot by site

# make the facility data into points
detection_rate <- detection_rate %>%
  st_as_sf(coords =  c("wwtp_longitude", "wwtp_latitude"),
           crs = "+proj=longlat +zone=18 +datum=WGS84") %>% 
  st_transform(st_crs(ny_counties))


# map
map_dr <- 
  ggplot()+
  geom_sf(data = ny_counties,
          #aes(fill = pcr_lab),
          fill = "white",
          color = "grey65"
          
  )+
  geom_sf(data = detection_rate,
          aes(color = as.numeric(`Detection rate`),
              size = capacity_mgd),
          alpha = 0.7)+
  theme_dth_maps +
  scale_color_gradientn(colors = rev(MetBrewer::met.brewer("Derain")),
                        name = "Detection\nrate")+
  theme(legend.position.inside = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.direction = "horizontal")+
  labs(size = "WWTP size\n(mgd)",
       title = "SARS-CoV-2 detection rate in wastewater")+
  theme(legend.position = c(0.5, 0.3),
        legend.background = element_blank(),
        legend.box.just = "left")
# guides(fill = guide_legend(nrow=2, byrow=TRUE))+
# guides(color = guide_legend(nrow=2, byrow=TRUE))

map_dr

# sample collected over time and samples sequenced stacked bar
samples_time <- dat_sewershed %>%
  filter(week >= ymd(min_week) & week <= ymd(max_week)) %>%
  filter(!is.na(mean_sars2_conc)) %>%
  group_by(week) %>%
  summarise(samples_collected = sum(conc_samples, na.rm = TRUE)
         ) %>%
  ungroup() 

# sequenced samples over time
seq_time <- dat_sewershed  %>%
  filter(week >= ymd(min_week) & week <= ymd(max_week)) %>%
  filter(!is.na(genomewide_pi)) %>%
  group_by(week) %>%
  summarise(samples_sequenced = sum(samples, na.rm = TRUE)
  ) %>%
  ungroup() 
  
# join
samples_time <- left_join(samples_time, seq_time, by = c("week"))
samples_time$samples_no_seq <- 
  samples_time$samples_collected - samples_time$samples_sequenced

samples_time <- samples_time %>%
  pivot_longer(cols = c(samples_sequenced, samples_no_seq))

# combine into one panel figure
samples_plot <- 
  ggplot(data = samples_time) +
  geom_bar(position = "stack",
           stat = "identity",
           aes(x = week, y = value,
               fill = name))+
  theme_dth_1+
  theme(legend.position = "bottom")+
  labs(x = "",
       y = "Samples collected",
       fill = "",
       title = "Samples collected over time")+
  scale_fill_manual(values = MetBrewer::met.brewer(n = 2,
                                                   name = "Egypt"),
                    labels = c("Samples not sent for sequencing",
                               "Samples sequenced"))

samples_plot 

# panel
ggarrange(map_dr, samples_plot , nrow = 1, labels = c("A","B"))

# save map
png("figures and tables/Figure-Detection rate map and samples over time.png",
    units = "in",
    width = 12, height = 6,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
ggarrange(map_dr, samples_plot , nrow = 1, labels = c("A","B"))
showtext_end()
dev.off()


# ------------------------------------------------
# SINGLE SAMPLE FIGURES
# ------------------------------------------------

# bird island
bird_pi <- readRDS("data/erie_NY029028410A.rds")%>%
  mutate(week = floor_date(ymd(date), unit = "weeks", week_start = 7)
  )

# meadowbrook and limestone
meadow_pi <- readRDS("data/onondaga_NY067027723A.rds")%>%
  mutate(week = floor_date(ymd(date), unit = "weeks", week_start = 7)
  )

# plot high transmission time period
# 2024-07-28
# select date and stack samples
# note using date for nyc since they sample twice weekly, both are analyzed
bird_high <- bird_pi %>%
  filter(week == "2023-12-10") %>%
  mutate(group = "Erie County - Bird Island WWTP")
meadow_high <- meadow_pi %>%
  filter(week == "2023-12-10") %>%
  mutate(group = "Onondaga County - Meadowbrook/Limestone WWTP")

# combine into one df
high_sample_df <- bind_rows(bird_high, meadow_high) %>%
  mutate(group = factor(group,
                        levels = c("Erie County - Bird Island WWTP",
                                   "Onondaga County - Meadowbrook/Limestone WWTP")))

shannon_high <- plot_bp_labels+
  geom_line(data = high_sample_df,
            aes(x = center,
                y = avg_h,
                color = group),
            lwd = 1)+
  theme_dth_1+
  theme(legend.position = "bottom",
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "Shannon H for select samples\nfrom the week of December 10, 2023",
       subtitle = "Average daily COVID-19 incidence > 75 cases per 100,000",
       y = "Shannon's H",
       x = "Base pair position")+
  scale_color_manual(values = met.brewer(name = "Hiroshige",
                                         n = 3),
                     labels = c("Large WWTP\nApprox pop served 431,000",
                                "Small WWTP\nApprox pop served 39,000"),
                     name = "")+
  scale_y_continuous(
    minor_breaks = NULL,
    breaks = seq(-0 , 0.05, 0.005) ,
    limits = c(-0.0025,0.03))

shannon_high

#pi high
pi_high <- plot_bp_labels+
  geom_line(data = high_sample_df,
            aes(x = center,
                y = avg_pi,
                color = group),
            lwd = 1)+
  theme_dth_1+
  theme(legend.position = "bottom",
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = expression("Pi"[w]* " for select samples from the week of December 10, 2023"),
       subtitle = "Average daily COVID-19 incidence > 75 cases per 100,000",
       y = expression("Pi"[w]),
       x = "Base pair position")+
  scale_color_manual(values = met.brewer(name = "Hiroshige",
                                         n = 3),
                     labels = c("Large WWTP\nApprox pop served 431,000",
                                "Small WWTP\nApprox pop served 39,000"),
                     name = "")+
  scale_y_continuous(
    minor_breaks = NULL,
    breaks = seq(-0 , 0.05, 0.005) ,
    limits = c(-0.0025,0.017))
pi_high

# low transmission time
# 2024-05-05

bird_low <- bird_pi %>%
  filter(ymd(date) ==  "2023-06-06") %>%
  mutate(group = "Erie County - Bird Island WWTP")
meadow_low <- meadow_pi %>%
  filter(week == "2023-06-11") %>%
  mutate(group = "Onondaga County - Meadowbrook/Limestone WWTP")

# combine into one df
low_sample_df <- bind_rows(bird_low, meadow_low)%>%
  mutate(group = factor(group,
                        levels = c("Erie County - Bird Island WWTP",
                                   "Onondaga County - Meadowbrook/Limestone WWTP")))

# shannon
shannon_low <- plot_bp_labels+
  geom_line(data = low_sample_df,
            aes(x = center,
                y = avg_h,
                color = group),
            lwd = 1)+
  theme_dth_1+
  theme(legend.position = "bottom",
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = "Shannon H for select samples\nfrom the week of June 11, 2023",
       subtitle = "Average daily COVID-19 incidence  < 10 cases per 100,000",
       y = "Shannon's H",
       x = "Base pair position")+
  scale_color_manual(values = met.brewer(name = "Hiroshige",
                                         n = 3),
                     labels = c("Large WWTP\nApprox pop served 431,000",
                                "Small WWTP\nApprox pop served 39,000"),
                     name = "")+
  scale_y_continuous(
    minor_breaks = NULL,
    breaks = seq(-0 , 0.05, 0.005) ,
    limits = c(-0.0025,0.03))

shannon_low

pi_low <- plot_bp_labels+
  geom_line(data = low_sample_df,
            aes(x = center,
                y = avg_pi,
                color = group),
            lwd = 1)+
  theme_dth_1+
  theme(legend.position = "bottom",
        plot.subtitle = element_text(hjust = 0.5))+
  labs(title = expression("Pi"[w]* " for select samples from the of June 11, 2023"),
                          subtitle = "Average daily COVID-19 incidence  < 10 cases per 100,000",
       y = expression("Pi"[w]),
       x = "Base pair position")+
  scale_color_manual(values = met.brewer(name = "Hiroshige",
                                         n = 3),
                     labels = c("Large WWTP\nApprox pop served 431,000",
                                "Small WWTP\nApprox pop served 39,000"),
                     name = "")+
  scale_y_continuous(
    minor_breaks = NULL,
    breaks = seq(-0 , 0.05, 0.005) ,
    limits = c(-0.0025,0.017))
pi_low

mylegend <- g_legend(pi_high)

ggarrange(ggarrange(pi_low+ theme(legend.position = "none"),
                    pi_high + theme(legend.position = "none"),
                    nrow = 1,
                    labels = c("A", "B")),
          ggarrange(shannon_low+ theme(legend.position = "none"),
                    shannon_high + theme(legend.position = "none"),
                    nrow = 1,
                    labels = c("C", "D")),
          mylegend, 
          nrow=3,
          heights=c(5,5, 1))

# save plots
png("figures and tables/Figure-sample pi and shannon values.png",
    units = "in",
    width = 9.5, height = 11,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
ggarrange(ggarrange(pi_low+ theme(legend.position = "none"),
                    pi_high + theme(legend.position = "none"),
                    nrow = 1,
                    labels = c("A", "B")),
          ggarrange(shannon_low+ theme(legend.position = "none"),
                    shannon_high + theme(legend.position = "none"),
                    nrow = 1,
                    labels = c("C", "D")),
          mylegend, 
          nrow=3,
          heights=c(5,5, 1))
dev.off()


# --------------------------------------
# TIME SERIES AND CORRELATIONS COMBINED
# --------------------------------------

# Statewide time series figure, all genome regions, panel for NYC
dat_state_long <- dat_state %>%
  pivot_longer(cols = c(
    genomewide_pi_state_3w,
    ntd_pi_state_3w,
    cov_mt_2_pi_state_3w,
    orf_pi_state_3w,
    spike_pi_state_3w,
    rbd_pi_state_3w
  ))%>%
  mutate(name_factor = factor(name,
                              levels = c("genomewide_pi_state_3w",
                                         "spike_pi_state_3w",
                                         "orf_pi_state_3w",
                                         "ntd_pi_state_3w",
                                         "rbd_pi_state_3w",
                                         "cov_mt_2_pi_state_3w"),
                              labels = c("Genome",
                                         "Spike",
                                         "ORF NSPs 5 and 6",
                                         "S1 NTD",
                                         "S1 RBD",
                                         "CoV-MT-2"))
  )

pal <- met.brewer(name = "Troy", n = 8)
pal

time_series_1 <- 
  ggplot(data = dat_state_long) +
  geom_bar(aes(x = week, y = case_incidence,
               fill = "grey60"),
           position = "dodge", stat = "identity")+
  geom_line(aes(x = week, y = value * 10000,
                color = name),
            lwd= 1.15)+
  facet_wrap(~name_factor, nrow = 6)+
  theme_dth_1+
  theme(legend.position = "none",
        legend.background = element_rect(color = "white"),
        axis.text.x = element_text(angle = 90, vjust = 0)
  )+
  scale_color_manual(values = MetBrewer::met.brewer(n = 6,
                                                    name = "Austria"),
                     name = "",
                     labels = c("CoV-MT-2",
                                "Genomewide",
                                "S1 NTD",
                                "ORF NSPs 5 and 6",
                                "S1 RBD",
                                "Spike"))+
  scale_x_date(date_labels = "%b %y",
               date_breaks = "1 month",
               limits = c(ymd("2023-01-01"),
                          ymd("2025-04-20"))
  )+
  scale_fill_manual(values = c("grey60" = "grey60"),
                    labels = c("Cases\nper 100k"),
                    name = ""
  )+
  labs(x = "",
       title = "Genome region Pi values")+
  scale_y_continuous(
    "COVID-19 cases per 100k", 
    sec.axis = sec_axis(~ ./ 10000, name = "Piw")
  )

time_series_1

# figure with correlations by genome region
cor_plot <- dat_state %>%
  pivot_longer(cols = c("genomewide_pi_state_3w",
                        "spike_pi_state_3w",
                        "orf_pi_state_3w",
                        "ntd_pi_state_3w",
                        "rbd_pi_state_3w",
                        "cov_mt_2_pi_state_3w")
  ) %>%
  mutate(name_factor = factor(name,
                              levels = c("genomewide_pi_state_3w",
                                         "spike_pi_state_3w",
                                         "orf_pi_state_3w",
                                         "ntd_pi_state_3w",
                                         "rbd_pi_state_3w",
                                         "cov_mt_2_pi_state_3w"),
                              labels = c("Genome",
                                         "Spike",
                                         "ORF NSPs 5 and 6",
                                         "S1 NTD",
                                         "S1 RBD",
                                         "CoV-MT-2"))
  ) %>%
  ggplot()+
  geom_point(aes(x = value, y = case_incidence, color = name))+
  facet_wrap(~name_factor, nrow = 6)+
  stat_cor(aes(x = value, y = case_incidence),
           method = "spearman",
           size = 4)+
  scale_color_manual(values = MetBrewer::met.brewer(n = 6,
                                                    name = "Austria"),
                     name = "",
                     labels = c("CoV-MT-2",
                                "Genomewide",
                                "S1 NTD",
                                "ORF NSPs 5 and 6",
                                "S1 RBD",
                                "Spike"))+
  theme_dth_1 +
  theme(legend.position = "none")+
  labs(title = "Spearman correlation:\nPi ~ case incidence",
       x = "Pi",
       y = "")

cor_plot


cowplot::plot_grid(time_series_1, cor_plot, nrow = 1,
                   rel_widths = c(2,1),
                   align = 'h', axis = 'tb')

# combined plot for pi
# save
png("figures and tables/Figure - pi time series and correlation.png",
    units = "in",
    width = 8.5, height = 11,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
cowplot::plot_grid(time_series_1, cor_plot, nrow = 1,
                   rel_widths = c(2,1),
                   align = 'h', axis = 'tb',
                   labels = c("A", "B"))
showtext_end()
dev.off()

# ---------------------------------------------------------------------------

# SHANNON

# Statewide time series figure, all genome regions, panel for NYC
dat_state_long <- dat_state %>%
  pivot_longer(cols = c(
    genomewide_h_state_3w,
    ntd_h_state_3w,
    cov_mt_2_h_state_3w,
    orf_h_state_3w,
    spike_h_state_3w,
    rbd_h_state_3w
  ))%>%
  mutate(name_factor = factor(name,
                              levels = c("genomewide_h_state_3w",
                                         "spike_h_state_3w",
                                         "orf_h_state_3w",
                                         "ntd_h_state_3w",
                                         "rbd_h_state_3w",
                                         "cov_mt_2_h_state_3w"),
                              labels = c("Genome",
                                         "Spike",
                                         "ORF NSPs 5 and 6",
                                         "S1 NTD",
                                         "S1 RBD",
                                         "CoV-MT-2"))
  )

time_series_1 <- 
  ggplot(data = dat_state_long) +
  geom_bar(aes(x = week, y = case_incidence,
               fill = "grey60"),
           position = "dodge", stat = "identity")+
  geom_line(aes(x = week, y = value * 10000,
                color = name),
            lwd= 1.15)+
  facet_wrap(~name_factor, nrow = 6)+
  theme_dth_1+
  theme(legend.position = "none",
        legend.background = element_rect(color = "white"),
        axis.text.x = element_text(angle = 90, vjust = 0)
  )+
  scale_color_manual(values = MetBrewer::met.brewer(n = 6,
                                                    name = "Austria"),
                     name = "",
                     labels = c("CoV-MT-2",
                                "Genomewide",
                                "S1 NTD",
                                "ORF NSPs 5 and 6",
                                "S1 RBD",
                                "Spike"))+
  scale_x_date(date_labels = "%b %y",
               date_breaks = "1 month",
               limits = c(ymd("2023-01-01"),
                          ymd("2025-04-20"))
  )+
  scale_fill_manual(values = c("grey60" = "grey60"),
                    labels = c("Cases\nper 100k"),
                    name = ""
  )+
  labs(x = "",
       title = "Genome region Shannon's H values")+
  scale_y_continuous(
    "COVID-19 cases per 100k", 
    sec.axis = sec_axis(~ ./ 10000, name = "Shannon's H")
  )

time_series_1

# figure with correlations by genome region
cor_plot <- dat_state %>%
  pivot_longer(cols = c("genomewide_h_state_3w",
                        "spike_h_state_3w",
                        "orf_h_state_3w",
                        "ntd_h_state_3w",
                        "rbd_h_state_3w",
                        "cov_mt_2_h_state_3w")
  ) %>%
  mutate(name_factor = factor(name,
                              levels = c("genomewide_h_state_3w",
                                         "spike_h_state_3w",
                                         "orf_h_state_3w",
                                         "ntd_h_state_3w",
                                         "rbd_h_state_3w",
                                         "cov_mt_2_h_state_3w"),
                              labels = c("Genome",
                                         "Spike",
                                         "ORF NSPs 5 and 6",
                                         "S1 NTD",
                                         "S1 RBD",
                                         "CoV-MT-2"))
  ) %>%
  ggplot()+
  geom_point(aes(x = value, y = case_incidence, color = name))+
  facet_wrap(~name_factor, nrow = 6)+
  stat_cor(aes(x = value, y = case_incidence),
           method = "spearman",
           size = 4)+
  scale_color_manual(values = MetBrewer::met.brewer(n = 6,
                                                    name = "Austria"),
                     name = "",
                     labels = c("CoV-MT-2",
                                "Genomewide",
                                "S1 NTD",
                                "ORF NSPs 5 and 6",
                                "S1 RBD",
                                "Spike"))+
  theme_dth_1 +
  theme(legend.position = "none")+
  labs(title = "Spearman correlation:\nShannon's H ~ case incidence",
       x = "Shannon's H",
       y = "")

cor_plot


cowplot::plot_grid(time_series_1, cor_plot, nrow = 1,
                   rel_widths = c(2,1),
                   align = 'h', axis = 'tb')

# combined plot for pi
# save
png("figures and tables/Figure - H time series and correlation.png",
    units = "in",
    width = 8.5, height = 11,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
cowplot::plot_grid(time_series_1, cor_plot, nrow = 1,
                   rel_widths = c(2,1),
                   align = 'h', axis = 'tb',
                   labels = c("A", "B"))
showtext_end()
dev.off()

# ----------------------------------------------------------------------------

# ----------------------------------------------------
# CORRELATION - LAG AND LEAD TIME
# ----------------------------------------------------

# correllation with hospital admissions

pi_lag <- lag_lead_function_hosp(dataframe = dat_state,
                                 columns = c("genomewide_pi_state_3w",
                                             "spike_pi_state_3w",
                                             "orf_pi_state_3w",
                                             "ntd_pi_state_3w",
                                             "rbd_pi_state_3w",
                                             "cov_mt_2_pi_state_3w"),
                                 lags = seq(1:10))


# plot
lag_plot <- 
  ggplot(data = pi_lag,
       aes(x = lag, y = spearman_cor, color = name))+
  geom_point(shape = 1, size = 2)+
  geom_point(shape = 17,
             data = pi_lag %>%
               group_by(name) %>%
               filter(spearman_cor == max(spearman_cor)),
             size = 3
  )+
  geom_line()+
  ylim(0,1)+
  theme_dth_1+
  theme(legend.position = "bottom")+
  scale_color_manual(values = MetBrewer::met.brewer(n = 6,
                                                    name = "Austria"),
                     name = "",
                     labels = c("CoV-MT-2",
                                "Genomewide",
                                "S1 NTD",
                                "ORF NSPs 5 and 6",
                                "S1 RBD",
                                "Spike"))+
  geom_vline(xintercept = 0)+
  scale_x_continuous(breaks = seq(-5, 5, by = 1),
                     limits = c(-5,5))+
  labs(x = "Lead time (weeks)",
       y = "Spearman correlation\ncoefficient",
       title= "Hospital incidence ~ Pi (wastewater) for\nstatewide data")+
  annotate("text",
           x = -3, y = 0.09,
           label = "signal lags incidence")+
  annotate("text",
           x = 3, y = 0.09,
           label = "signal leads incidence")

# repeat for shannon

lag_h <- lag_lead_function_hosp(dataframe = dat_state,
                                columns = c("genomewide_h_state_3w",
                                            "spike_h_state_3w",
                                            "orf_h_state_3w",
                                            "ntd_h_state_3w",
                                            "rbd_h_state_3w",
                                            "cov_mt_2_h_state_3w"),
                                lags = seq(1:10))

lag_plot_2 <- ggplot(data = lag_h ,
       aes(x = lag, y = pearson_cor, color = name))+
  geom_point(shape = 1, size = 2)+
  geom_point(shape = 17,
             data = lag_h %>%
               group_by(name) %>%
               filter(pearson_cor == max(pearson_cor)),
             size = 3
  )+
  geom_line()+
  ylim(0,1)+
  theme_dth_1+
  theme(legend.position = "bottom")+
  scale_color_manual(values = MetBrewer::met.brewer(n = 6,
                                                    name = "Austria"),
                     name = "",
                     labels = c("CoV-MT-2",
                                "Genomewide",
                                "S1 NTD",
                                "ORF NSPs 5 and 6",
                                "S1 RBD",
                                "Spike"))+
  geom_vline(xintercept = 0)+
  scale_x_continuous(breaks = seq(-5, 5, by = 1),
                     limits = c(-5,5))+
  labs(x = "Lead time (weeks)",
       y = "Spearman correlation\ncoefficient",
       title= "Hospital incidence ~ Shannon H (wastewater) for\nstatewide data")+
  annotate("text",
           x = -3, y = 0.09,
           label = "signal lags incidence")+
  annotate("text",
           x = 3, y = 0.09,
           label = "signal leads incidence")

# add  cases
pi_lag_2 <- lag_lead_function(dataframe = dat_state,
                            columns = c("genomewide_pi_state_3w",
                                        "spike_pi_state_3w",
                                        "orf_pi_state_3w",
                                        "ntd_pi_state_3w",
                                        "rbd_pi_state_3w",
                                        "cov_mt_2_pi_state_3w"),
                            lags = seq(1:10))


# plot
lag_plot_case <- 
  ggplot(data = pi_lag_2,
         aes(x = lag, y = spearman_cor, color = name))+
  geom_point(shape = 1, size = 2)+
  geom_point(shape = 17,
             data = pi_lag_2 %>%
               group_by(name) %>%
               filter(spearman_cor == max(spearman_cor)),
             size = 3
  )+
  geom_line()+
  ylim(0,1)+
  theme_dth_1+
  theme(legend.position = "bottom")+
  scale_color_manual(values = MetBrewer::met.brewer(n = 6,
                                                    name = "Austria"),
                     name = "",
                     labels = c("CoV-MT-2",
                                "Genomewide",
                                "S1 NTD",
                                "ORF NSPs 5 and 6",
                                "S1 RBD",
                                "Spike"))+
  geom_vline(xintercept = 0)+
  scale_x_continuous(breaks = seq(-5, 5, by = 1),
                     limits = c(-5,5))+
  labs(x = "Lead time (weeks)",
       y = "Spearman correlation\ncoefficient",
       title= "Case incidence ~ Pi (wastewater) for\nstatewide data")+
  annotate("text",
           x = -3, y = 0.09,
           label = "signal lags incidence")+
  annotate("text",
           x = 3, y = 0.09,
           label = "signal leads incidence")


# repeat for shannon

lag_h_2 <- lag_lead_function(dataframe = dat_state,
                           columns = c("genomewide_h_state_3w",
                                       "spike_h_state_3w",
                                       "orf_h_state_3w",
                                       "ntd_h_state_3w",
                                       "rbd_h_state_3w",
                                       "cov_mt_2_h_state_3w"),
                           lags = seq(1:10))

lag_plot_2_case <- ggplot(data = lag_h_2 ,
                     aes(x = lag, y = pearson_cor, color = name))+
  geom_point(shape = 1, size = 2)+
  geom_point(shape = 17,
             data = lag_h_2 %>%
               group_by(name) %>%
               filter(pearson_cor == max(pearson_cor)),
             size = 3
  )+
  geom_line()+
  ylim(0,1)+
  theme_dth_1+
  theme(legend.position = "bottom")+
  scale_color_manual(values = MetBrewer::met.brewer(n = 6,
                                                    name = "Austria"),
                     name = "",
                     labels = c("CoV-MT-2",
                                "Genomewide",
                                "S1 NTD",
                                "ORF NSPs 5 and 6",
                                "S1 RBD",
                                "Spike"))+
  geom_vline(xintercept = 0)+
  scale_x_continuous(breaks = seq(-5, 5, by = 1),
                     limits = c(-5,5))+
  labs(x = "Lead time (weeks)",
       y = "Spearman correlation\ncoefficient",
       title= "Case incidence ~ Shannon H (wastewater) for\nstatewide data")+
  annotate("text",
           x = -3, y = 0.09,
           label = "signal lags incidence")+
  annotate("text",
           x = 3, y = 0.09,
           label = "signal leads incidence")

# extract legend function
mylegend <- g_legend(lag_plot)



# save
png("figures and tables/Figure - lag plot.png",
    units = "in",
    width = 10, height = 11,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
ggarrange(ggarrange(lag_plot_case+ theme(legend.position = "none"),
                    lag_plot_2_case + theme(legend.position = "none"),
                    nrow = 1,
                    labels = c("A", "B")),
          ggarrange(lag_plot+ theme(legend.position = "none"),
                    lag_plot_2 + theme(legend.position = "none"),
                    nrow = 1,
                    labels = c("C", "D")),
          mylegend, 
          nrow=3,
          heights=c(5,5, 1))
dev.off()

# -----------------------------------------------------------------------------

# ------------------------------------------------
# VARIANT COUNT ANALYSIS
# ------------------------------------------------

# New figure
ggplot(data = dat_state)+
  geom_bar(data = dat_state,
           aes(x = week,
               y = case_incidence),
           position = "stack",
           stat = "identity",
           fill = "grey60")+
  geom_line(aes(x = week,
                y = ntd_pi_state_3w*10000),
            lwd = 1
  )+
  geom_line(aes(x = week,
                n_variants_no_thresh_3w_mean*3),
            color = "yellow",
            lwd = 1)+
  geom_line(aes(x = week,
                n_variants_5_3w_mean*10),
            color = "dodgerblue",
            lwd = 1)

# basic, smoothed values by sewershed for now
# statewide reported  variant counts
var_plot_no_thresh <- 
  ggplot()+
  geom_bar(data = dat_state,
           aes(x = week,
               y = case_incidence,
               fill = "grey60"),
           position = "stack",
           stat = "identity")+
  geom_line(data = dat_state,
            aes(x = week, y= n_variants_no_thresh_3w_mean*6,
                color = "darkblue"),
            linewidth = 1
            )+
  theme_dth_1+
  theme(legend.position = "bottom",
        legend.background = element_blank())+
  labs(
       x = "")+
  scale_y_continuous(
    "COVID-19 cases per 100k", 
    sec.axis = sec_axis(~ ./ 6, name = "Lineage count")
  )+
  scale_fill_manual(values = c("grey60" = "grey60"),
                    label = "Cases",
                    name = "")+
  scale_color_manual(values = c("darkblue" = "darkblue"),
                     label = "Freyja lineages",
                     name = "")

var_plot_no_thresh

var_plot_thresh <- 
  ggplot()+
  geom_bar(data = dat_state,
           aes(x = week,
               y = case_incidence,
               fill = "grey60"),
           position = "stack",
           stat = "identity")+
  geom_line(data = dat_state,
            aes(x = week, y= n_variants_5_3w_mean*20,
                color = "darkblue"),
            linewidth = 1
  )+
  theme_dth_1+
  theme(legend.position = "bottom",
        legend.background = element_blank())+
  labs(
       x = "")+
  scale_y_continuous(
    "COVID-19 cases per 100k", 
    sec.axis = sec_axis(~ ./ 20, name = "Lineage count")
  )+
  scale_fill_manual(values = c("grey60" = "grey60"),
                    label = "Cases",
                    name = "")+
  scale_color_manual(values = c("darkblue" = "darkblue"),
                     label = "Freyja lineages",
                     name = "")

var_plot_thresh

ggarrange(var_plot_no_thresh, var_plot_thresh,
          nrow = 2, common.legend = TRUE, legend="bottom")

freyja_cor_1 <- 
  ggplot(data = dat_state,
       aes(x = n_variants_no_thresh_3w_mean,
           y = case_incidence))+
  geom_point(color = "darkblue")+
  stat_cor(method = "spearman"
           )+
  theme_dth_1+
  labs(x = "Mean number of Freyja lineages",
       y = "COVID-19 cases per 100k")
freyja_cor_1

freyja_cor_2 <- 
  ggplot(data = dat_state,
       aes(x = n_variants_5_3w_mean,
           y = case_incidence))+
  geom_point(color = "darkblue")+
  stat_cor(method = "spearman")+
  theme_dth_1+
  labs(x = "Mean number of Freyja lineages",
       y = "COVID-19 cases per 100k")

# extract legend
mylegend<-g_legend(var_plot_no_thresh)

# create each row of the panel
row_1 <- plot_grid(
  var_plot_no_thresh+ theme(legend.position = "none"),
  freyja_cor_1,
  nrow = 1,
  rel_widths = c(2,1),
  align = 'h', axis = 'tb',
  labels = c("A", "B")
)
row_1

# add common title
title <- ggdraw() + 
  draw_label("Mean count of Freyja lineages (no threshold cutoff)", 
             fontface='bold')
row_1_final <- plot_grid(title, row_1, ncol=1, rel_heights=c(0.1, 1))


row_2 <- plot_grid(
  var_plot_thresh+ theme(legend.position = "none"),
  freyja_cor_2,
  nrow = 1,
  rel_widths = c(2,1),
  align = 'h', axis = 'tb',
  labels = c("C", "D")
)

row_2

# add title
title <- ggdraw() + 
  draw_label("Mean count of Freyja lineages (5% cutoff)", 
             fontface='bold')
row_2_final <- plot_grid(title, row_2, ncol=1, rel_heights=c(0.1, 1))


# save
png("figures and tables/Figure - variant counts.png",
    units = "in",
    width = 11, height = 9,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
plot_grid(row_1_final, row_2_final,
          nrow = 3,
          mylegend,
          rel_heights = c(5,5,1)
)

dev.off()

# ------------------------------------------------
# PREDICTIVE ANALSYSIS - CASES
# ------------------------------------------------

# STATE MODEL 

# remove nas
s <- dat_state %>%
  filter(!is.na(case_incidence)) %>%
  filter(!is.na(ntd_pi_state_3w)) %>%
  filter(!is.na(mean_sars2_conc_state_3w)) %>%
  filter(!is.na(depth_state_3w)) %>%
  filter(!is.na(n_variants_no_thresh_3w_mean))

# model to fit
model_gaussian <- glm(
  case_incidence ~ 
    scale(ntd_pi_state_3w)+
    scale(mean_sars2_conc_state_3w)+
    scale(depth_state_3w)+
    scale(n_variants_no_thresh_3w_mean),
  family = gaussian(),
  data = s
)

# create table
state_table_gaus <- model_summary_function(
  dataframe = s,
  outcome = "case_incidence",
  model = model_gaussian,
  glmer_option = "No"
) %>%
  mutate(group = "Full model")

# ntd only
model_gaussian <- glm(
  case_incidence ~ 
    scale(ntd_pi_state_3w),
  family = gaussian(),
  data = s
) 

# table
state_table_ntd <- model_summary_function(
  dataframe = s,
  outcome = "case_incidence",
  model = model_gaussian,
  glmer_option = "No"
)%>%
  mutate(group = "S1 NTD only model")

# conc only
model_gaussian <- glm(
  case_incidence ~ 
    scale(mean_sars2_conc_state_3w),
  family = gaussian(),
  data = s
)

# table
state_table_conc <- model_summary_function(
  dataframe = s,
  outcome = "case_incidence",
  model = model_gaussian,
  glmer_option = "No"
)%>%
  mutate(group = "Concentration only model")

# CREATE TABLES FOR PAPER

# state model table
table <- full_join(
  state_table_gaus, state_table_conc, by = c("variable")
) %>%
  left_join(state_table_ntd, by = c("variable"))

table <- bind_rows(state_table_gaus,
                   state_table_conc,
                   state_table_ntd) %>%
  dplyr::select(group, everything())
table$`Pr(>|t|)` <- ifelse(table$`Pr(>|t|)` < 0.0001,
                           "<0.0001",
                           table$`Pr(>|t|)`)

# edit estimate column
table$variable<- c(
  "Intercept",
  "S1 NTD Pi",
  "SARS-CoV-2 concentration",
  "Mean sample depth of read",
  "Number of Freyja variants (no threshold)",
  "n",
  "Efron R2",
  "AIC",
  "Intercept",
  "SARS-CoV-2 concentration",
  "n",
  "Efron R2",
  "AIC",
  "Intercept",
  "S1 NTD Pi",
  "n",
  "Efron R2",
  "AIC"
)


# make it an ft table
t <- table_as_flex_function(dataframe = table,
                       title = "Tale: Generalized liner model results for S1 NTD")

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
             path = paste("figures and tables/",
                          "Table - state model results.docx")
)
