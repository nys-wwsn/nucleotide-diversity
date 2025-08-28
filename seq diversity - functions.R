# Functions for sequencing diversity project

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

# calculation correlations with confidence intervals
cor_ci_function <- function(dataframe,
                            Region,
                            predictor_name,
                            predictor_group){
  
  # docstring information
  #' Output a Spearman correlation with 95% CI
  #' 
  #' @description This function outputs Spearman correlation values for each
  #' predictor of interest with cases and hospitalizations for COVID and 
  #' provides 95% LL and UL for each.
  #' 
  #' @param dataframe dataframe. Object with data for correlation tests.
  #' @param Region character. If a genome region is the predictor, the name
  #' of the region
  #' @param predictor_name character. Variable name as written in dataframe
  #' @param predictor_group character. Pi or H or blank
  #' 
  #' @return Dataframe with correlations.
  #' 
  
  # add variables to the dataframe
  #dataframe$comparison_variable <- dataframe[[comparison_variable]]
  dataframe$predictor_name <- dataframe[[predictor_name]]
  
  set.seed(100)
  # case correlation
  case_cor <-  ci_cor(
    x = as.data.frame(cbind(dataframe$predictor_name, dataframe$case_incidence)),
    probs = c(0.025, 0.95),
    method = "spearman",
    type = "bootstrap",
    R = 1000
  )
  spearman_statistic <- case_cor$estimate
  ll <- case_cor$interval[1]
  ul <- case_cor$interval[2]
  
  # hosp correlation
  hosp_cor <-  ci_cor(
    x = as.data.frame(cbind(dataframe$predictor_name, dataframe$hosp_incidence)),
    probs = c(0.025, 0.95),
    method = "spearman",
    type = "bootstrap",
    R = 1000
  )
  spearman_statistic_h <- hosp_cor$estimate
  ll_h <- hosp_cor$interval[1]
  ul_h <- hosp_cor$interval[2]
  
  observations <- length(dataframe$predictor_name)
  row_df <- as.data.frame(cbind(predictor_group,
                                Region, 
                                paste(
                                  signif(spearman_statistic, 3), 
                                  paste("[", signif(ll, 3), ", ", 
                                        signif(ul, 3), "]", sep = ""),
                                  sep = " "), 
                                paste(
                                  signif(spearman_statistic_h, 3), 
                                  paste("[", signif(ll_h, 3), ", ", 
                                        signif(ul_h, 3), "]", sep = ""),
                                  sep = " "),
                                observations
  )
  ) %>%
    rename(
      "Case incidence, 95% CI" = V3,
      "Hospitalization incidence, 95% CI" = V4)
  return(row_df)
}

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

# granger causality function - outputs F test for time series prediction test
granger_casuality_function <- function(dataframe,
                                       y,
                                       x,
                                       lag,
                                       group,
                                       outcome){
  # docstring information
  #' Export results of Granger causality test
  #' 
  #' @description This function compares two time series using the grangertest
  #' function from the lmtest package. Each time series is adjusted to be
  #' stationary using a lag of 1 week.
  #' 
  #' @param dataframe dataframe. Object with the two time series to test.
  #' @param y character. outcome variable, the one being predicted by the other
  #' time series
  #' @param x character. predictor variable
  #' @param lag integer. Number of lags for optimal correlation
  #' @param group character. Character string for predictor name.
  #' @param outcome character. Character string for outcome name.
  #' 
  #' @return Dataframe with output of the grangertest function including
  #' f test statistic and p value
  #' 
  
  # create the lagged variables using diff
  # case incidence is 1 week based on pearson correlation coefficients
  # hosp incidence is 2 weeks based on pearson correlation coefficients
  d_x <- diff(dataframe[[x]])
  d_y <- diff(dataframe[[y]])
  
  # output results
  g1 <- as.data.frame(grangertest(d_y ~ d_x, order = lag))
  colnames(g1) <- c("ResDF", "DF", "F statistic", "P value")
  
  # keep first row f statistic and p value
  g1 <- g1 %>%
    mutate(group = group,
           outcome = outcome,
           test = "x predicts y") %>%
    select(group, outcome, test, `F statistic`, `P value`) %>%
    na.omit()
  
  g2 <- as.data.frame(grangertest(d_x ~ d_y, order = lag))
  colnames(g2) <- c("ResDF", "DF", "F statistic", "P value")
  
  # keep first row f statistic and p value
  g2 <- g2 %>%
    mutate(group = group,
           outcome = outcome,
           test = "y predicts x") %>%
    select(group, outcome, test, `F statistic`, `P value`) %>%
    na.omit()
  g <- bind_rows(g1, g2)
  return(g)
}

# function for making time series/correlation plots for depth subsamples
time_cor_plot_function <- function(dataframe,
                                   depth_start,
                                   depth_end,
                                   title){
  
  set.seed(23)
  sample_3_3w <- dat_sewershed  %>%
    filter(depth >= depth_start & depth <= depth_end) %>%
    sample_n(1000) %>%
    arrange(week) %>%
    group_by(week)%>%
    summarise(
      # statewide pi
      genomewide_pi_state_3w = weighted.mean(x = genomewide_pi_ma3, 
                                             w = population_served, 
                                             na.rm = TRUE),
      ntd_pi_state_3w = weighted.mean(x = ntd_pi_ma3, 
                                      w = population_served, 
                                      na.rm = TRUE)
    )%>%
    ungroup() %>%
    left_join(cases, by = c("week"))
  
  # correlation plot
  p_cor <- 
    ggplot(data = sample_3_3w, aes(x = ntd_pi_state_3w,
                                   y = case_incidence))+
    geom_point()+
    stat_cor(method = "spearman")+
    theme_dth_1+
    labs(x = "S1 NTD Pi",
         y = "",
         title = "")
  
  # time series plot
  p_time <- 
    ggplot(data = sample_3_3w)+
    geom_bar(data = dat_state %>%
               filter(group == "State"),
             aes(x = week,
                 y = case_incidence,
                 fill = "grey60"),
             position = "dodge",
             stat = "identity")+
    geom_line(aes(x = week,
                  y = ntd_pi_state_3w*10000,
                  color = "darkblue"),
              lwd = 1.5)+
    theme_dth_1+
    labs(x = "",
         title = "")+
    scale_y_continuous(
      "COVID-19 cases\nper 100k", 
      sec.axis = sec_axis(~ ./ 10000, name = "S1 NTD Pi"))+
    scale_color_manual(values = c("darkblue"),
                       labels = c("Pi"),
                       name = "")+
    scale_fill_manual(values = c("grey60"),
                      labels = "Cases",
                      name = "")+
    theme(legend.background = element_blank())+
    theme(legend.position = "bottom")
  
  # make the panel
  panel <- ggarrange(p_time+
                       theme(legend.position = "none"), p_cor,
                     align = "h",
                     widths = c(2,1))
  
  # add title
  panel_final <- annotate_figure(panel, title)
  
  # return the panel plot
  return(panel_final)
}
# -----------------------------------------------------------------------------