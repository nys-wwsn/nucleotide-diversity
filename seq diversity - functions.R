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
library(DescTools)


# --------------------------------------
# FUNCTIONS
# --------------------------------------

# DIVERSITY CALCULATION FUNCTIONS


# Functions

# calculate pi function

# Functions

# calculate diversity functions -> pi and shannon
calculate_diversity_function <- function(input, pass_option){
  
  # docstring information
  #' Calculate nucleotide and Shannon diversity
  #' 
  #' @description This function calculates nucleotide diversity for wastewater
  #' samples and Shannon's H. 
  #' 
  #' @param input dataframe. Freyja variants file from the Freyja variant
  #' pipeline.
  #' @param pass_option character string. Whether to include all records or only
  #' Freyja pass records. Options are "all" or "pass only".
  #' 
  #' @details Nucleotide diversity is calculated using the following equation: 
  #' \deqn{\pi_{s} =( \frac{n}{n-1} ) (1 - \sum{f^{2}})}
  #' Where \eqn{n} is the total number of reads spanning that position, \eqn{f}
  #'  is the frequency of a variant, and the sum is over all variants at that 
  #'  position.
  #' 
  #' Shannon's H is calculated in two parts. The proportion of each base/allele
  #' observed is calculated using the following equation:
  #' \deqn{H_{s} = \sum{(log(prop) * prop) * -1}}
  #' 
  #' Shannon's H full equation is:
  #' \deqn{$prop = \frac{\text{frequency of allele}}{\text{total alleles observed}} }
  #' 
  #' @usage calculate_diversity_function(input, pass_option)
  #' @return Dataframe with nucleotide diversity, Shannon's H, depth of read
  #' for the one sample.
  
  #' @note This function is wrapped inside the 
  #' "parallel_diversity_per_base_function".
  #' @examples
  #' 
  
  library(dplyr)
  library(readr)
  
  # filter for reads that pass
  if(pass_option == "pass only"){
    input <- input %>% filter(PASS == TRUE)
  } else if(pass_option == "all"){
    input <- input
  }
  
  # select reference counts and alternate counts
  reference_counts <- input %>% 
    dplyr::select(POS, REF_DP) %>% 
    distinct(POS, .keep_all = TRUE) %>% 
    rename(k = REF_DP) %>%
    mutate(freq = k)
  alternate_counts <- input %>% 
    dplyr::select(POS, ALT_DP) %>% 
    rename(k = ALT_DP) %>%
    mutate(freq = k)
  
  # bind together and add row ids
  input <- bind_rows(reference_counts, alternate_counts) %>%
    arrange(POS)
  input$row_id <- 1:nrow(input)
  
  # create depth dataset
  depths <- input %>% group_by(POS) %>% summarise(sum_k = sum(k))
  
  # proportion that was observed by base
  prop <- function(freq){
    total <- sum(freq)
    prop_value <- freq/total
    return(prop_value)
  }
  
  input <- input %>%
    group_by(POS) %>%
    mutate(
      prop = prop(freq)
    )
  
  # shannon index
  h <- function(prop){
    h_value <- sum((log(prop) * prop ) * -1)
  }
  
  # pi function
  pi <- function(s) {
    n <- sum(s)
    pi_value <- (1 - sum((s / n) ^ 2)) * (n / (n - 1))
    return(pi_value)
  }
  
  # calculate pi
  result <- input %>%
    group_by(POS) %>%
    summarise(pi = pi(k),
              h = h(prop)) %>%
    left_join(depths, by = "POS") %>%
    rename(depth = sum_k)
  
  # impute 0 for NA shannons
  result$h <- ifelse(is.na(result$h),
                     0,
                     result$h)
  
  # save output
  return(result)
}

# windowed diversity function
windowed_diversity_function <- function(sample_results, 
                                        genome_size, 
                                        window_size, 
                                        window_step){
  
  # docstring information
  #' Calculate windowed Pi and Shannon values across groups of base pairs
  #' 
  #' @description This function uses the Pi per base and Shannon per base data
  #' output from the `calculate_diversity_function` to estimate a windowed
  #' value for each with sliding windows that are adjustable.
  #' 
  #' @param sample_results dataframe. Pi and Shannon per base output from
  #' `calculate_diversity_function`.
  #' @param genome_size integer. Length of the genome being analyzed. For 
  #' SARS-CoV-2, we use 30,000 base pairs.
  #' @param window_size numeric. Length of the window size to use. For example,
  #' our study uses 1000 bp windows.
  #' @param window_step numeric. Length of the step between window calculations
  #' for the function to slide over. For example, a 1000 bp window with 100 bp
  #' step would mean that a windowed Pi is calculated for every 1000 bps moving 
  #' ahead 100 bps for each calculation
  #' 
  #' @details Windowed Pi is calculated using the following equation:
  #' 
  #'\deqn{\pi_{w} = \frac{1}{L} \sum{\pi_s}$}
  #'
  #'Where \eqn{L} is the window size in base pairs. Positions with no variation in 
  #'the sample are considered to have \eqn{\pi_s} = 0. 
  #' 
  #' Windowed Shannon's H is calculated using a similar equation:
  #' \deqn{H_{w} = \frac{1}{L} \sum{H_{s}}}
  #' 
  #' @return Dataframe with windowed nucleotide diversity, windowed Shannon's H, 
  #' windowed mean depth of read for the one sample.
  
  #' @note This function is wrapped inside the 
  #' `parallel_mean_diversity_function`.
  #' @examples
  #' 
  
  library(dplyr)
  library(readr)
  library(tidyr)
  
  genome_size <- genome_size
  window_size <- window_size
  window_step <- window_step
  
  pi <- sample_results
  
  # stop running the code if the max position is greater than the genome size
  if (nrow(pi) > 0) {
    stopifnot(max(pi$POS) <= genome_size)
  }
  
  # empty list object
  records <- list()
  
  window_start <- 0
  window_end <- window_start + window_size
  window_center <- window_size %/% 2 + window_start
  
  # windowed pi and shannon and mean depth
  while (window_end < genome_size) {
    data <- pi %>% filter(POS >= window_start & POS < window_end)
    windowed_pi <- sum(data$pi, na.rm = TRUE) / window_size
    windowed_h = sum(data$h, na.rm = TRUE) / window_size
    mean_depth <- mean(data$depth, na.rm = TRUE)
    if (is.na(mean_depth)) {
      mean_depth <- 0
    }
    records[[length(records) + 1]] <- c(window_start, window_center, window_end, windowed_pi, windowed_h, mean_depth)
    window_start <- window_start + window_step
    window_end <- window_start + window_size
    window_center <- window_size %/% 2 + window_start
  }
  
  result <- as.data.frame(do.call(rbind, records))
  colnames(result) <- c("start", "center", "end", "avg_pi","avg_h", "avg_depth")
  
  # save down the new file
  return(result)
  
}

# for loop 1 function - pi ber base per file in parallel
# for each file (each file is a sample)
# you ned to set up your directory for saving the output files

parallel_diversity_per_base_function <- function(data_dir,
                                                 save_path,
                                                 pass_option){
  
  # docstring information
  #' Calculate diversity measures per base in parallel
  #' 
  #' @description Estimate Pi and Shannon per base in parallel for multiple 
  #' wastewater genetic sequencing samples.
  #' 
  #' @param data_dir file path. Directory path for where the genetic data
  #' are saved.
  #' @param save_path file path. Directory path to save the output.
  #' @param pass_option character. Options are "all" or "pass only". Indicates
  #' if the Pi per base / shannon per base should use all records or only those
  #' that pass Freyja QC.
  #' 
  #' @details This function uses the `calculate_diversity_function` 
  #' to estimate Pi and Shannon per base in parallel for multiple wastewater
  #' genetic sequencing data files. The function is designed to query each file
  #' separately, where each file is a sample.
  #' 
  #' @return Tab separated file (TSV) for each sample.
  #' @examples
  #' parallel_diversity_per_base_function(data_dir = "variants/",
  #' save_path = "data/pi/",
  #' pass_option = "pass only")
  
  library(foreach)
  library(doParallel)
  
  parallel::detectCores()
  n.cores <- parallel::detectCores()/2
  
  #create the cluster
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  
  #check cluster definition (optional)
  print(my.cluster)
  
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  #check if it is registered (optional)
  foreach::getDoParRegistered()
  
  #how many workers are available? (optional)
  foreach::getDoParWorkers()
  
  # directory
  data_dir <- data_dir
  # files to iterate over
  list_of_files <- list.files(path = data_dir,
                              recursive = TRUE,
                              pattern = "\\.tsv$",
                              full.names = TRUE)
  
  # run process in parallel
  foreach(
    i = list_of_files
  ) %dopar% {
    
    library(dplyr)
    library(ggplot2)
    library(stringr)
    library(tidyr)
    library(purrr)
    library(readr)
    
    # functions
    calculate_diversity_function <- function(input, pass_option){
      library(dplyr)
      library(readr)
      
      # filter for reads that pass
      if(pass_option == "pass only"){
        input <- input %>% filter(PASS == TRUE)
      } else if(pass_option == "all"){
        input <- input
      }
      
      # select reference counts and alternate counts
      reference_counts <- input %>% 
        dplyr::select(POS, REF_DP) %>% 
        distinct(POS, .keep_all = TRUE) %>% 
        rename(k = REF_DP) %>%
        mutate(freq = k)
      alternate_counts <- input %>% 
        dplyr::select(POS, ALT_DP) %>% 
        rename(k = ALT_DP) %>%
        mutate(freq = k)
      
      # bind together and add row ids
      input <- bind_rows(reference_counts, alternate_counts) %>%
        arrange(POS)
      input$row_id <- 1:nrow(input)
      
      # create depth dataset
      depths <- input %>% group_by(POS) %>% summarise(sum_k = sum(k))
      
      # proportion that was observed by base
      prop <- function(freq){
        total <- sum(freq)
        prop_value <- freq/total
        return(prop_value)
      }
      
      input <- input %>%
        group_by(POS) %>%
        mutate(
          prop = prop(freq)
        )
      
      # shannon index
      
      h <- function(prop){
        h_value <- sum((log(prop) * prop ) * -1)
      }
      
      # pi function
      pi <- function(s) {
        n <- sum(s)
        pi_value <- (1 - sum((s / n) ^ 2)) * (n / (n - 1))
        return(pi_value)
      }
      
      # calculate pi
      result <- input %>%
        group_by(POS) %>%
        summarise(pi = pi(k),
                  h = h(prop)) %>%
        left_join(depths, by = "POS") %>%
        rename(depth = sum_k)
      
      # impute 0 for NA shannons
      result$h <- ifelse(is.na(result$h),
                         0,
                         result$h)
      # save output
      return(result)
    }
    
    # path to save files
    path <- save_path
    
    # file for the input of the function
    df_i <- i %>%
      set_names() %>%  
      map_df(read_tsv, .id = "file_name") %>% 
      mutate(name_rev = sapply(str_split((sapply((str_split(file_name, "/")), tail, 1)), "-"), head, 1)) %>%
      separate(name_rev, into = c('date', 'cdc_id'), sep = 8) %>%
      #mutate(cdc_id = str_split(cdc_id, "_", 2))
      separate_wider_delim(cdc_id, delim = "_", names = c("cdc_id", "drop"))
    
    # run the function
    pi_file <- calculate_diversity_function(input = df_i,
                                            pass_option = pass_option)
    
    # save the output in a folder
    write_tsv(pi_file, file = paste(path, head(df_i$date, 1), head(df_i$cdc_id, 1), ".tsv", sep =""))
  }
  
  
  # stop cluster
  parallel::stopCluster(cl = my.cluster)
}

# for loop 2 function - windowed pi per file
parallel_windowed_diversity_function <- function(data_dir,
                                                 save_path){
  
  # docstring information
  #' Calculate windowed diversity in parallel
  #' 
  #' @description Estimate mean Pi and mean Shannon per window in parallel for 
  #' multiple wastewater genetic sequencing samples.
  #' 
  #' @param data_dir file path. Directory path for where the genetic data
  #' are saved.
  #' @param save_path file path. Directory path to save the output.
  #' 
  #' @details This function uses the `windowed_diversity_function` 
  #' to estimate windowed Pi and windowed Shannon in parallel for multiple 
  #' wastewater genetic sequencing data files. The function is designed to query 
  #' each file separately, where each file is a sample.
  #' 
  #' @return Tab separated file (TSV) for each sample.
  #' @examples
  #' parallel_windowed_diversity_function(data_dir = "data/pi/",
  #' save_path = "data/windowed_pi/")
  #' 
  library(foreach)
  library(doParallel)
  
  parallel::detectCores()
  n.cores <- parallel::detectCores()/2
  
  #create the cluster
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  
  #check cluster definition (optional)
  print(my.cluster)
  
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  #check if it is registered (optional)
  foreach::getDoParRegistered()
  
  #how many workers are available? (optional)
  foreach::getDoParWorkers()
  
  # files to iterate over
  list_of_files <- list.files(path = data_dir,
                              recursive = TRUE,
                              pattern = "\\.tsv$",
                              full.names = TRUE)
  
  # run process in parallel
  foreach(
    i = list_of_files
  ) %dopar% {
    
    library(dplyr)
    library(ggplot2)
    library(stringr)
    library(tidyr)
    library(purrr)
    library(readr)
    
    # functions
    windowed_diversity_function <- function(sample_results, 
                                            genome_size, 
                                            window_size, 
                                            window_step){
      
      library(dplyr)
      library(readr)
      library(tidyr)
      
      genome_size <- genome_size
      window_size <- window_size
      window_step <- window_step
      
      pi <- sample_results
      
      # stop running the code if the max position is greater than the genome size
      if (nrow(pi) > 0) {
        stopifnot(max(pi$POS) <= genome_size)
      }
      
      # empty list object
      records <- list()
      
      window_start <- 0
      window_end <- window_start + window_size
      window_center <- window_size %/% 2 + window_start
      
      # windowed pi and shannon and mean depth
      while (window_end < genome_size) {
        data <- pi %>% filter(POS >= window_start & POS < window_end)
        windowed_pi <- sum(data$pi, na.rm = TRUE) / window_size
        windowed_h = sum(data$h, na.rm = TRUE) / window_size
        mean_depth <- mean(data$depth, na.rm = TRUE)
        if (is.na(mean_depth)) {
          mean_depth <- 0
        }
        records[[length(records) + 1]] <- c(window_start, window_center, window_end, windowed_pi, windowed_h, mean_depth)
        window_start <- window_start + window_step
        window_end <- window_start + window_size
        window_center <- window_size %/% 2 + window_start
      }
      
      result <- as.data.frame(do.call(rbind, records))
      colnames(result) <- c("start", "center", "end", "avg_pi","avg_h", "avg_depth")
      
      # save down the new file
      return(result)
      
    }
    
    # path to save files
    path <- save_path
    
    # file for the input of the function
    df_i <- i %>%
      set_names() %>%  
      map_df(read_tsv, .id = "file_name") %>% 
      mutate(name_rev = sapply(str_split((sapply((str_split(file_name, "/")), tail, 1)), "-"), head, 1)) %>%
      separate(name_rev, into = c('date', 'cdc_id'), sep = 8) %>%
      #mutate(cdc_id = str_split(cdc_id, "_", 2))
      separate_wider_delim(cdc_id, delim = ".", names = c("cdc_id", "drop"))
    
    # run the function
    window_pi_file <- windowed_diversity_function(sample_results = df_i,
                                                  genome_size = 30000,
                                                  window_size = 1000,
                                                  window_step = 100)
    
    # save the output in a folder
    write_tsv(window_pi_file, file = paste(path, head(df_i$date, 1), head(df_i$cdc_id, 1), ".tsv", sep =""))
  }
  
  
  # stop cluster
  parallel::stopCluster(cl = my.cluster)
}

# for loop 3 function - mean pi per sample
parallel_mean_diversity_function <- function(data_dir_1,
                                             save_path,
                                             data_dir_2,
                                             bp_start,
                                             bp_end,
                                             final_file_path){
  
  # docstring information
  #' Calculate mean diversity in parallel
  #' 
  #' @description Estimate mean Pi and mean Shannon per window in parallel for
  #' specific genome regions.
  #' 
  #' @param data_dir_1 file path. Directory path for where the genetic data
  #' are saved.
  #' @param save_path file path. Directory path to save the output.
  #' @param data_dir_2 file path. Directory path for the processed genetic data
  #' for the region of interest.
  #' @param bp_start numeric. Base pair start position for the genome region of 
  #' interest.
  #' @param bp_end numeric. Base pair end position for the genome region of 
  #' interest.
  #' @param final_file_path. Direcotry path and file name for the saved output.
  #' 
  #' @details This function uses the results of `windowed_diversity_function` 
  #' to estimate mean Pi values and mean Shannon H values per sample for 
  #' specific genome regions. The function will process the individual variant
  #' data files in parallel and output one file with one row per sample.
  #' 
  #' @return Comma separated file (CSV).
  #' @examples
  #' # Genomewide Mean pi per sample
  #' parallel_mean_diversity_function(data_dir_1 = "data/windowed_pi/",
  #'                                  save_path = "data/genomewide_pi/",
  #'                                  data_dir_2 = "data/genomewide_pi/",
  #'                                  bp_start = 1,
  #'                                  bp_end = 30000,
  #'                                  final_file_path = 
  #'                                    paste(
  #'                                      "data/genomwide_pi_",
  #'                                      Sys.Date(),
  #'                                      ".csv",
  #'                                      sep = "")
  #' )
  #' 
  #' 
  # parallel
  library(foreach)
  library(doParallel)
  
  parallel::detectCores()
  n.cores <- parallel::detectCores()/2
  
  #create the cluster
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  
  #check cluster definition (optional)
  print(my.cluster)
  
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  #check if it is registered (optional)
  foreach::getDoParRegistered()
  
  #how many workers are available? (optional)
  foreach::getDoParWorkers()
  
  # files to iterate over
  data_dir <- data_dir_1
  list_of_files <- list.files(path = data_dir,
                              recursive = TRUE,
                              pattern = "\\.tsv$",
                              full.names = TRUE)
  
  # run process in parallel
  foreach(
    i = list_of_files
  ) %dopar% {
    
    library(dplyr)
    library(ggplot2)
    library(stringr)
    library(tidyr)
    library(purrr)
    library(readr)
    
    # path for files
    path <- save_path
    
    # file for the input of the function
    df_i <- i %>%
      set_names() %>%  
      map_df(read_tsv, .id = "file_name") %>% 
      mutate(name_rev = sapply(str_split((sapply((str_split(file_name, "/")), tail, 1)), "-"), head, 1)) %>%
      separate(name_rev, into = c('date', 'cdc_id'), sep = 8) %>%
      #mutate(cdc_id = str_split(cdc_id, "_", 2))
      separate_wider_delim(cdc_id, delim = ".", names = c("cdc_id", "drop"))
    
    # filter for the S1 RBD only -> note, might need to redo step 2 to better analyze this one
    df_i <- df_i %>%
      filter(start >= bp_start) %>%
      filter(start <= bp_end)
    
    
    # function for the iteration 
    # mean pi value, calculate outside of function so we can save the intermediate data
    g_diversity <- mean(df_i$avg_pi, na.rm = TRUE)
    shannon <- mean(df_i$avg_h, na.rm = TRUE)
    depth <- mean(df_i$avg_depth, na.rm = TRUE)
    
    # final dataset
    g_div_df <- as.data.frame(cbind(head(df_i$date, 1), 
                                    head(df_i$cdc_id, 1), 
                                    g_diversity, 
                                    shannon, 
                                    depth)
    )
    colnames(g_div_df) <- c("date", "cdc_id", "genomewide_pi",
                            "genomewide_h",
                            "depth")
    
    
    # save the output in a folder
    write_tsv(g_div_df, file = paste(path, head(df_i$date, 1), head(df_i$cdc_id, 1), ".tsv", sep =""))
  }
  
  
  # stop cluster
  parallel::stopCluster(cl = my.cluster)
  
  # Load the mean files and combine into one final flat file
  
  # load and bind the files together into one 
  data_dir <- data_dir_2
  
  # parallel
  library(foreach)
  library(doParallel)
  
  parallel::detectCores()
  n.cores <- parallel::detectCores()/2
  
  #create the cluster
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  
  #check cluster definition (optional)
  print(my.cluster)
  
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  #check if it is registered (optional)
  foreach::getDoParRegistered()
  
  #how many workers are available? (optional)
  foreach::getDoParWorkers()
  
  # files to iterate over
  list_of_files <- list.files(path = data_dir,
                              recursive = TRUE,
                              pattern = "\\.tsv$",
                              full.names = TRUE)
  # run process in parallel
  genomwide_pi <-
    foreach(
      i = list_of_files
    ) %dopar% {
      
      library(dplyr)
      library(ggplot2)
      library(stringr)
      library(tidyr)
      library(purrr)
      library(readr)
      
      # file for the input of the function
      df_i <- i %>%
        set_names() %>%  
        map_df(read_tsv, .id = "file_name") %>% 
        mutate(name_rev = sapply(str_split((sapply((str_split(file_name, "/")), tail, 1)), "-"), head, 1)) %>%
        separate(name_rev, into = c('date', 'cdc_id'), sep = 8) %>%
        #mutate(cdc_id = str_split(cdc_id, "_", 2))
        separate_wider_delim(cdc_id, delim = ".", names = c("cdc_id", "drop"))
    }
  
  
  # stop cluster
  parallel::stopCluster(cl = my.cluster)
  
  genomwide_pi <- do.call(rbind, genomwide_pi)
  
  write.csv(genomwide_pi, final_file_path, row.names = FALSE)
  
}


# -----------------------------------------------------------------------------

# MAIN ANALYSIS AND PLOT FUNCTIONS

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
  } else if(glmer_option == "glmmtmb"){
    # Effron R2 - compares the actual to the predicted values 
    # against the residuals
    dataframe$outcome_name <- dataframe[[outcome]]
    
    Actual <- dataframe$outcome_name
    Predicted = predict(model , type="response")
    Residuals = residuals(model )
    
    effron_r <- efronRSquared(residual = Residuals, 
                              predicted = Predicted, 
                              statistic = "EfronRSquared")
    effron_r
    
    r2 <- effron_r
    
    
    r2 <- as.data.frame(cbind(names(r2), r2))
    colnames(r2) <- c("variable", "Est")
    aic <- AIC(model)
    aic <- as.data.frame(cbind("AIC", aic))
    colnames(aic) <- c("variable", "Est")
    
    # model coefficients, se, and p values
    cc <- coef(summary(model))
    cc <- within(as.data.frame(cc$cond),
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
                                         `Pr(>|z|)`)
    # add number of observations
    n <- nobs(model)
    n <- as.data.frame(cbind("n", n))
    colnames(n) <- c("variable", "Est")
    
    tab_export$Est <- as.character(tab_export$Est)
    tab_export_final <- bind_rows(tab_export, n, r2, aic)
    
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
    dplyr::select(group, outcome, test, `F statistic`, `P value`) %>%
    na.omit()
  
  g2 <- as.data.frame(grangertest(d_x ~ d_y, order = lag))
  colnames(g2) <- c("ResDF", "DF", "F statistic", "P value")
  
  # keep first row f statistic and p value
  g2 <- g2 %>%
    mutate(group = group,
           outcome = outcome,
           test = "y predicts x") %>%
    dplyr::select(group, outcome, test, `F statistic`, `P value`) %>%
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
    labs(x = expression(paste("S1 NTD ",Pi[ww])),
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
      sec.axis = sec_axis(~ ./ 10000, name = expression(paste("S1 NTD ",Pi[ww]))))+
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

# boxplot function

# make it a function and create all the plots
spatial_boxplot_function <- function(dataframe,
                                     aggregation_name,
                                     outcome_var,
                                     aggregation_label,
                                     pi_name,
                                     h_name,
                                     conc_name,
                                     freyja_name){
  
  dataframe$aggregation <- dataframe[[aggregation_name]]
  dataframe$outcome <- dataframe[[outcome_var]]
  dataframe$pi <- dataframe[[pi_name]]
  dataframe$h <- dataframe[[h_name]]
  dataframe$conc <- dataframe[conc_name]
  dataframe$freyja_name <- dataframe[[freyja_name]]
  
  # correlation by region
  region_cor <- dataframe %>%
    group_by(aggregation) %>%
    summarize(
      pi_cor = cor(pi,
                   outcome,
                   method = c("spearman"),
                   use = "na.or.complete"),
      h_cor = cor(h,
                  outcome,
                  method = c("spearman"),
                  use = "na.or.complete"),
      freyja_cor = cor(freyja_name,
                       outcome,
                       method = c("spearman"),
                       use = "na.or.complete"),
      conc_cor = cor(conc,
                     outcome,
                     method = c("spearman"),
                     use = "na.or.complete")
    ) %>%
    ungroup()
  
  # pivot longer
  region_cor <- region_cor %>%
    pivot_longer(
      cols = c("pi_cor",
               "h_cor",
               "freyja_cor",
               "conc_cor")
    )
  
  # factor the labels
  region_cor$group <- factor(region_cor$name,
                             levels = c("pi_cor",
                                        "h_cor",
                                        "freyja_cor",
                                        "conc_cor"),
                             labels = c(
                               "Concentration",
                               "Freyja variant count",
                               expression("S1 NTD H"[ww]),
                               expression("S1 NTD "~ pi[ww])))
  
  # panel - top for cases, bottome for hospitalizations
  boxplot <- 
    ggplot(data = region_cor, aes(x = group, y = value, fill = name))+
    geom_boxplot()+
    ylim(0,1)+
    labs(title = paste(aggregation_label, " weekly aggregation"),
         x = "",
         y = "Spearman correlation",
         fill = "Group")+
    theme_dth_1+
    scale_fill_manual(
      values = c(freyja_cor ="#e76254",
                 h_cor = "#72bcd5",
                 pi_cor = "#376795",
                 conc_cor = "#ef8a47"),
      labels = c(
        "Concentration",
        "Freyja variant count",
        expression("S1 NTD H"[ww]),
        expression("S1 NTD "~ pi[ww]))
    )+
    theme(legend.position = "right",
          axis.text.x = element_blank())
  
  return(boxplot)
}

# ------------------------------------------------------------------------------

# county scatterplot function
# county cor plot function
county_cor_plot_function <- function(dataframe,
                                     county_name,
                                     var_name,
                                     outcome_name,
                                     cor_value,
                                     pal_color,
                                     pop_value,
                                     outcome_text,
                                     p_val,
                                     var_text){
  
  df <- dataframe %>%
    filter(county == county_name)
  df$var_name = df[[var_name]]
  df$outcome_name = df[[outcome_name]]
  
  plot <-
    ggplot(data = df,
           aes(x = var_name,
               y = outcome_name))+
    geom_point(
      color = pal_color)+
    theme_dth_1+
    theme(legend.background = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.ticks = element_line(size = 0.5),
          plot.subtitle = element_text(hjust = 0.5))+
    labs(title = paste("\u03c1 = ", cor_value, ", p value ", p_val),
         x = var_text,
         y = outcome_text)
  
  return(plot)
}

#-------------------------------------------------------------------------------


# function for comparing the spearman correlation coefficients across spatial
# scales
spatial_cor_compare_function <- function(
    dataframe,
    aggregation,
    outcome_var,
    pi_var,
    h_var,
    freyja_var,
    conc_var){
  
  dataframe$pi_var <- dataframe[[pi_var]]
  dataframe$h_var <- dataframe[[h_var]]
  dataframe$freyja_var <- dataframe[[freyja_var]]
  dataframe$conc_var <- dataframe[[conc_var]]
  dataframe$outcome_var <- dataframe[[outcome_var]]
  
  # correlation across all sewersheds instead of one by one
  cor <- cor.test(dataframe$pi_var,
                  dataframe$outcome_var,
                  method = "spearman")
  
  pval_1 <- cor$p.value
  est_1 <- cor$estimate
  n_1 <- sum(!is.na(dataframe$pi_var))
  z1 <- FisherZ(est_1)
  
  # comparison to concentration
  cor <- cor.test(dataframe$conc_var,
                  dataframe$outcome_var,
                  method = "spearman")
  
  pval_2 <- cor$p.value
  est_2 <- cor$estimate
  n_2 <- sum(!is.na(dataframe$conc_var))
  z2 <- FisherZ(est_2)
  
  fisher_z <- (z1 - z2)/sqrt((1/(n_1-3)) + (1/(n_2-3)))
  
  aggregation <- aggregation
  comparison <- paste0("Pi ww ~ ", outcome_var)
  
  row1 <- as.data.frame(
    cbind(aggregation,
          comparison,
          n_1,
          round_2(est_1),
          round_2(fisher_z)))
  
  colnames(row1) <- c("Aggregation",
                      "Comparison",
                      "n",
                      "Spearman correlation coefficient",
                      "Fisher's Z")
  # ------------------------------------------
  # repeat for h
  # ------------------------------------------
  # correlation across all sewersheds instead of one by one
  cor <- cor.test(dataframe$h_var,
                  dataframe$outcome_var,
                  method = "spearman")
  
  pval_1 <- cor$p.value
  est_1 <- cor$estimate
  n_1 <- sum(!is.na(dataframe$h_var))
  z1 <- FisherZ(est_1)
  
  # comparison to concentration
  cor <- cor.test(dataframe$conc_var,
                  dataframe$outcome_var,
                  method = "spearman")
  
  pval_2 <- cor$p.value
  est_2 <- cor$estimate
  n_2 <- sum(!is.na(dataframe$conc_var))
  z2 <- FisherZ(est_2)
  
  fisher_z <- (z1 - z2)/sqrt((1/(n_1-3)) + (1/(n_2-3)))
  
  aggregation <- aggregation
  comparison <- paste0("H ww ~ ", outcome_var)
  
  row2 <- as.data.frame(
    cbind(aggregation,
          comparison,
          n_1,
          round_2(est_1),
          round_2(fisher_z)))
  
  colnames(row2) <- c("Aggregation",
                      "Comparison",
                      "n",
                      "Spearman correlation coefficient",
                      "Fisher's Z")
  # ------------------------------------------
  # repeat for freyja
  # ------------------------------------------
  # correlation across all sewersheds instead of one by one
  cor <- cor.test(dataframe$freyja_var,
                  dataframe$outcome_var,
                  method = "spearman")
  
  pval_1 <- cor$p.value
  est_1 <- cor$estimate
  n_1 <- sum(!is.na(dataframe$freyja_var))
  z1 <- FisherZ(est_1)
  
  # comparison to concentration
  cor <- cor.test(dataframe$conc_var,
                  dataframe$outcome_var,
                  method = "spearman")
  
  pval_2 <- cor$p.value
  est_2 <- cor$estimate
  n_2 <- sum(!is.na(dataframe$conc_var))
  z2 <- FisherZ(est_2)
  
  fisher_z <- (z1 - z2)/sqrt((1/(n_1-3)) + (1/(n_2-3)))
  
  aggregation <- aggregation
  comparison <- paste0("Variant count ~ ", outcome_var)
  
  row3 <- as.data.frame(
    cbind(aggregation,
          comparison,
          n_1,
          round_2(est_1),
          round_2(fisher_z)))
  
  colnames(row3) <- c("Aggregation",
                      "Comparison",
                      "n",
                      "Spearman correlation coefficient",
                      "Fisher's Z")
  # ------------------------------------------
  # repeat for conc
  # ------------------------------------------
  # comparison to concentration
  cor <- cor.test(dataframe$conc_var,
                  dataframe$outcome_var,
                  method = "spearman")
  
  pval_2 <- cor$p.value
  est_2 <- cor$estimate
  n_2 <- sum(!is.na(dataframe$conc_var))
  
  
  aggregation <- aggregation
  comparison <- paste0("Concentration ~ ", outcome_var)
  
  row4 <- as.data.frame(
    cbind(aggregation,
          comparison,
          n_2,
          round_2(est_2),
          ""))
  
  colnames(row4) <- c("Aggregation",
                      "Comparison",
                      "n",
                      "Spearman correlation coefficient",
                      "Fisher's Z")
  
  ###
  # combine dataframe
  ###
  df <- bind_rows(row1, row2, row3, row4)
  
  return(df)
}