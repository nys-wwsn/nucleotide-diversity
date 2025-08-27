###############
###############
## SEQUENCING OF SARS-COV-2 GENOME DIVERSITY PROJECT
###############
###############

# Author - Dustin T Hill

# Created - May 19, 2025
# Last updated - August 27, 2025, 2025

# Description: This script loads each Freyja summary file (TSV format) and 
# processes to calculate nucleotide diveristy (Pi) and Shannon's H. Each
# file from Freyja represents one wastewater sample and the file name includes
# the sample date and location ID. 

# Docstring R packages was utilized to provide help results for R functions 
# created in this script. After loading docstring with library(docstring),
# users can query the help window wit ?function_name.

# Packages
# install.packages("dplyr")
# install.packages("readr")
# install.packages("tidyr")
# install.packages("foreach")
# install.packages("doParallel")
# install.packages("stringr")
# install.packages("purrr")
# install.packages("docstring)

library(docstring)

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
    select(POS, REF_DP) %>% 
    distinct(POS, .keep_all = TRUE) %>% 
    rename(k = REF_DP) %>%
    mutate(freq = k)
  alternate_counts <- input %>% 
    select(POS, ALT_DP) %>% 
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
        select(POS, REF_DP) %>% 
        distinct(POS, .keep_all = TRUE) %>% 
        rename(k = REF_DP) %>%
        mutate(freq = k)
      alternate_counts <- input %>% 
        select(POS, ALT_DP) %>% 
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

# GENOMEWIDE PI
# calculate genomewide pi for each TSV file

# data from July 2024 to April 2025
parallel_diversity_per_base_function(data_dir = "variants/",
                              save_path = "data/pi/",
                              pass_option = "pass only")

# Windowed pi for the genome
parallel_windowed_diversity_function(data_dir = "data/pi/",
                              save_path = "data/windowed_pi/")

# Genomewide Mean pi per sample
parallel_mean_diversity_function(data_dir_1 = "data/windowed_pi/",
                          save_path = "data/genomewide_pi/",
                          data_dir_2 = "data/genomewide_pi/",
                          bp_start = 1,
                          bp_end = 30000,
                          final_file_path = 
                            paste(
                              "data/genomwide_pi_",
                          Sys.Date(),
                          ".csv",
                          sep = "")
                          )

# From here, we calculate mean pi values for specific regions using the windowed
# pi files, not going back to the beginning

# Region specific steps

# orf 5 and 6
parallel_mean_diversity_function(data_dir_1 = "data/windowed_pi/",
                          save_path = "data/orf_region_nsp5_6/",
                          data_dir_2 = "data/orf_region_nsp5_6/",
                          bp_start = 10063,
                          bp_end = 11842,
                          final_file_path = 
                            paste(
                              "data/orf_nsp5_6_pi_",
                              Sys.Date(),
                              ".csv",
                              sep = "")
)

# spike protein
parallel_mean_diversity_function(data_dir_1 = "data/windowed_pi/",
                          save_path = "data/complete_spike_pi/",
                          data_dir_2 = "data/complete_spike_pi/",
                          bp_start = 21563,
                          bp_end = 25384,
                          final_file_path = 
                            paste(
                              "data/spike_pi_",
                              Sys.Date(),
                              ".csv",
                              sep = "")
)

# cov methyl 2
parallel_mean_diversity_function(data_dir_1 = "data/windowed_pi/",
                          save_path = "data/cov_mt_2/",
                          data_dir_2 = "data/cov_mt_2/",
                          bp_start = 20661,
                          bp_end = 21549,
                          final_file_path = 
                            paste(
                              "data/cov_mt_2_",
                              Sys.Date(),
                              ".csv",
                              sep = "")
)

# s1 ntd
parallel_mean_diversity_function(data_dir_1 = "data/windowed_pi/",
                          save_path = "data/s1_ntd_pi/",
                          data_dir_2 = "data/s1_ntd_pi/",
                          bp_start = 21598,
                          bp_end = 22474,
                          final_file_path = 
                            paste(
                              "data/s1_ntd_pi_",
                              Sys.Date(),
                              ".csv",
                              sep = "")
)


# s1 RBD
parallel_mean_diversity_function(data_dir_1 = "data/windowed_pi/",
                          save_path = "data/s1_rbd_pi/",
                          data_dir_2 = "data/s1_rbd_pi/",
                          bp_start = 22516,
                          bp_end = 23185,
                          final_file_path = 
                            paste(
                              "data/s1_rbd_pi_",
                              Sys.Date(),
                              ".csv",
                              sep = "")
)

