###############
###############
## SEQUENCING OF SARS-COV-2 GENOME DIVERSITY PROJECT
###############
###############

# Author - Dustin T Hill

# Created - May 19, 2025
# Last updated - August 27, 2025

# Description: This script loads each Freyja summary file (TSV format) and 
# processes to calculate nucleotide diveristy (Pi) and Shannon's H. Each
# file from Freyja represents one wastewater sample and the file name includes
# the sample date and location ID. 

# Docstring R package was utilized to provide help results for R functions 
# created in this script. After loading docstring with library(docstring),
# users can query the help window witH ?function_name.

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

# LOAD FUNCTIONS
source("seq diversity - functions.R")

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
                          save_path = "data/diversity_output/genomewide_pi/",
                          data_dir_2 = "data/diversity_output/genomewide_pi/",
                          bp_start = 1,
                          bp_end = 30000,
                          final_file_path = 
                            paste(
                              "data/diversity_output/genomwide_pi_",
                          Sys.Date(),
                          ".csv",
                          sep = "")
                          )

# From here, we calculate mean pi values for specific regions using the windowed
# pi files, not going back to the beginning

# Region specific steps

# orf 5 and 6
parallel_mean_diversity_function(data_dir_1 = "data/windowed_pi/",
                          save_path = "data/diversity_output/orf_region_nsp5_6/",
                          data_dir_2 = "data/diversity_output/orf_region_nsp5_6/",
                          bp_start = 10063,
                          bp_end = 11842,
                          final_file_path = 
                            paste(
                              "data/diversity_output/orf_nsp5_6_pi_",
                              Sys.Date(),
                              ".csv",
                              sep = "")
)

# spike protein
parallel_mean_diversity_function(data_dir_1 = "data/diversity_output/windowed_pi/",
                          save_path = "data/diversity_output/complete_spike_pi/",
                          data_dir_2 = "data/diversity_output/complete_spike_pi/",
                          bp_start = 21563,
                          bp_end = 25384,
                          final_file_path = 
                            paste(
                              "data/diversity_output/spike_pi_",
                              Sys.Date(),
                              ".csv",
                              sep = "")
)

# cov methyl 2
parallel_mean_diversity_function(data_dir_1 = "data/windowed_pi/",
                          save_path = "data/diversity_output/cov_mt_2/",
                          data_dir_2 = "data/diversity_output/cov_mt_2/",
                          bp_start = 20661,
                          bp_end = 21549,
                          final_file_path = 
                            paste(
                              "data/diversity_output/cov_mt_2_",
                              Sys.Date(),
                              ".csv",
                              sep = "")
)

# s1 ntd
parallel_mean_diversity_function(data_dir_1 = "data/windowed_pi/",
                          save_path = "data/diversity_output/s1_ntd_pi/",
                          data_dir_2 = "data/diversity_output/s1_ntd_pi/",
                          bp_start = 21598,
                          bp_end = 22474,
                          final_file_path = 
                            paste(
                              "data/diversity_output/s1_ntd_pi_",
                              Sys.Date(),
                              ".csv",
                              sep = "")
)


# s1 RBD
parallel_mean_diversity_function(data_dir_1 = "data/windowed_pi/",
                          save_path = "data/diversity_output/s1_rbd_pi/",
                          data_dir_2 = "data/diversity_output/s1_rbd_pi/",
                          bp_start = 22516,
                          bp_end = 23185,
                          final_file_path = 
                            paste(
                              "data/diversity_output/s1_rbd_pi_",
                              Sys.Date(),
                              ".csv",
                              sep = "")
)

