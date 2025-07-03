###########
###########
#### GENOMEWIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-06-30

# DATA PREP SCRIPT
# Prepare genome sequence data, case data, hospitalization data, quantification
# data for analysis. This script pulls from individual files and combines
# into a final output file for use in the "seq diversity- analysis.R" file

# This script also creates the descriptive statistics tables and figures

# --------------------------------------
# PACKAGES
# --------------------------------------

library(zoo)
library(stringr)
library(lubridate)
library(tidyr)
library(dplyr)

# --------------------------------------
# FUNCTIONS
# --------------------------------------

# --------------------------------------
# DATA PREPARATION FOR GENOME
# --------------------------------------

# genomewide
genome <- read.csv("data/genomewide_pi.csv", stringsAsFactors = FALSE)

# edit date filed
genome$date <- ymd(genome$date)

# weekly mean
genome_weekly <- genome %>%
  group_by(cdc_id, week = floor_date(date, unit = "weeks", week_start = 7)) %>%
  summarize(genomewide_pi = mean(genomewide_pi, na.rm = TRUE),
            genomewide_h = mean(genomewide_h, na.rm = TRUE),
            depth = mean(depth, na.rm = TRUE),
            samples = n()
  ) %>%
  ungroup()

# spike

# spike data
spike <- read.csv("data/spike_pi.csv")

# edit date filed
spike$date <- ymd(spike$date)

# weekly mean
spike_weekly <- spike %>%
  group_by(cdc_id, week = floor_date(date, unit = "weeks", week_start = 7)) %>%
  summarize(spike_pi = mean(genomewide_pi, na.rm = TRUE),
            spike_h = mean(genomewide_h, na.rm = TRUE),
            spike_depth = mean(depth, na.rm = TRUE)
  ) %>%
  ungroup()

diversity_df <- left_join(genome_weekly, spike_weekly, by = c("cdc_id", "week"))

# s1 ntd
ntd <- read.csv("data/s1_ntd_pi.csv")

# edit date filed
ntd$date <- ymd(ntd$date)

# weekly mean
ntd_weekly <- ntd %>%
  group_by(cdc_id, week = floor_date(date, unit = "weeks", week_start = 7)) %>%
  summarize(ntd_pi = mean(genomewide_pi, na.rm = TRUE),
            ntd_h = mean(genomewide_h, na.rm = TRUE),
            ntd_depth = mean(depth, na.rm = TRUE)
  ) %>%
  ungroup()

# merge to main file
diversity_df <- left_join(diversity_df, ntd_weekly, by = c("cdc_id", "week"))

# s1 rbd
rbd <- read.csv("data/s1_rbd_pi.csv")

# edit date filed
rbd$date <- ymd(rbd$date)

# weekly mean
rbd_weekly <- rbd %>%
  group_by(cdc_id, week = floor_date(date, unit = "weeks", week_start = 7)) %>%
  summarize(rbd_pi = mean(genomewide_pi, na.rm = TRUE),
            rbd_h = mean(genomewide_h, na.rm = TRUE),
            rbd_depth = mean(depth, na.rm = TRUE)
  ) %>%
  ungroup()

# merge to main file
diversity_df <- left_join(diversity_df, rbd_weekly, by = c("cdc_id", "week"))

# orf
orf <- read.csv("data/orf_nsp5_6_pi.csv")

# edit date filed
orf$date <- ymd(orf$date)

# weekly mean
orf_weekly <- orf %>%
  group_by(cdc_id, week = floor_date(date, unit = "weeks", week_start = 7)) %>%
  summarize(orf_pi = mean(genomewide_pi, na.rm = TRUE),
            orf_h = mean(genomewide_h, na.rm = TRUE),
            orf_depth = mean(depth, na.rm = TRUE)
  ) %>%
  ungroup()

# merge to main file
diversity_df <- left_join(diversity_df, orf_weekly, by = c("cdc_id", "week"))

# cov mt 2
cov_mt_2 <- read.csv("data/cov_mt_2_pi.csv")

# edit date filed
cov_mt_2$date <- ymd(cov_mt_2$date)

# weekly mean
cov_mt_2_weekly <- cov_mt_2 %>%
  group_by(cdc_id, week = floor_date(date, unit = "weeks", week_start = 7)) %>%
  summarize(cov_mt_2_pi = mean(genomewide_pi, na.rm = TRUE),
            cov_mt_2_h = mean(genomewide_h, na.rm = TRUE),
            cov_mt_2_depth = mean(depth, na.rm = TRUE)
  ) %>%
  ungroup()

# add to main file
diversity_df <- left_join(diversity_df, cov_mt_2_weekly, 
                          by = c("cdc_id", "week")
                          )

diversity_df <- diversity_df %>%
  rename(facility_id = cdc_id)

# -----------------------------------------------------------------------------

# --------------------------------------
# DATA PREPARATION FOR QUANT DATA
# --------------------------------------

# data load - meta
meta <- read.csv("data/nys-wws-sewersheds.csv", stringsAsFactors = FALSE)

# data load - quant data
quant.data <- readRDS("data/sars2_concentration.rds")

# Calculating copies of SARS2 RNA #
quant.data$copies <- 3.5
quant.data$copies <- ifelse(
  !is.na(quant.data$sars2_copies_ml), 
  as.numeric(as.character(quant.data$sars2_copies_ml)), 
  quant.data$copies)
quant.data$copies <- ifelse(quant.data$sars_pos==0, 1, quant.data$copies)
quant.data$copies <- ifelse(
  is.na(quant.data$copies), 
  as.numeric(as.character(quant.data$sars2_copies_ml)), 
  quant.data$copies)
quant.data$copies <- ifelse(is.na(quant.data$copies), 1, quant.data$copies)
quant.data$copies[quant.data$sars2_copies_ml == 0] <- 1
quant.data$copies <- ifelse(quant.data$copies == 0, 1, quant.data$copies)
quant.data$copies[quant.data$copies < 1] <- 1

# remove values > 100000
quant.data  <- quant.data %>%
  filter(copies < 100000)

# weekly mean of the quant data for now
conc_mean <- quant.data %>%
  group_by(week = floor_date(sample_collect_date, 
                             unit = "weeks", 
                             week_start = 7), 
           facility_id
           ) %>%
  summarize(mean_sars2_conc = mean(copies, na.rm = TRUE),
            conc_samples = n()
  ) %>%
  ungroup()

# -----------------------------------------------------------------------------

# --------------------------------------
# DATA PREPARATION FOR LINEAGE DATA
# --------------------------------------

# load data
load(file = "data/lineage name files.R")

# count lineages per week
var_unique_sewer <- var.data %>%
  group_by(facility_id,
    week = lubridate::floor_date(date, 
                                 unit = 'weeks', 
                                 week_start = 7),
    region,
    county
  ) %>%
  mutate(n_variants_no_thresh = length(unique(lineage))
  ) %>%
  ungroup() %>%
  dplyr::filter(variant_pct > 0.05)%>%
  group_by(facility_id,
    week = lubridate::floor_date(date, 
                                 unit = 'weeks', 
                                 week_start = 7),
    region,
    county
  ) %>%
  mutate(n_variants_5 = length(unique(lineage))
  ) %>%
  ungroup() %>%
  filter(week <= as.Date("2025-04-21")) %>%
  select(week, 
         facility_id, 
         county, 
         region, 
         n_variants_5, 
         n_variants_no_thresh) %>%
  distinct()

# number of samples  by region each week
sample_count <- var.data %>%
  group_by(region, week = floor_date(date, 
                                     unit = 'weeks', 
                                     week_start = 7) ) %>%
  summarize(samples_by_region = length(unique(facility_id))
  ) %>%
  ungroup()%>%
  filter(week <= as.Date("2025-04-21"))

saveRDS(sample_count, "data/sample_count.rds")

# count lineages per week by county
var_unique_county <- var.data %>%
  group_by(
           week = lubridate::floor_date(date, 
                                        unit = 'weeks', 
                                        week_start = 7),
           region,
           county
  ) %>%
  mutate(n_variants_no_thresh = length(unique(lineage))
  ) %>%
  ungroup() %>%
  dplyr::filter(variant_pct > 0.05)%>%
  group_by(
    week = lubridate::floor_date(date, 
                                 unit = 'weeks', 
                                 week_start = 7),
    region,
    county
  ) %>%
  mutate(n_variants_5 = length(unique(lineage))
  ) %>%
  ungroup() %>%
  filter(week <= as.Date("2025-04-21")) %>%
  select(week, region, county, n_variants_5, n_variants_no_thresh) %>%
  distinct()

# variant county by region
var_unique_region <- var.data %>%
  group_by(
    week = lubridate::floor_date(date, 
                                 unit = 'weeks', 
                                 week_start = 7),
    region
  ) %>%
  mutate(n_variants_no_thresh = length(unique(lineage))
  ) %>%
  ungroup() %>%
  dplyr::filter(variant_pct > 0.05)%>%
  group_by(
    week = lubridate::floor_date(date, 
                                 unit = 'weeks', 
                                 week_start = 7),
    region
  ) %>%
  mutate(n_variants_5 = length(unique(lineage))
  ) %>%
  ungroup() %>%
  filter(week <= as.Date("2025-04-21")) %>%
  select(week, region,  n_variants_5, n_variants_no_thresh) %>%
  distinct()

# variant county statewide
nyc <- c("New York", "Bronx", "Queens", "Kings", "Richmond")

var_unique_state <- var.data %>%
  mutate(group = ifelse(county %in% nyc, "NYC", "State")
         ) %>%
  group_by(
    week = lubridate::floor_date(date, 
                                 unit = 'weeks', 
                                 week_start = 7),
    group
  ) %>%
  mutate(n_variants_no_thresh = length(unique(lineage))
  ) %>%
  ungroup() %>%
  dplyr::filter(variant_pct > 0.05)%>%
  group_by(
    week = lubridate::floor_date(date, 
                                 unit = 'weeks', 
                                 week_start = 7),
    group
  ) %>%
  mutate(n_variants_5 = length(unique(lineage))
  ) %>%
  ungroup() %>%
  filter(week <= as.Date("2025-04-21")) %>%
  select(week, group, n_variants_5, n_variants_no_thresh) %>%
  distinct()

# -----------------------------------------------------------------------------

# --------------------------------------
# DATA PREPARATION FOR HOSPITALIZATIONS
# --------------------------------------

# load data
hosp_data <- readRDS("data/covid_hosp.rds")
county_region <- read.csv("data/County_Region.csv")
county_region$county[county_region$county == "St. Lawrence"] <- "St Lawrence"

# fix date
hosp_data$date <- strptime(hosp_data$as_of_date, "%m/%d/%Y")
hosp_data$date <- ymd(hosp_data$date)
hosp_data$county <- str_to_title(hosp_data$facility_county)
hosp_data$county <- ifelse(hosp_data$county == "St. Lawrence", 
                           "St Lawrence", 
                           hosp_data$county)

# sum by county
hosp_county <- hosp_data %>%
  group_by(week = floor_date(date, 
                             unit = "weeks", 
                             week_start = 7), county) %>%
  summarise(
    sum_Total.New.Admissions.Reported = sum(total_new_admissions_reported, 
                                            na.rm = TRUE)
  ) %>%
  ungroup()

# Add per 100k
per100 <- county_region %>% dplyr::select(county, county_pop)
hosp_county <- left_join(hosp_county, per100, by = c("county"))
hosp_county$sum_Total.New.Admissions.Reported_100 <- 
  hosp_county$sum_Total.New.Admissions.Reported * 
  (100000/hosp_county$county_pop)

# drop columns we do not need and rename columns
hosp_county_weekly <- hosp_county %>%
  dplyr::select(week, 
                county, 
                county_pop, 
                sum_Total.New.Admissions.Reported, 
                sum_Total.New.Admissions.Reported_100
                ) %>%
  rename(hospitalizations = sum_Total.New.Admissions.Reported,
         hosp_incidence = sum_Total.New.Admissions.Reported_100)

# filter for the dates we have variant data for plus three extra weeks for lags
hosp_county_weekly <- hosp_county_weekly %>%
  filter(week >= min(diversity_df$week) & week <= max(diversity_df$week)+weeks(3))


# -----------------------------------------------------------------------------

# --------------------------------------
# DATA PREPARATION FOR CASE DATA
# --------------------------------------

# load data
all.cases <- readRDS("data/covid_cases.rds")

# data prep (date variable change and make case data numeric)
all.cases$date <- mdy(all.cases$test_date)
all.cases$cases <- as.numeric(gsub(",", "", all.cases$total_new_positives))
all.cases$county <- all.cases$geography_description
all.cases$county[all.cases$county == "St. Lawrence"] <- "St Lawrence"

# sum per week
all.cases_weekly <- all.cases %>%
  arrange(date) %>%
  group_by(county, week = floor_date(date, unit = "weeks", week_start = 7)) %>%
  summarize(cases = sum(cases, na.rm = TRUE)
  ) %>%
  ungroup()

# remove region names from county list
names_remove <- c("Capital Region ", 
                  "Central New York", 
                  "Long Island", 
                  "Mid-Hudson", 
                  "Mohawk Valley", 
                  "New York City", 
                  "North Country",
                  "Southern Tier", 
                  "STATEWIDE", 
                  "Western New York", 
                  "Finger Lakes")
cases_county_weekly <- all.cases_weekly %>% filter(!county %in% names_remove)

# add case incidence
cases_county_weekly <- left_join(cases_county_weekly, per100, by = c("county"))
cases_county_weekly$case_incidence <- cases_county_weekly$cases * 
  (100000/cases_county_weekly$county_pop)

# filter to var data start date and end date plus three weeks
cases_county_weekly <- cases_county_weekly %>%
  filter(week >= min(diversity_df$week) & week <= max(diversity_df$week)+weeks(3))

# st lawrence county fix


# -----------------------------------------------------------------------------

# --------------------------------------
# SEWERSHED LEVEL COMBINED DATA FILE
# --------------------------------------

dat_sewershed <- full_join(diversity_df, 
                           conc_mean, 
                           by = c("facility_id", "week")) %>%
  filter(week >= (as.Date("2023-01-01"))) %>%
  filter(week <= "2025-04-21")

# if missing conc, drop the data
dat_sewershed <- dat_sewershed %>%
  filter(!is.na(mean_sars2_conc)) %>%
  distinct()

# check for duplicates
dups <- dat_sewershed %>%
  group_by(facility_id, week) %>%
  filter(n() > 1) %>%
  distinct()

# fill out the data with missing weeks
week_max <- max(diversity_df$week, na.rm = TRUE)
week_min <- min(dat_sewershed$week, na.rm = TRUE)
dat_sewershed <- dat_sewershed %>%
  group_by(facility_id) %>%
  complete(week = seq.Date(from = ymd(week_min), 
                           to = ymd(week_max), 
                           by = "weeks")
  ) %>%
  ungroup() %>%
  filter(week <= week_max)

# linear interpolate missing pi/h values
dat_sewershed <- dat_sewershed %>%
  group_by(facility_id) %>%
  dplyr::arrange(week)%>%
  dplyr::mutate(# missing pi values
                genomewide_pi_approx= na.approx(genomewide_pi, 
                                                maxgap = 3, 
                                                na.rm = FALSE),
                spike_pi_approx= na.approx(spike_pi, 
                                           maxgap = 3, 
                                           na.rm = FALSE),
                orf_pi_approx= na.approx(orf_pi, 
                                         maxgap = 3, na.rm = FALSE),
                ntd_pi_approx= na.approx(ntd_pi, maxgap = 3, 
                                         na.rm = FALSE),
                rbd_pi_approx= na.approx(rbd_pi, 
                                         maxgap = 3, 
                                         na.rm = FALSE),
                cov_mt_2_pi_approx = na.approx(cov_mt_2_pi, 
                                               maxgap = 3, 
                                               na.rm = FALSE),
                # missing concentration
                mean_sars2_conc_approx = na.approx(mean_sars2_conc, 
                                                   maxgap = 3, 
                                                   na.rm = FALSE),
                # missing h values
                genomewide_h_approx= na.approx(genomewide_h, 
                                                maxgap = 3, 
                                                na.rm = FALSE),
                spike_h_approx= na.approx(spike_h, 
                                           maxgap = 3, 
                                           na.rm = FALSE),
                orf_h_approx= na.approx(orf_h, 
                                         maxgap = 3, na.rm = FALSE),
                ntd_h_approx= na.approx(ntd_h, maxgap = 3, 
                                         na.rm = FALSE),
                rbd_h_approx= na.approx(rbd_h, 
                                         maxgap = 3, 
                                         na.rm = FALSE),
                cov_mt_2_h_approx = na.approx(cov_mt_2_h, 
                                               maxgap = 3, 
                                               na.rm = FALSE),
                
                # missing depth
                depth_approx = na.approx(depth,
                                         maxgap = 3,
                                         na.rm = FALSE),
                depth_ntd_approx = na.approx(ntd_depth,
                                             maxgap = 3,
                                             na.rm = FALSE)
                
  )%>%
  dplyr::ungroup()

# review na's
st <- dat_sewershed %>%
  filter(is.na(genomewide_pi_approx))

# 3 week rolling average

# rolling average per site
dat_sewershed <- dat_sewershed %>%
  group_by(facility_id) %>%
  dplyr::arrange(week)%>%
  dplyr::mutate(
    genomewide_pi_ma3 = rollmean(genomewide_pi_approx, 
                                 3, 
                                 align = "right", 
                                 na.pad = TRUE, 
                                 na.rm = TRUE),
    spike_pi_ma3 = rollmean(spike_pi_approx, 
                            3, 
                            align = "right", 
                            na.pad = TRUE, 
                            na.rm = TRUE),
    orf_pi_ma3 = rollmean(orf_pi_approx, 3, 
                          align = "right", 
                          na.pad = TRUE, 
                          na.rm = TRUE),
    ntd_pi_ma3 = rollmean(ntd_pi_approx, 
                          3, 
                          align = "right", 
                          na.pad = TRUE, 
                          na.rm = TRUE),
    rbd_pi_ma3 = rollmean(rbd_pi_approx, 
                          3, 
                          align = "right", 
                          na.pad = TRUE, 
                          na.rm = TRUE),
    cov_mt_2_pi_ma3 = rollmean(cov_mt_2_pi_approx, 
                               3, 
                               align = "right", 
                               na.pad = TRUE, 
                               na.rm = TRUE),
    mean_sars2_conc_ma3 = rollmean(mean_sars2_conc_approx, 
                                   3, 
                                   align = "right", 
                                   na.pad = TRUE, 
                                   na.rm = TRUE),
    # h rolling average
    genomewide_h_ma3 = rollmean(genomewide_h_approx, 
                                 3, 
                                 align = "right", 
                                 na.pad = TRUE, 
                                 na.rm = TRUE),
    spike_h_ma3 = rollmean(spike_h_approx, 
                            3, 
                            align = "right", 
                            na.pad = TRUE, 
                            na.rm = TRUE),
    orf_h_ma3 = rollmean(orf_h_approx, 3, 
                          align = "right", 
                          na.pad = TRUE, 
                          na.rm = TRUE),
    ntd_h_ma3 = rollmean(ntd_h_approx, 
                          3, 
                          align = "right", 
                          na.pad = TRUE, 
                          na.rm = TRUE),
    rbd_h_ma3 = rollmean(rbd_h_approx, 
                          3, 
                          align = "right", 
                          na.pad = TRUE, 
                          na.rm = TRUE),
    cov_mt_2_h_ma3 = rollmean(cov_mt_2_h_approx, 
                               3, 
                               align = "right", 
                               na.pad = TRUE, 
                               na.rm = TRUE),
    
    # rolling avg depth
    depth_ma3 = rollmean(depth_approx,
                         3,
                         align = "right",
                         na.pad = TRUE,
                         na.rm = TRUE),
    depth_ntd_ma3 = rollmean(depth_ntd_approx,
                             3,
                             align = "right",
                             na.pad = TRUE,
                             na.rm = TRUE)
    
    
  )%>%
  dplyr::ungroup()

# add variant counts
dat_sewershed <- left_join(dat_sewershed, var_unique_sewer, 
                           by = c("week", "facility_id"))

# remove na sewersheds
dat_sewershed <- dat_sewershed %>%
  filter(!is.na(facility_id))

# check for dups
dups <- dat_sewershed %>%
  group_by(facility_id, week) %>%
  filter(n() > 1)

# --------------------------------------
# COUNTY LEVEL COMBINED DATA FILE
# --------------------------------------

meta_c <- meta %>%
  select(facility_id, county, population_served)

# now average by county
dat_county <- dat_sewershed %>%
  left_join(meta_c, by = c("facility_id", "county")) %>%
  arrange(week) %>%
  group_by(week, county)%>%
  summarise(
    # pi values interpolate
    genomewide_pi_county_3w = weighted.mean(x = genomewide_pi_ma3, 
                                            w = population_served, 
                                            na.rm = TRUE),
    spike_pi_county_3w = weighted.mean(x = spike_pi_ma3, 
                                       w = population_served, 
                                       na.rm = TRUE),
    orf_pi_county_3w = weighted.mean(x = orf_pi_ma3, 
                                     w = population_served, 
                                     na.rm = TRUE),
    ntd_pi_county_3w = weighted.mean(x = ntd_pi_ma3, 
                                     w = population_served, 
                                     na.rm = TRUE),
    rbd_pi_county_3w = weighted.mean(x = rbd_pi_ma3, 
                                     w = population_served, 
                                     na.rm = TRUE),
    cov_mt_2_pi_county_3w = weighted.mean(x = cov_mt_2_pi_ma3, 
                                          w = population_served, 
                                          na.rm = TRUE),
    # concentration interpolate
    mean_sars2_conc_county_3w = weighted.mean(x = mean_sars2_conc_ma3, 
                                              w = population_served, 
                                              na.rm = TRUE),
    # h values interpolate
    genomewide_h_county_3w = weighted.mean(x = genomewide_h_ma3, 
                                            w = population_served, 
                                            na.rm = TRUE),
    spike_h_county_3w = weighted.mean(x = spike_h_ma3, 
                                       w = population_served, 
                                       na.rm = TRUE),
    orf_h_county_3w = weighted.mean(x = orf_h_ma3, 
                                     w = population_served, 
                                     na.rm = TRUE),
    ntd_h_county_3w = weighted.mean(x = ntd_h_ma3, 
                                     w = population_served, 
                                     na.rm = TRUE),
    rbd_h_county_3w = weighted.mean(x = rbd_h_ma3, 
                                     w = population_served, 
                                     na.rm = TRUE),
    cov_mt_2_h_county_3w = weighted.mean(x = cov_mt_2_h_ma3, 
                                          w = population_served, 
                                          na.rm = TRUE),
    
    # depth
    depth_county_3w = weighted.mean(x = depth_ma3,
                                    w = population_served,
                                    na.rm = TRUE),
    
    depth_ntd_county_3w = weighted.mean(x = depth_ntd_ma3,
                                        w = population_served,
                                        na.rm = TRUE)
  )%>%
  ungroup() %>%
  filter(week <= week_max)

# add case / hosp data - full join
dat_county <- full_join(dat_county, 
                        hosp_county_weekly, 
                        by = c("week", "county"))
dat_county <- full_join(dat_county, 
                        cases_county_weekly, 
                        by = c("week", "county"))

# fill na hosp incidence and hospital admissions with 0
dat_county$hosp_incidence <- ifelse(is.na(dat_county$hosp_incidence), 
                                    0, 
                                    dat_county$hosp_incidence)
dat_county$hospitalizations <- ifelse(is.na(dat_county$hospitalizations), 
                                      0, 
                                      dat_county$hospitalizations)

# remove duplicate pop column
dat_county <- dat_county %>%
  select(-county_pop.x) %>%
  rename(county_pop = county_pop.y)

# add variant counts
dat_county <- left_join(dat_county, var_unique_county, by = c("week", "county"))

# check for dups
dups <- dat_county %>%
  group_by(county, week) %>%
  filter(n() > 1)

# --------------------------------------
# REGIONAL LEVEL COMBINED DATA FILE
# --------------------------------------

meta_r <- meta %>%
  select(facility_id, county, region, population_served) 

# now average by region
dat_region <- dat_sewershed %>%
  left_join(meta_r, by = c("facility_id", "region", "county")) %>%
  arrange(week) %>%
  group_by(week, region)%>%
  summarise(
    # regional pi values
    genomewide_pi_region_3w = weighted.mean(x = genomewide_pi_ma3, 
                                            w = population_served, 
                                            na.rm = TRUE),
    spike_pi_region_3w = weighted.mean(x = spike_pi_ma3, 
                                       w = population_served, 
                                       na.rm = TRUE),
    orf_pi_region_3w = weighted.mean(x = orf_pi_ma3, 
                                     w = population_served, 
                                     na.rm = TRUE),
    ntd_pi_region_3w = weighted.mean(x = ntd_pi_ma3, 
                                     w = population_served, 
                                     na.rm = TRUE),
    rbd_pi_region_3w = weighted.mean(x = rbd_pi_ma3, 
                                     w = population_served, 
                                     na.rm = TRUE),
    cov_mt_2_pi_region_3w = weighted.mean(x = cov_mt_2_pi_ma3, 
                                          w = population_served, 
                                          na.rm = TRUE),
    # regional conc values
    mean_sars2_conc_region_3w = weighted.mean(x = mean_sars2_conc_ma3, 
                                              w = population_served, 
                                              na.rm = TRUE),
    # regional h values
    genomewide_h_region_3w = weighted.mean(x = genomewide_h_ma3, 
                                            w = population_served, 
                                            na.rm = TRUE),
    spike_h_region_3w = weighted.mean(x = spike_h_ma3, 
                                       w = population_served, 
                                       na.rm = TRUE),
    orf_h_region_3w = weighted.mean(x = orf_h_ma3, 
                                     w = population_served, 
                                     na.rm = TRUE),
    ntd_h_region_3w = weighted.mean(x = ntd_h_ma3, 
                                     w = population_served, 
                                     na.rm = TRUE),
    rbd_h_region_3w = weighted.mean(x = rbd_h_ma3, 
                                     w = population_served, 
                                     na.rm = TRUE),
    cov_mt_2_h_region_3w = weighted.mean(x = cov_mt_2_h_ma3, 
                                          w = population_served, 
                                          na.rm = TRUE),
    # depth
    depth_region_3w = weighted.mean(x = depth_ma3,
                                    w = population_served,
                                    na.rm = TRUE),
    depth_ntd_region_3w = weighted.mean(x = depth_ntd_ma3,
                                        w = population_served,
                                        na.rm = TRUE)
  )%>%
  ungroup() %>%
  filter(week <= week_max)

# add clinical data
meta_r <- meta %>%
  select(region, county) %>%
  distinct()

clinical_weekly <- left_join(hosp_county_weekly, 
                             cases_county_weekly, 
                             by = c("county", "week"))


clinical_region <- clinical_weekly %>%
  left_join(meta_r, by = c("county")) %>%
  group_by(region, week) %>%
  summarize(cases = sum(cases, na.rm = TRUE),
            hospitalizations = sum(hospitalizations, na.rm = TRUE)
  ) %>%
  ungroup()

# add regional incidence
per100 <- county_region %>%
  select(region, region_pop) %>%
  distinct()
per100$region[per100$region == "Capital Region"] <- "Capital"

clinical_region <- left_join(clinical_region, per100, by = c("region"))
clinical_region$case_incidence <- clinical_region$cases * 
  (100000/clinical_region$region_pop)
clinical_region$hosp_incidence <- clinical_region$hospitalizations * 
  (100000/clinical_region$region_pop)

dat_region <- left_join(dat_region, clinical_region, by = c("week", "region"))

# replace na hosp with 0
dat_region$hospitalizations <- ifelse(is.na(dat_region$hospitalizations), 
                                      0, 
                                      dat_region$hospitalizations)
dat_region$hosp_incidence<- ifelse(is.na(dat_region$hosp_incidence), 
                                    0, 
                                    dat_region$hospitalizations)

# merge variant count
dat_region <- left_join(dat_region, var_unique_region, by = c("week", "region"))
dat_region <- dat_region %>%
  filter(!is.na(region))

# check for dups
dups <- dat_region %>%
  group_by(region, week) %>%
  filter(n() > 1)

# --------------------------------------
# STATEWIDE LEVEL COMBINED DATA FILE
# --------------------------------------

# note, grouping without nyc, so add new field for state or nyc
nyc <- c("New York", "Bronx", "Queens", "Kings", "Richmond")
meta_s <- meta %>%
  select(facility_id, population_served)
dat_state <- left_join(dat_sewershed, meta_s, 
                           by = c("facility_id"))
dat_state$group <- ifelse(dat_state$county %in% nyc, "NYC", "State")

# statewide average
dat_state <- dat_state %>%
  arrange(week) %>%
  group_by(week, group)%>%
  summarise(
    # statewide pi
    genomewide_pi_state_3w = weighted.mean(x = genomewide_pi_ma3, 
                                           w = population_served, 
                                           na.rm = TRUE),
    spike_pi_state_3w = weighted.mean(x = spike_pi_ma3, 
                                      w = population_served, 
                                      na.rm = TRUE),
    orf_pi_state_3w = weighted.mean(x = orf_pi_ma3, 
                                    w = population_served, 
                                    na.rm = TRUE),
    ntd_pi_state_3w = weighted.mean(x = ntd_pi_ma3, 
                                    w = population_served, 
                                    na.rm = TRUE),
    rbd_pi_state_3w = weighted.mean(x = rbd_pi_ma3, 
                                    w = population_served, 
                                    na.rm = TRUE),
    cov_mt_2_pi_state_3w = weighted.mean(x = cov_mt_2_pi_ma3, 
                                         w = population_served, 
                                         na.rm = TRUE),
    # statewide concentration
    mean_sars2_conc_state_3w = weighted.mean(x = mean_sars2_conc_ma3, 
                                             w = population_served, 
                                             na.rm = TRUE),
    # statewide h values
    genomewide_h_state_3w = weighted.mean(x = genomewide_h_ma3, 
                                           w = population_served, 
                                           na.rm = TRUE),
    spike_h_state_3w = weighted.mean(x = spike_h_ma3, 
                                      w = population_served, 
                                      na.rm = TRUE),
    orf_h_state_3w = weighted.mean(x = orf_h_ma3, 
                                    w = population_served, 
                                    na.rm = TRUE),
    ntd_h_state_3w = weighted.mean(x = ntd_h_ma3, 
                                    w = population_served, 
                                    na.rm = TRUE),
    rbd_h_state_3w = weighted.mean(x = rbd_h_ma3, 
                                    w = population_served, 
                                    na.rm = TRUE),
    cov_mt_2_h_state_3w = weighted.mean(x = cov_mt_2_h_ma3, 
                                         w = population_served, 
                                         na.rm = TRUE),
    # depth
    depth_state_3w = weighted.mean(x = depth_ma3,
                                   w = population_served,
                                   na.rm = TRUE),
    depth_ntd_state_3w = weighted.mean(x = depth_ntd_ma3,
                                       w = population_served,
                                       na.rm = TRUE)
  )%>%
  ungroup() %>%
  filter(week <= week_max)

# add case and hosp data
clinical_weekly <- left_join(hosp_county_weekly, 
                             cases_county_weekly, 
                             by = c("county", "week"))
clinical_weekly$group <- ifelse(clinical_weekly$county %in% nyc,
                                "NYC", "State")
clinical_state <- clinical_weekly %>%
  group_by(week,group) %>%
  summarize(cases = sum(cases, na.rm = TRUE),
            hospitalizations = sum(hospitalizations, na.rm = TRUE)
  ) %>%
  ungroup()

# merge
dat_state <- left_join(dat_state, clinical_state, by = c("week", "group"))

# add incidence
dat_state$hosp_incidence <- dat_state$hospitalizations * (100000/19870000)
dat_state$case_incidence <- dat_state$cases * (100000 / 19870000)

# add variant count
dat_state <- left_join(dat_state, var_unique_state, by = c("week", "group"))

# check for dups
dups <- dat_state %>%
  group_by(group, week) %>%
  filter(n() > 1)

# ------------------------------------------------------------------------

# Add essential metadata to the sewershed data (e.g., pcr lab)
meta_s <- meta %>%
  select(facility_id, county, region, wwtp_name, population_served, 
         capacity_mgd, 
         wwtp_longitude,
         wwtp_latitude)

pcr_lab <- quant.data %>%
  select(date, pcr_lab, facility_id) %>%
  mutate(week  = floor_date(date, unit = "weeks", week_start = 7)) %>%
  select(-date) %>%
  distinct() %>%
  group_by(facility_id, week) %>%
  slice(1)

dups <- pcr_lab %>%
  group_by(facility_id, week) %>%
  filter(n() > 1)

dat_sewershed <- dat_sewershed %>%
  select(-county, - region)

# join to dat_sewershed
dat_sewershed <- left_join(dat_sewershed, pcr_lab, by = c("week", "facility_id"))
dat_sewershed <- left_join(dat_sewershed, meta_s, by = c("facility_id"))

dups <- dat_sewershed %>%
  group_by(facility_id, week) %>%
  filter(n() > 1)

# filter out nyc counties and sewersheds
dat_sewershed <- dat_sewershed %>%
  filter(region != "New York City")
dat_county <- dat_county %>%
  filter(region != "New York City")
dat_state <- dat_state %>%
  filter(group != "NYC")

# --------------------------------------
# SAVE THE MERGED DATA FILES
# --------------------------------------

save(dat_sewershed, dat_county, dat_region, dat_state,
     file = "data/combined_data.R")


