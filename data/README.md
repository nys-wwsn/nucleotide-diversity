# Data description

This folder contains the data used in the analysis and production of figures. There are several sub-folders that contain similar sets of data. Descriptions for these data are available in the ReadMe files within the folders.

## clinical data

This folder contains the data for COVID-19 case reports by county and hospital admissions.
The data were retrieved from the [NYS DOH COVID tracker](https://coronavirus.health.ny.gov/covid-19-data-new-york). 

## Counties-Shoreline

New York State County boundaries from the [US Census Tiger/Line Shapefile database](https://www.census.gov/geographies/mapping-files/time-series/geo/tiger-line-file.html).

## diversity_output

This folder contains a file for each genome region that had diversity measures calculated for them.

## genomewide_pi

This folder contains per sample diversity measures for the genome. The output is produced by the script `seq diversity-1-Pi and Shannon calculation.R`. Only the example output files are included here for five samples.

## other data

This folder contains additional data files related to mapping the sewershed spatial data, aggregated lineage names from the Freyja output, concentration data, and other miscellaneous files. See the folder's ReadMe for mor specific information.

## pi

This folder contains per sample diversity measures per base. The output is produced by the script `seq diversity-1-Pi and Shannon calculation.R`. Only the example output files are included here for five samples.

## windowed_pi

This folder contains per sample windowed diversity measures. The output is produced by the script `seq diversity-1-Pi and Shannon calculation.R`. Only the example output files are included here for five samples.

## combined_data.Rdata

This R data file contains four data files. Use `load("data/combined_data.Rdata")` to read in the files. Once loaded into your R session, 4 data files will appear in your environment. Each file has the same underlying structure for different geographic regions indicated in the name.

| R data object name | Geography |
| :------- | :------: |
| `dat_sewershed` | Sewershed level data |
| `dat_county` | NYS county level data |
| `dat_region` | NYS region level data |
| `dat_state` | NYS statewide data |

Below is a description of the data in the `dat_state` data file, which is the most commonly used in the analysis script.

| Field name | Structure| Description |
| :------- | :------: | :------: |
| week | date | Calendar week that a sample was collected |
| group | character | Indicates the group aggregation of the data |
| genomewide_pi_state_3w | numeric | Three week rolling average Pi value for the genome for the state. Other data files are indicated after the underscore. e.g., _county |
| spike_pi_state_3w | numeric | Three week rolling average Pi value for the spike region for the state. Other data files are indicated after the underscore. e.g., _county |
| ntd_pi_state_3w | numeric | Three week rolling average Pi value for the S1 NTD region for the state. Other data files are indicated after the underscore. e.g., _county |
| rbd_pi_state_3w | numeric | Three week rolling average Pi value for the S1 RBD region for the state. Other data files are indicated after the underscore. e.g., _county |
| cov_mt_2_pi_state_3w | numeric | Three week rolling average Pi value for the 2’ O-Mtase  region for the state. Other data files are indicated after the underscore. e.g., _county |
| mean_sars_conc_state_3w | numeric | Three week rolling average concentration measured in copies per mL |
| genomewide_h_state_3w | numeric | Three week rolling average Shannon H value for the genome for the state. Other data files are indicated after the underscore. e.g., _county |
| spike_h_state_3w | numeric | Three week rolling average Shannon H value for the spike region for the state. Other data files are indicated after the underscore. e.g., _county |
| ntd_h_state_3w | numeric | Three week rolling average Shannon H value for the S1 NTD region for the state. Other data files are indicated after the underscore. e.g., _county |
| rbd_h_state_3w | numeric | Three week rolling average Shannon H value for the S1 RBD region for the state. Other data files are indicated after the underscore. e.g., _county |
| cov_mt_2_h_state_3w | numeric | Three week rolling average Shannon H value for the 2’ O-Mtase  region for the state. Other data files are indicated after the underscore. e.g., _county |
| depth_state_3w | numeric | Three week rolling average depth of read from the variants file value averaged across all samples across the state. Other data files are indicated after the underscore. e.g., _county |
| depth_state_ntd_3w | numeric | Three week rolling average depth of read from S1 NTD region from the variants file value averaged across all samples across the state. Other data files are indicated after the underscore. e.g., _county |
| n_variants_no_thresh_3w_mean | numeric | Three week rolling average number of Freyja variants reported in the state. No prevalence/abundance threshold was used, and all variants are included. |
| n_variants_5_3w_mean | numeric | Three week rolling average number of Freyja variants reported in the state. A five percent threshold was used. |
| cases | numeric | Case counts across the state for positive COVID-19 PCR tests |
| hospitalizations | numeric | Hospitalization counts across the state for individuals admitted to the hospital that also tested positive for COVID-19 |
| case_incidence | numeric | Case counts per 100,000 population across the state for positive COVID-19 PCR tests |
| hosp_incidence | numeric | Hospitalization counts per 100,000 population across the state for individuals admitted to the hospital that also tested positive for COVID-19 |
| hosp_incidence | numeric | Hospitalization counts per 100,000 population across the state for individuals admitted to the hospital that also tested positive for COVID-19 |
