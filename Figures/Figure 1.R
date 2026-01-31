###########
###########
#### GENOMEWIDE DIVERSITY OF SARS-COV-2 SEQUENCES OF WASTEWATER SAMPLES        #
###########
###########

# Script author: Dustin T. Hill

# Created 2025-04-25
# Last updated 2025-09-25

# Main analysis - Figure 1

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

# Load functions
source("seq diversity - functions.R")

# Load plot themes
source("seq diversity - plot themes.R")

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


# --------------
# Figure 1
# --------------
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
          color = "grey65",
          lwd = 1
          
  )+
  geom_sf(data = detection_rate,
          aes(color = as.numeric(`Detection rate`),
              size = capacity_mgd),
          alpha = 0.7)+
  geom_sf(data = detection_rate,
          aes(
            size = capacity_mgd),
          shape = 1)+
  theme_dth_maps +
  scale_color_gradientn(colors = rev(MetBrewer::met.brewer("Hiroshige")),
                        name = "Detection\nrate")+
  theme(legend.position.inside = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.direction = "horizontal")+
  labs(size = "WWTP size\n(mgd)",
       title = "SARS-CoV-2 detection rate in wastewater")+
  theme(legend.position = c(0.6, 0.3),
        legend.background = element_blank(),
        legend.box.just = "left")+
  scale_size(range = c(0, 20))

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

pal_1 <- met.brewer(name = "Hiroshige", n = 6)

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
  scale_fill_manual(values = c("samples_no_seq" = pal_1[1],
                               "samples_sequenced" = pal_1[6]),
                    labels = c("Samples not sequenced",
                               "Samples sequenced"))

samples_plot 



# panel
ggarrange(map_dr, samples_plot , nrow = 1, labels = c("A","B"))

# save map
png("Figures/Figure 1.png",
    units = "in",
    width = 12, height = 6,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
ggarrange(map_dr, samples_plot , nrow = 1, labels = c("A","B"))
showtext_end()
dev.off()

