#
# COMPARE DIVERSITY TO POPULATION DENSITY PER COUNTY
#


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
library(units)

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

# load county spatial data
ny_counties <- st_read("data/Counties-Shoreline/Counties_Shoreline.shp") %>%
  mutate(county = NAME)
ny_counties$county[ny_counties$county == "St Lawrence"] <- "St Lawrence"

# calculate area
ny_counties$area <- st_area(ny_counties)

# change to km2
ny_counties$area <- units::set_units(ny_counties$area, km^2)

# select columns and drop unites
ny_counties <- ny_counties %>%
  st_drop_geometry() %>%
  select(county, area) %>%
  drop_units()

# COUNTY DATA - POPULATION/SQ KM ~ DIVERSITY

# join area
dat_county <- left_join(dat_county,
                        ny_counties,
                        by = c("county"))

# pop density
dat_county$pop_density <- dat_county$county_pop / dat_county$area

# plot
plot <- 
  ggplot(data = dat_county,
       aes(x = pop_density,
           y = ntd_pi_county_3w))+
  geom_point()+
  stat_cor(method ="spearman")+
  theme_dth_1+
  labs(title = "Population density ~ diversity",
       x = "Population density (persons per square kilometer)",
       y = expression(Pi[ww]))

# save
png("Supplemental files/Figure S unknown 2.png",
    units = "in",
    width = 6, height = 4,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
plot 
showtext_end()
dev.off()