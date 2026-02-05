###
### Supplemental figure - map showing the geographic scales
###
###


### Created February 5, 2026
### Updated February 5, 2026

# Packages
library(dplyr)
library(sf)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(ggpubr)

# --------------------------------------
# PLOT THEMES
# --------------------------------------
source("seq diversity - plot themes.R")

# Data
# combined data files for each geography
load(file = "data/combined_data.Rdata")

# sewershed locations from our study
sewersheds <- dat_sewershed %>%
  select(wwtp_latitude, wwtp_longitude, facility_id, capacity_mgd) %>%
  distinct() %>%
  st_as_sf(coords =  c("wwtp_longitude", "wwtp_latitude"),
           crs = "+proj=longlat +zone=18 +datum=WGS84") %>% 
  st_transform(st_crs(ny_counties))

# county boundaries
ny_counties <- st_read("data/Counties-Shoreline/Counties_Shoreline.shp") %>%
  mutate(county = NAME)
ny_counties$county[ny_counties$county == "St. Lawrence"] <- "St Lawrence"

# create region boundaries object
county_region <- read.csv("data/other data/County_Region.csv")
county_region$county[county_region$county == "St. Lawrence"] <- "St Lawrence"
regions <- left_join(ny_counties, county_region, by = c("county"))

# indicator for NYC or not
regions$nyc <- ifelse(regions$region == "New York City",
                      "Not represented in the current study",
                      "Represented in this study")

# state boundary object
state <- ny_counties %>%
  summarize()


# Plot 1 - sewersheds
p1 <- 
  ggplot()+
  geom_sf(data = state, fill = MetBrewer::met.brewer("Hiroshige", 1),
          color = "grey")+
  geom_sf(data = sewersheds, aes(size = capacity_mgd),
          alpha = 0.5)+
  theme_dth_maps+
  theme(legend.position = "bottom",
        plot.subtitle = element_text(hjust = 0.5))+
  labs(size = "WWTP maximum capacity (mgd)",
       subtitle = "194 WWTP sampling sites",
       title = "Sewershed aggregation")

# plot 2 - county boundaries
p2 <- ggplot()+
  geom_sf(data = regions, aes(fill = nyc),
          color = "white")+
  theme_dth_maps+
  labs(subtitle = "57 counties with data in this study",
       title = "County aggregation",
       fill = "")+
  scale_fill_manual(values = MetBrewer::met.brewer(name = "Hiroshige", n = 2))+
  theme(legend.position = "bottom",
        plot.subtitle = element_text(hjust = 0.5))
p2

# regions
p3 <- ggplot()+
  geom_sf(data = regions, aes(fill = region),
          color = "white")+
  theme_dth_maps+
  labs(subtitle = "The 10 regions in New York State",
       title = "Regional aggregation")+
  theme(legend.position = "none",
        plot.subtitle = element_text(hjust = 0.5))+
  scale_fill_manual(values = MetBrewer::met.brewer(name = "Hiroshige",
                                                   n = 10))
p3

# state boundary
p4 <- ggplot()+
  geom_sf(data = state, fill = MetBrewer::met.brewer("Hiroshige", 1), color = "grey")+
  theme_dth_maps+
  labs(subtitle = "New York State",
       title  = "Statewide aggregation")+
  theme(plot.subtitle = element_text(hjust =0.5))
p4

ggarrange(p1,p2,p3,p4, nrow = 2, ncol = 2,
             labels = c("A", "B", "C", "D"))

# save
png("Supplemental files/Figure S8.png",
    units = "in",
    width = 8, height = 6,
    res = 600)
showtext::showtext_auto()
showtext::showtext_opts(dpi = 600)
ggarrange(p1,p2,p3,p4, nrow = 2, ncol = 2,
          labels = c("A", "B", "C", "D"))
showtext_end()
dev.off()