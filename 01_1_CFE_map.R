#####################################################################
# Last updated 2023-03-21 - Alexander Van Nynatten
# Plots the eDNA sampling locations and the Cape Fold Eco-region boundaries 
# Plots locations where reference tissue samples were collected for sequence database
#####################################################################

#####################################################################
## Loads required packages and data
#####################################################################

library(osmdata)
library(sf)
library(tidyverse)
library(ggspatial)
library(rmapshaper)

#####################################################################
## Makes bounding box of the CFE
#####################################################################

big_bbox <- st_bbox(c(xmin = 15.0, xmax = 35.0, ymin = -35, ymax = -20), 
	crs = st_crs(4326)) %>%
	st_as_sfc()

#####################################################################
## Loads location data where eDNA samples were collected
#####################################################################

sample_sf <- read_csv('./00_raw_data/sample_metadata/location_metadata.csv') %>%
	st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
	mutate(Project_code = ifelse(grepl('BK', Location_code), 'BK', 'Bos'))

## Loads river and drainage polygon data 
# data from: https://www.dws.gov.za/iwqs/gis_data/river/rivs500k.aspx

# Reads map data files
river_sf <- read_sf('./00_raw_data/map_data/wriall500/wriall500.shp')
drainage_sf <- read_sf('./00_raw_data/map_data/hca_4/hca_4.shp')

#####################################################################
## Subsamples the data for the CFE region in South Africa
#####################################################################

# Subsamples the drainages in the Cape Fold Ecoregion
crop.drainage_sf <- drainage_sf %>%
	filter(PRIMARY_ %in% c('G', 'E', 'H', 'J', 'K', 'M', 'N', 'L')) %>%
	st_transform(crs = 4326)

# Joins the sub-drainages for more simplified plotting
all.crop.drainage_sf <- crop.drainage_sf %>%
	group_by(PRIMARY_) %>%
    summarize(geometry = st_union(st_buffer(geometry, 100)))

# Subsamples just the sub-drainages in this study
sampled.crop.drainage_sf <- crop.drainage_sf %>%
	filter(QUATERNARY %in% unique(st_intersection(crop.drainage_sf, sample_sf)$QUATERNARY))

# Crops rivers in CFE
crop.river_sf <- river_sf %>%
	st_transform(crs = 4326) %>%
	st_crop(big_bbox) %>%
	st_intersection(crop.drainage_sf)

#####################################################################
## Small map of Africa to show location of CFE
#####################################################################

# Gets global data for countries (small scale)
sm_world <- rnaturalearth::ne_countries(scale = 'small', returnclass = 'sf')

# Subsamples Africa
Af_map <- sm_world %>%
	filter(continent %in% 'Africa') %>%
	st_make_valid() %>%
	st_transform(crs = 4326)

#####################################################################
# Loads location data where reference sequence samples were collected
#####################################################################

reference_sf <- read_csv('./00_raw_data/sample_metadata/reference_metadata.csv') %>%
  drop_na() %>%
  st_as_sf(coords = c("long", "lat"), crs = 4326) %>%
  st_intersection(crop.drainage_sf)

#####################################################################
## Plots Figure 1
#####################################################################

# Map of the location of the CFE relative the African continent
ggplot() +
   geom_sf(data = Af_map, fill = 'grey90', size = 0.5) +
   geom_sf(data = all.crop.drainage_sf, fill = '#D1EFDF', size = 1) +
   theme_test()

ggsave('./01_sample_map/Fig1a_Africa.pdf')

# Map of the CFE
ggplot() +
   geom_sf(data = all.crop.drainage_sf, aes(fill = PRIMARY_), size = 1) +
   geom_sf(data = crop.river_sf, colour = '#497186') +
   geom_sf(data = sampled.crop.drainage_sf, fill = NA, colour = 'grey20', size = 1) +
   geom_sf(data = reference_sf, aes(fill = Species), size = 3, shape = 22) +
   theme_test() +
  theme(panel.background = element_rect(fill = 'white'),
	panel.grid = element_line(color = alpha("black", 0.25), size = 0.25)) +
	annotation_scale(location = "br", width_hint = 0.4) +
	annotation_north_arrow(location = "tr", which_north = "true",
		pad_x = unit(0.2, "in"), pad_y = unit(0.2, "in"),
		style = north_arrow_fancy_orienteering)

ggsave('./01_sample_map/Fig1a_CFE.pdf')

#####################################################################
## Obtains data for smaller rivers for inset map from open street maps
# Not available in the 500k scale rivers of South Africa
#####################################################################

# Downloads waterway elements
OSM_rivers <- opq(bbox = st_bbox(sampled.crop.drainage_sf), timeout = 100) %>%
add_osm_feature('waterway') %>%
osmdata_sf() %>% 
unname_osmdata_sf()

# Drops waterway linestring elements from the drainages of interest
crop.OSM_rivers <- OSM_rivers$osm_lines[c('geometry', 'waterway')] %>%
	filter(!is.na(waterway)) %>%
	st_make_valid() %>%
	mutate(waterway = 'waterway') %>%
	rename(Type = waterway) %>%
	st_transform(crs = 4326) %>%
	st_intersection(sampled.crop.drainage_sf)

#####################################################################
## Plots Figure 1bc
#####################################################################

Map of Bos and Blindekloof Rivers
ggplot() +
   geom_sf(data = sampled.crop.drainage_sf, fill = NA, colour = 'grey20', size = 1) +
   geom_sf(data = crop.OSM_rivers, colour = '#497186', size = 0.2) +
   geom_sf(data = sample_sf, colour = "black", size = 0.8) +
   theme_test() + 
   coord_sf(expand = 0) + 
  annotation_scale(location = "bl", width_hint = 0.025) + 
  annotation_scale(location = "br", width_hint = 0.025)
  
ggsave('./01_sample_map/Fig1bc_DrainageMap.pdf')

#####################################################################

save.image(file="./01_sample_map/Fig1.RData")
