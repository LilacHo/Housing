library(here)
library(sf)

here::i_am("code/0_1_prepare_tract.R")

## Add area to tract ####
# Import shapefile
# Data source: https://ofm.wa.gov/data-research/population-demographics/estimates/small-area/
tract <- st_read(here::here("data", "saep_tract20", "tract20.shp"))

# check CRS
st_crs(tract) # NAD83(HARN) / Washington South (ftUS)

# reproject to meter (EPSG:2856; NAD83(HARN) / Washington South (m))
tract_proj <- st_transform(tract, 2856)

tract_proj$area_m2 <- st_area(tract_proj)
tract_proj$area_ha <- tract_proj$area_m2 / 10000

# output csv with geometry dropped
tract_df <- st_drop_geometry(tract_proj)

write.csv(tract_df,
          here::here("data", "saep_tract20", "tract20_with_area.csv"),
          row.names = FALSE)
