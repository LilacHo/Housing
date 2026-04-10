library(here)
library(sf)

here::i_am("code/0_1_prepare_block.R")

## Add area to block ####
# Import shapefile
# Data sourse: https://geo.wa.gov/datasets/wa-ofm::ofm-saep-block20/about
block <- st_read(here::here("data", "saep_block20", "saep_block20_shapefile","block20.shp"))

# check CRS
st_crs(block) # -> NAD83(HARN) / Washington South (ftUS) 

# reproject to meter
block_proj <- st_transform(block, 2856) # EPSG:2856; NAD83(HARN) / Washington South (m)

block_proj$area_m2 <- st_area(block_proj)
block_proj$area_ha <- block_proj$area_m2 / 10000

# output csv
block_df <- st_drop_geometry(block_proj)

write.csv(block_df,
          here::here("data", "saep_block20", "block20_with_area.csv"),
          row.names = FALSE)

