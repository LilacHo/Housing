library(here)
library(tidyverse)

here::i_am("code/0_2_prepare_ndvi_tract.R")

## NDVI investigation (tract level) ####
ndvi_tract <- read.csv(here::here("data", "ndvi_tract_final.csv"))

# Normalize column names to mirror block dataset conventions
# (source file uses "weighted_nvdi" [sic] and "weighted_built")
ndvi_tract <- ndvi_tract %>%
  rename(NDVI      = weighted_nvdi,
         PCT_BUILT = weighted_built)

# Built / greenness correlation
cor(ndvi_tract$PCT_BUILT, ndvi_tract$NDVI, use = "complete.obs")

# compare NDVI across buffers (500, 1000, 2000) within the same tract-year-season
df_ndvi_wide <- ndvi_tract %>%
  select(tract, year, season, buffer, NDVI) %>%
  pivot_wider(names_from = buffer, values_from = NDVI)

df_ndvi_wide_drop_na <- df_ndvi_wide %>% drop_na(`500`, `1000`, `2000`)

# Pairwise comparisons (paired t-tests)
t.test(df_ndvi_wide_drop_na$`500`,  df_ndvi_wide_drop_na$`1000`, paired = TRUE)
t.test(df_ndvi_wide_drop_na$`500`,  df_ndvi_wide_drop_na$`2000`, paired = TRUE)
t.test(df_ndvi_wide_drop_na$`1000`, df_ndvi_wide_drop_na$`2000`, paired = TRUE)


## Tract ####
tract_df <- read.csv(here::here("data", "saep_tract20", "tract20_with_area.csv"))

# Keep identifiers, area, and the two decennial housing & population columns
tract_df <- tract_df %>%
  mutate(tract = as.numeric(GEOID20)) %>%
  select(tract, area_ha,
         POP2010, POP2020,
         HU2010,  HU2020)

# Reshape housing/population to long form so each tract-year has one row
tract_long <- tract_df %>%
  pivot_longer(
    cols = c(POP2010, POP2020, HU2010, HU2020),
    names_to  = c(".value", "year"),
    names_pattern = "(POP|HU)(\\d{4})"
  ) %>%
  mutate(year = as.integer(year))

# Merge NDVI with tract HU/POP for the matching year
ndvi_tract2 <- ndvi_tract %>%
  left_join(tract_long, by = c("tract", "year"))

# Calculate density of POP & HU
ndvi_tract2 <- ndvi_tract2 %>%
  mutate(
    POP_density = POP / area_ha,
    HU_density  = HU  / area_ha
  )

# examine the distribution
hist(ndvi_tract2$POP_density)
hist(log(ndvi_tract2$POP_density))
hist(ndvi_tract2$HU_density)
hist(log(ndvi_tract2$HU_density))

# output
write.csv(ndvi_tract2, here::here("data", "tract_final.csv"), row.names = FALSE)
