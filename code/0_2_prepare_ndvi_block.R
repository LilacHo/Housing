library(here)
library(tidyverse)

here::i_am("code/0_2_prepare_ndvi_block.R")

## NDVI investigation ####
ndvi_block <- read.csv(here::here("data","ndvi_block_final.csv"))

cor(ndvi_block$POP, ndvi_block$HU)
# 0.937783

cor(ndvi_block$PCT_BUILT, ndvi_block$NDVI, use = "complete.obs") 
# -0.4047123

# compare NDVI across buffers (500, 1000, 2000) within the same GEOID20
df_ndvi_wide <- ndvi_block %>%
  select(GEOID20, year, season, buffer, NDVI) %>% 
  pivot_wider(names_from = buffer, values_from = NDVI)

# Drop missing values
df_ndvi_wide_drop_na <- df_ndvi_wide %>% drop_na(`500`, `1000`, `2000`)

# Pairwise comparisons (paired t-tests)
t.test(df_ndvi_wide_drop_na$`500`, df_ndvi_wide_drop_na$`1000`, paired = TRUE)
# reject the null hypothesis
# NDVI differs between buffer sizes 500 and 1000

t.test(df_ndvi_wide_drop_na$`500`, df_ndvi_wide_drop_na$`2000`, paired = TRUE)
# reject the null hypothesis
# NDVI differs between buffer sizes 500 and 2000

t.test(df_ndvi_wide_drop_na$`1000`, df_ndvi_wide_drop_na$`2000`, paired = TRUE)
# reject the null hypothesis
# NDVI differs between buffer sizes 1000 and 2000

## Block ####
block_df <- read.csv(here::here("data", "saep_block20", "block20_with_area.csv"))

block_df <- block_df %>%
  select(GEOID20, area_ha)

# Add area to ndvi dataframe
ndvi_block2 <- ndvi_block %>%
  left_join(block_df, by = "GEOID20")

# Calculate density of POP & HU
ndvi_block2 <- ndvi_block2 %>%
  mutate(
    POP_density = POP / area_ha,
    HU_density  = HU  / area_ha
  )

# examine the distribution
hist(ndvi_block2$POP_density)
hist(log(ndvi_block2$POP_density))
hist(ndvi_block2$HU_density)
hist(log(ndvi_block2$HU_density))

# output
write.csv(ndvi_block2, here::here("data", "block_final.csv"), row.names = FALSE)
