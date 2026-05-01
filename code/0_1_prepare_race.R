library(here)
library(tidyverse)
library(foreign)

here::i_am("code/0_1_prepare_race.R")


race2010 <- read.dbf(here::here("data", "race", "SADE2010_race.dbf"), as.is = TRUE)

race2010 <- race2010 %>%
  select(-YEAR) %>%
  rename(
    AIAN2010 = AIAN,
    ASIAN2010 = ASIAN,
    Black2010 = Black,
    Two2010 = Two,
    White2010 = White,
    NHOPI2010 = NHOPI
  )

race2020 <- read.dbf(here::here("data", "race", "SADE2020_race.dbf"), as.is = TRUE)

race2020 <- race2020 %>%
  select(-YEAR) %>%
  rename(
    AIAN2020 = AIAN,
    ASIAN2020 = ASIAN,
    Black2020 = Black,
    Two2020 = Two,
    White2020 = White,
    NHOPI2020 = NHOPI
  )

race_10_20 <- full_join(race2010, race2020, by = "tract20l")
write.csv(race_10_20, here::here("data", "race", "race_10_20.csv"), row.names = FALSE)
