library(here)
library(tidyverse)
library(foreign)

here::i_am("code/0_prepare_age.R")


age2010 <- read.dbf(here::here("data", "age", "SADE2010_Age.dbf"), as.is = TRUE)

age2010 <- age2010 %>%
  select(-YEAR) %>%
  rename_with(
    ~ paste0(.x, "_2010"),
    starts_with("age")
  )

age2020 <- read.dbf(here::here("data", "age", "SADE2020_Age.dbf"), as.is = TRUE)

age2020 <- age2020 %>%
  select(-YEAR) %>%
  rename_with(
    ~ paste0(.x, "_2020"),
    starts_with("age")
  )

age_10_20 <- full_join(age2010, age2020, by = "tract20l")
write.csv(age_10_20, here::here("data", "age", "age_10_20.csv"), row.names = FALSE)
