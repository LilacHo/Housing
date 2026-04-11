library(here)
library(tidyverse)

here::i_am("code/1_block.R")

# ============================================================
# PSEUDO-RESIDUAL CHANGE ANALYSIS, 2010 -> 2020
# Greenness and imperviousness relative to housing density
# ============================================================

# -----------------------------
# 0. Packages
# -----------------------------
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(stringr)

# -----------------------------
# 1. USER INPUTS / ASSUMPTIONS
# -----------------------------

ndvi_block2 <- read.csv(here::here("data", "block_final.csv"))

# Dataset contains:
# `GEOID20`
# `year`: 2010, 2020
# `season`: Dormant_Season, Growing_Season
# `buffer`: 500, 1000, 2000
# greenness: `NDVI`
# impervious: `PCT_BUILT`
# `POP_density`
# `HU_density`

# -----------------------------
# 2. BASIC CLEANING
# -----------------------------

df <- ndvi_block2 %>%
  mutate(
    across(all_of(c(PCT_BUILT, NDVI, POP_density, HU_density)), as.numeric),
    year = as.integer(year)
  ) %>%
  filter(year %in% c(2010, 2020)) %>%
  filter(!is.na(.data[[PCT_BUILT]]),
         !is.na(.data[[NDVI]]),
         !is.na(.data[[POP_density]]),
         !is.na(.data[[HU_density]])) %>%
  mutate(
    log_POP_density = log(POP_density),
    log_HU_density  = log(HU_density)
  )

id_vars <- c("GEOID20", "season", "buffer")

housing_var    <- "log_HU_density"
population_var <- "log_POP_density"
greenness_var  <- "NDVI"   
impervious_var <- "PCT_BUILT"       


# -----------------------------
# 3. WIDE DATA FOR 2010 / 2020
# -----------------------------
# This makes one row per GEOID20-season-buffer combination.
wide <- df2 %>%
  select(all_of(id_vars), year, 
         !!sym(housing_var), !!sym(population_var),
         !!sym(greenness_var), !!sym(impervious_var)) %>%
  rename(
    housing     = !!sym(housing_var),
    population  = !!sym(population_var),
    greenness   = !!sym(greenness_var),
    impervious  = !!sym(impervious_var)
  ) %>%
  pivot_wider(
    names_from  = year,
    values_from = c(housing, population, greenness, impervious),
    names_sep   = "_"
  ) %>%
  filter(
    !is.na(housing_2010), !is.na(housing_2020),
    !is.na(population_2010), !is.na(population_2020),
    !is.na(greenness_2010), !is.na(greenness_2020),
    !is.na(impervious_2010), !is.na(impervious_2020)
  )

