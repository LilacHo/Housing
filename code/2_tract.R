library(here)
library(tidyverse)

here::i_am("code/2_tract.R")

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

# -----------------------------
# 4. FUNCTION: FIT 2010 MODEL, GET 2010 RESIDUALS, 2020 EXPECTED, 2020 PSEUDO-RESIDUALS
# -----------------------------
# We do this separately within each season x buffer combination.
# That is usually the safest move if your "above" analysis was also split by these dimensions.

run_baseline_change_models <- function(dat) {
  
  # Need enough rows to estimate the 2010 models
  if (nrow(dat) < 10) {
    return(
      dat %>%
        mutate(
          exp_greenness_2020   = NA_real_,
          exp_impervious_2020  = NA_real_,
          resid_greenness_2010 = NA_real_,
          resid_impervious_2010 = NA_real_,
          pseudo_resid_greenness_2020 = NA_real_,
          pseudo_resid_impervious_2020 = NA_real_,
          greenness_change = NA_real_,
          impervious_change = NA_real_
        )
    )
  }
  
  # 2010 baseline models
  m_green_2010 <- lm(greenness_2010 ~ housing_2010, data = dat)
  m_imperv_2010 <- lm(impervious_2010 ~ housing_2010, data = dat)
  
  # 2010 actual residuals
  dat$resid_greenness_2010  <- resid(m_green_2010)
  dat$resid_impervious_2010 <- resid(m_imperv_2010)
  
  # 2020 expected values using 2010 models
  dat$exp_greenness_2020 <- predict(
    m_green_2010,
    newdata = dat %>% transmute(housing_2010 = housing_2020)
  )
  
  dat$exp_impervious_2020 <- predict(
    m_imperv_2010,
    newdata = dat %>% transmute(housing_2010 = housing_2020)
  )
  
  # 2020 pseudo-residuals using 2010 relationship
  dat$pseudo_resid_greenness_2020  <- dat$greenness_2020  - dat$exp_greenness_2020
  dat$pseudo_resid_impervious_2020 <- dat$impervious_2020 - dat$exp_impervious_2020
  
  # Change variable = 2020 pseudo-residual - 2010 residual
  dat$greenness_change  <- dat$pseudo_resid_greenness_2020  - dat$resid_greenness_2010
  dat$impervious_change <- dat$pseudo_resid_impervious_2020 - dat$resid_impervious_2010
  
  dat
}

change_df <- wide %>%
  group_by(season, buffer) %>%
  group_modify(~ run_baseline_change_models(.x)) %>%
  ungroup()

# -----------------------------
# 5. SAVE MODEL SUMMARIES FOR THE 2010 BASELINE MODELS
# -----------------------------
baseline_model_summaries <- wide %>%
  group_by(season, buffer) %>%
  group_modify(~ {
    dat <- .x
    if (nrow(dat) < 10) return(tibble())
    
    m_green <- lm(greenness_2010 ~ housing_2010, data = dat)
    m_imperv <- lm(impervious_2010 ~ housing_2010, data = dat)
    
    bind_rows(
      tidy(m_green) %>% mutate(outcome = "greenness"),
      tidy(m_imperv) %>% mutate(outcome = "impervious")
    ) %>%
      select(outcome, everything())
  }) %>%
  ungroup()

glance_model_summaries <- wide %>%
  group_by(season, buffer) %>%
  group_modify(~ {
    dat <- .x
    if (nrow(dat) < 10) return(tibble())
    
    m_green <- lm(greenness_2010 ~ housing_2010, data = dat)
    m_imperv <- lm(impervious_2010 ~ housing_2010, data = dat)
    
    bind_rows(
      glance(m_green) %>% mutate(outcome = "greenness"),
      glance(m_imperv) %>% mutate(outcome = "impervious")
    )
  }) %>%
  ungroup()

# -----------------------------
# 6. HYPOTHETICAL DEMOGRAPHIC CHANGE SETUP
# -----------------------------
# You said you do not yet have racial and age information.
# Below is example code assuming that for each GEOID20 you later merge in columns such as:
#
# pct_white_2010, pct_white_2020
# pct_black_2010, pct_black_2020
# pct_latino_2010, pct_latino_2020
# pct_asian_2010, pct_asian_2020
# pct_age65plus_2010, pct_age65plus_2020
# median_age_2010, median_age_2020
#
# These can live in a separate dataframe called demo_df.

# Example skeleton:
# demo_df <- data.frame(
#   GEOID20 = unique(change_df$GEOID20),
#   pct_white_2010 = runif(length(unique(change_df$GEOID20)), 20, 80),
#   pct_white_2020 = runif(length(unique(change_df$GEOID20)), 20, 80),
#   pct_black_2010 = runif(length(unique(change_df$GEOID20)), 0, 30),
#   pct_black_2020 = runif(length(unique(change_df$GEOID20)), 0, 30),
#   pct_latino_2010 = runif(length(unique(change_df$GEOID20)), 0, 40),
#   pct_latino_2020 = runif(length(unique(change_df$GEOID20)), 0, 40),
#   pct_asian_2010 = runif(length(unique(change_df$GEOID20)), 0, 40),
#   pct_asian_2020 = runif(length(unique(change_df$GEOID20)), 0, 40),
#   pct_age65plus_2010 = runif(length(unique(change_df$GEOID20)), 5, 25),
#   pct_age65plus_2020 = runif(length(unique(change_df$GEOID20)), 5, 30),
#   median_age_2010 = runif(length(unique(change_df$GEOID20)), 28, 50),
#   median_age_2020 = runif(length(unique(change_df$GEOID20)), 28, 55)
# )

# Merge demographic data to the change metrics.
# If demographic data are at GEOID20 level only, this duplicates the same demographic change
# across season/buffer rows. That is OK if you later choose one season/buffer combination
# or run separate models by season/buffer.
analysis_df <- change_df %>%
  left_join(demo_df, by = "GEOID20") %>%
  mutate(
    change_pct_white      = pct_white_2020      - pct_white_2010,
    change_pct_black      = pct_black_2020      - pct_black_2010,
    change_pct_latino     = pct_latino_2020     - pct_latino_2010,
    change_pct_asian      = pct_asian_2020      - pct_asian_2010,
    change_pct_age65plus  = pct_age65plus_2020  - pct_age65plus_2010,
    change_median_age     = median_age_2020     - median_age_2010
  )

# -----------------------------
# 7. RUN DEMOGRAPHIC-CHANGE REGRESSIONS
# -----------------------------
# These regress demographic change on:
#   greenness_change + impervious_change
#
# Interpretation:
#   - more negative greenness_change = more loss of greenness relative to what 2020 housing density
#     would have predicted under the 2010 relationship
#   - more positive impervious_change = more imperviousness relative to that same baseline expectation
#
# You may want controls too (baseline demographics, 2010 housing density, metro fixed effects, etc.),
# but below is the simple version matching your request.

demo_outcomes <- c(
  "change_pct_white",
  "change_pct_black",
  "change_pct_latino",
  "change_pct_asian",
  "change_pct_age65plus",
  "change_median_age"
)

run_demo_models <- function(dat, outcomes) {
  map_dfr(outcomes, function(y) {
    f <- as.formula(
      paste0(y, " ~ greenness_change + impervious_change")
    )
    
    mod <- lm(f, data = dat)
    
    tidy(mod) %>%
      mutate(outcome = y) %>%
      select(outcome, everything())
  })
}

# Option A: pooled across all rows
demo_model_results_pooled <- run_demo_models(
  dat = analysis_df %>%
    filter(
      !is.na(greenness_change),
      !is.na(impervious_change)
    ),
  outcomes = demo_outcomes
)

# Option B: separate by season and buffer
demo_model_results_by_group <- analysis_df %>%
  group_by(season, buffer) %>%
  group_modify(~ run_demo_models(.x, demo_outcomes)) %>%
  ungroup()

# -----------------------------
# 8. OPTIONAL: ADD CONTROLS
# -----------------------------
# Often a stronger specification would include baseline demographics and baseline housing.
# Example:
#
# change_pct_white ~ greenness_change + impervious_change + pct_white_2010 + housing_2010
#
# Here is a helper if you want that approach.

run_demo_models_with_controls <- function(dat, outcomes) {
  map_dfr(outcomes, function(y) {
    
    baseline_var <- case_when(
      y == "change_pct_white"     ~ "pct_white_2010",
      y == "change_pct_black"     ~ "pct_black_2010",
      y == "change_pct_latino"    ~ "pct_latino_2010",
      y == "change_pct_asian"     ~ "pct_asian_2010",
      y == "change_pct_age65plus" ~ "pct_age65plus_2010",
      y == "change_median_age"    ~ "median_age_2010",
      TRUE                        ~ NA_character_
    )
    
    rhs <- c("greenness_change", "impervious_change", baseline_var, "housing_2010")
    rhs <- rhs[!is.na(rhs)]
    
    f <- as.formula(paste(y, "~", paste(rhs, collapse = " + ")))
    mod <- lm(f, data = dat)
    
    tidy(mod) %>%
      mutate(outcome = y) %>%
      select(outcome, everything())
  })
}

demo_model_results_controls <- analysis_df %>%
  filter(!is.na(greenness_change), !is.na(impervious_change)) %>%
  run_demo_models_with_controls(outcomes = demo_outcomes)

# -----------------------------
# 9. OPTIONAL: SIMPLE INTERPRETATION EXAMPLE
# -----------------------------
# Example: extract coefficient for greenness_change from pooled white-share model

white_model <- lm(
  change_pct_white ~ greenness_change + impervious_change,
  data = analysis_df
)

summary(white_model)

# Interpretation template:
# If coef(greenness_change) < 0:
#   neighborhoods with more negative greenness_change
#   (i.e., they became greener less than expected, or lost greenness relative to expected)
#   tended to have increases in white share if the outcome is change_pct_white and signs line up accordingly.
#
# More concretely:
# - Negative x-values on greenness_change = "greenness dimming" relative to the 2010 baseline relationship
# - Positive coefficient on greenness_change for an outcome means demographic change rises as greenness_change rises
# - Negative coefficient means demographic change rises as greenness_change falls

# -----------------------------
# 10. EXPORT RESULTS
# -----------------------------
write.csv(change_df, "greenness_impervious_change_metrics.csv", row.names = FALSE)
write.csv(baseline_model_summaries, "baseline_2010_model_coefficients.csv", row.names = FALSE)
write.csv(glance_model_summaries, "baseline_2010_model_fitstats.csv", row.names = FALSE)
write.csv(demo_model_results_pooled, "demographic_change_models_pooled.csv", row.names = FALSE)
write.csv(demo_model_results_by_group, "demographic_change_models_by_season_buffer.csv", row.names = FALSE)
write.csv(demo_model_results_controls, "demographic_change_models_with_controls.csv", row.names = FALSE)