library(here)
library(tidyverse)

here::i_am("code/0_3_prepare_tract_age_race.R")

# ============================================================
# Prepare tract-level age and race data for analysis.
#
# Inputs:
#   data/age_10_20.csv   (from 0_1_prepare_age.R)
#   data/race_10_20.csv  (from 0_1_prepare_race.R)
#
# Outputs:
#   data/age_pct_tract.csv  — % in 4 age groups per tract-year
#   data/race_pct_tract.csv    — % in 6 race categories per tract-year
#
# Age groups:  0-19, 20-29, 30-59, 60+
# Race groups: AIAN, Asian, Black, NHOPI, Two or More, White
# ============================================================

# -------------------------------------------------------
# 1. AGE: 18 five-year bands → 4 groups
# -------------------------------------------------------
age <- read.csv(here::here("data", "age", "age_10_20.csv"))

age_bands <- c("age0_4","age5_9","age10_14","age15_19",
               "age20_24","age25_29",
               "age30_34","age35_39","age40_44","age45_49","age50_54","age55_59",
               "age60_64","age65_69","age70_74","age75_79","age80_84","age85_up")

# Band indices for each group:
#   0-19:  1-4   (age0_4 .. age15_19)
#   20-29: 5-6   (age20_24, age25_29)
#   30-59: 7-12  (age30_34 .. age55_59)
#   60+:   13-18 (age60_64 .. age85_up)

compute_age_groups <- function(df, yr) {
  cols <- paste0(age_bands, "_", yr)
  A <- as.matrix(df[, cols])
  total <- rowSums(A)
  total_safe <- ifelse(total == 0, NA_real_, total)
  data.frame(
    tract     = as.numeric(gsub("[^0-9]", "", as.character(df$tract20l))),
    year      = yr,
    total_pop = total,
    pct_0_19  = rowSums(A[, 1:4])   / total_safe * 100,
    pct_20_29 = rowSums(A[, 5:6])   / total_safe * 100,
    pct_30_59 = rowSums(A[, 7:12])  / total_safe * 100,
    pct_60p   = rowSums(A[, 13:18]) / total_safe * 100
  )
}

age_groups <- bind_rows(
  compute_age_groups(age, 2010),
  compute_age_groups(age, 2020)
)

cat("Age groups: ", nrow(age_groups), "rows,",
    n_distinct(age_groups$tract), "tracts\n")
cat("Spot check (first tract, 2010):\n")
print(head(age_groups %>% filter(year == 2010), 3))

write.csv(age_groups,
          here::here("data", "age_pct_tract.csv"),
          row.names = FALSE)


# -------------------------------------------------------
# 2. RACE: 6 individual categories → percentages
# -------------------------------------------------------
race <- read.csv(here::here("data", "race", "race_10_20.csv"))

race_cats   <- c("AIAN", "ASIAN", "Black", "NHOPI", "Two", "White")
race_labels <- c("AIAN", "Asian", "Black", "NHOPI", "Two_or_More", "White")

compute_race_pct <- function(df, yr) {
  cols <- paste0(race_cats, yr)
  R <- as.matrix(df[, cols])
  total <- rowSums(R)
  total_safe <- ifelse(total == 0, NA_real_, total)
  P <- R / total_safe * 100
  out <- data.frame(
    tract     = as.numeric(gsub("[^0-9]", "", as.character(df$tract20l))),
    year      = yr,
    total_pop = total
  )
  for (i in seq_along(race_labels)) {
    out[[paste0("pct_", race_labels[i])]] <- P[, i]
  }
  out
}

race_pct <- bind_rows(
  compute_race_pct(race, 2010),
  compute_race_pct(race, 2020)
)

cat("\nRace pct: ", nrow(race_pct), "rows,",
    n_distinct(race_pct$tract), "tracts\n")
cat("Spot check (first tract, 2020):\n")
print(head(race_pct %>% filter(year == 2020), 3))

write.csv(race_pct,
          here::here("data", "race_pct_tract.csv"),
          row.names = FALSE)

cat("\nOutputs written:\n")
cat("  data/tract_age_groups.csv\n")
cat("  data/tract_race_pct.csv\n")
cat("DONE\n")
