library(here)
library(tidyverse)

here::i_am("code/2_2_tract_ndvi.R")

# ============================================================
# TRACT-LEVEL PSEUDO-RESIDUAL CHANGE ANALYSIS, 2010 -> 2020
# Greenness (NDVI) only, relative to housing density.
#
# Companion to 1_2_block_ndvi.R at the census tract scale.
# Revised scenario set based on 2_1_tract_scenario_comparison:
#   - Season is NOT collapsible (|Cohen's d| > 1.3;
#     cross-season ΔGHI correlation ~0).
#   - Buffer IS collapsible at the tract level (nested F-tests
#     non-significant, cross-buffer r = 0.75–0.99, d < 0.07).
#
# Final scenario set (2, down from 6):
#   - Dormant_Season x 1000 m
#   - Growing_Season x 1000 m
#
# Modeled after Lambert et al. (2025) Conservation Science and
# Practice: "Building the neighborhood for the trees"
#
# GUIDING QUESTIONS:
#   1. Where is housing being built and who lives there
#      (and how have demographics changed over time)?
#   2. Who is experiencing changes in exposure to green
#      environments?
#
# DEMOGRAPHIC CATEGORIES:
#   Race: AIAN, Asian, Black, NHOPI, Two or More, White
#         (individual categories, not collapsed)
#   Age:  0-19, 20-29, 30-59, 60+
# ============================================================

library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(ggplot2)

out_dir <- here::here("output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------
# 1. READ DATA
# -----------------------------
ndvi_tract <- read.csv(here::here("data", "tract_final.csv"))

# -----------------------------
# 2. BASIC CLEANING
# -----------------------------
df <- ndvi_tract %>%
  mutate(
    across(all_of(c("NDVI", "HU_density")), as.numeric),
    year   = as.integer(year),
    buffer = as.integer(buffer)
  ) %>%
  filter(year %in% c(2010, 2020)) %>%
  filter(!is.na(NDVI), !is.na(HU_density), HU_density > 0) %>%
  mutate(log_HU_density = log(HU_density))

# -----------------------------
# 3. PIVOT TO WIDE
# -----------------------------
wide <- df %>%
  select(tract, season, buffer, year, log_HU_density, NDVI, HU, POP, area_ha) %>%
  pivot_wider(names_from  = year,
              values_from = c(log_HU_density, NDVI, HU, POP),
              names_sep   = "_") %>%
  drop_na() %>%
  mutate(delta_log_HU = log_HU_density_2020 - log_HU_density_2010,
         delta_NDVI   = NDVI_2020 - NDVI_2010,
         delta_HU     = HU_2020 - HU_2010,
         delta_POP    = POP_2020 - POP_2010,
         HU_density_2010 = HU_2010 / area_ha,
         HU_density_2020 = HU_2020 / area_ha,
         delta_HU_density = HU_density_2020 - HU_density_2010)

# -----------------------------
# 4. RESTRICT TO UPDATED SCENARIO SET
# -----------------------------
scenarios <- tibble::tribble(
  ~season,           ~buffer,
  "Dormant_Season",     1000,
  "Growing_Season",     1000
)

wide <- wide %>% semi_join(scenarios, by = c("season", "buffer"))

cat("Retained rows:", nrow(wide), "\n")
cat("Unique tracts:", n_distinct(wide$tract), "\n\n")


# ============================================================
# 5. FIT 2010 BASELINE REGRESSIONS
# ============================================================
results <- scenarios %>%
  mutate(data = map2(season, buffer, ~ wide %>%
                       filter(season == .x, buffer == .y))) %>%
  mutate(mod_green = map(data, ~ lm(NDVI_2010 ~ log_HU_density_2010, data = .x)))

reg_table <- results %>%
  mutate(
    tidy_green = map(mod_green, ~ {
      s <- summary(.x)
      tibble(
        intercept = coef(.x)[1],
        slope     = coef(.x)[2],
        slope_se  = coef(s)[2, 2],
        slope_p   = coef(s)[2, 4],
        R2        = s$r.squared,
        n         = nobs(.x)
      )
    })
  ) %>%
  select(season, buffer, tidy_green) %>%
  unnest(tidy_green)

cat("\n================================================================\n")
cat("TABLE 1: 2010 OLS Regression Coefficients (Tract-Level Greenness)\n")
cat("  NDVI_2010 ~ log(HU_density_2010)\n")
cat("================================================================\n")
print(as.data.frame(reg_table), digits = 4, row.names = FALSE)

write.csv(reg_table,
          file.path(out_dir, "tract_table1_regression_coefficients.csv"),
          row.names = FALSE)


# ============================================================
# 6. COMPUTE RESIDUALS AND PSEUDO-RESIDUALS
# ============================================================
results <- results %>%
  mutate(
    data = map2(data, mod_green, function(d, mg) {
      d %>%
        mutate(
          resid_green_2010 = NDVI_2010 - predict(mg, newdata = d),
          pred_green_2020  = predict(mg, newdata = data.frame(
            log_HU_density_2010 = d$log_HU_density_2020)),
          resid_green_2020 = NDVI_2020 - pred_green_2020,
          change_green     = resid_green_2020 - resid_green_2010
        )
    })
  )

all_data <- results %>%
  select(data) %>%
  unnest(data) %>%
  mutate(
    state_green_2010 = ifelse(resid_green_2010 >= 0, "Bright", "Dim"),
    state_green_2020 = ifelse(resid_green_2020 >= 0, "Bright", "Dim"),
    transition_green = paste0(state_green_2010, " -> ", state_green_2020),
    direction_green  = case_when(
      change_green > 0 ~ "Brightened",
      change_green < 0 ~ "Dimmed",
      TRUE ~ "No change"
    ),
    scenario = factor(paste(ifelse(season == "Dormant_Season",
                                   "Dormant", "Growing"),
                            "1000 m"),
                      levels = c("Dormant 1000 m", "Growing 1000 m"))
  )


# ============================================================
# TABLE 2: Change summary by scenario
# ============================================================
change_summary <- all_data %>%
  group_by(season, buffer) %>%
  summarise(
    n                    = n(),
    mean_change_green    = mean(change_green),
    sd_change_green      = sd(change_green),
    median_change_green  = median(change_green),
    pct_brightened_green = mean(change_green > 0) * 100,
    pct_dimmed_green     = mean(change_green < 0) * 100,
    .groups = "drop"
  )

cat("\n================================================================\n")
cat("TABLE 2: Tract-Level Greenness Change Summary by Scenario\n")
cat("================================================================\n")
print(as.data.frame(change_summary), digits = 4, row.names = FALSE)

write.csv(change_summary,
          file.path(out_dir, "tract_table2_change_summary.csv"),
          row.names = FALSE)


# ============================================================
# TABLE 3: Bright/Dim transitions
# ============================================================
transition_green <- all_data %>%
  group_by(season, buffer, transition_green) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = transition_green, values_from = n, values_fill = 0)

cat("\n================================================================\n")
cat("TABLE 3: Greenness Bright/Dim Transitions (2010 -> 2020)\n")
cat("================================================================\n")
print(as.data.frame(transition_green), row.names = FALSE)

write.csv(transition_green,
          file.path(out_dir, "tract_table3_transitions.csv"),
          row.names = FALSE)


# ============================================================
# 7. STATISTICAL COMPARISONS
# ============================================================
cat("\n================================================================\n")
cat("PAIRED t-TEST: Dormant_Season_1000 vs Growing_Season_1000\n")
cat("================================================================\n")

paired <- all_data %>%
  select(tract, season, change_green) %>%
  pivot_wider(names_from = season, values_from = change_green) %>%
  drop_na()

print(t.test(paired$Dormant_Season, paired$Growing_Season, paired = TRUE))


# ============================================================
# 8. VISUALIZATIONS (GREENNESS)
# ============================================================
theme_lambert <- theme_minimal(base_size = 11) +
  theme(strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom")

# FIGURE 1: 2010 baseline
fig1 <- ggplot(all_data,
               aes(x = log_HU_density_2010, y = NDVI_2010)) +
  geom_point(alpha = 0.2, size = 0.5, color = "grey40") +
  geom_smooth(method = "lm", se = TRUE, color = "#2c7bb6", linewidth = 0.8) +
  facet_wrap(~ scenario) +
  labs(title    = "Tract-Level 2010 Baseline: Greenness vs. Housing Density",
       subtitle = "OLS regression of NDVI on log(HU density/ha), 2 retained scenarios",
       x = "log(Housing Unit Density per ha), 2010",
       y = "NDVI, 2010") +
  theme_lambert

ggsave(file.path(out_dir, "tract_fig1_baseline_green_regression.png"),
       fig1, width = 10, height = 5, dpi = 300)

# FIGURE 2: ΔGHI vs housing density change
fig2 <- ggplot(all_data,
               aes(x = delta_log_HU, y = change_green)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_point(aes(color = direction_green), alpha = 0.4, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7) +
  scale_color_manual(
    values = c("Brightened" = "#fdae61", "Dimmed" = "#878787", "No change" = "grey50"),
    name = "Direction") +
  facet_wrap(~ scenario) +
  labs(title    = "Tract-Level ΔGHI vs. Housing Density Change",
       subtitle = "Positive = brightening (greener than expected); Negative = dimming",
       x = expression(Delta ~ "log(HU density/ha), 2010" %->% "2020"),
       y = expression(Delta ~ "Greenness Index")) +
  theme_lambert

ggsave(file.path(out_dir, "tract_fig2_change_green_vs_housing.png"),
       fig2, width = 10, height = 5, dpi = 300)

# FIGURE 3: ΔGHI vs raw NDVI change
fig3 <- ggplot(all_data,
               aes(x = delta_NDVI, y = change_green)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_point(aes(color = direction_green), alpha = 0.4, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7) +
  scale_color_manual(
    values = c("Brightened" = "#fdae61", "Dimmed" = "#878787", "No change" = "grey50"),
    name = "Direction") +
  facet_wrap(~ scenario) +
  labs(title    = "Tract-Level ΔGHI vs. Raw NDVI Change",
       x = expression(Delta ~ "NDVI, 2010" %->% "2020"),
       y = expression(Delta ~ "Greenness Index")) +
  theme_lambert

ggsave(file.path(out_dir, "tract_fig3_change_green_vs_ndvi.png"),
       fig3, width = 10, height = 5, dpi = 300)

# FIGURE 4: Arrow plot of trajectories
set.seed(42)
arrow_sample <- all_data %>%
  filter(abs(delta_NDVI) > 0.005 | abs(delta_log_HU) > 0.02) %>%
  group_by(scenario) %>%
  slice_sample(n = min(400, n())) %>%
  ungroup()

fig4 <- ggplot(arrow_sample) +
  geom_segment(
    aes(x = log_HU_density_2010, xend = log_HU_density_2020,
        y = NDVI_2010, yend = NDVI_2020,
        color = direction_green),
    arrow = arrow(length = unit(0.07, "inches"), type = "closed"),
    alpha = 0.5, linewidth = 0.35) +
  scale_color_manual(
    values = c("Brightened" = "#fdae61", "Dimmed" = "#878787", "No change" = "grey50"),
    name = "Index Direction") +
  facet_wrap(~ scenario) +
  labs(title    = "Tract-Level Trajectories: Greenness vs. Housing Density, 2010-2020",
       subtitle = "Arrows show tract-level change (subsample); color = index brightening/dimming",
       x = "log(Housing Unit Density per ha)",
       y = "NDVI") +
  theme_lambert

ggsave(file.path(out_dir, "tract_fig4_arrow_green_trajectories.png"),
       fig4, width = 10, height = 5, dpi = 300)

# FIGURE 5: Boxplot
fig5 <- ggplot(all_data, aes(x = scenario, y = change_green, fill = scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5, width = 0.5) +
  scale_fill_manual(values = c("Dormant 1000 m" = "#abd9e9",
                                "Growing 1000 m" = "#66c164"),
                    guide = "none") +
  labs(title    = "Distribution of Tract-Level ΔGHI by Scenario",
       x = NULL,
       y = "ΔGHI (pseudo-residual change)") +
  theme_lambert

ggsave(file.path(out_dir, "tract_fig5_boxplot_change.png"),
       fig5, width = 7, height = 5, dpi = 300)

# FIGURE 6: Mean ΔGHI with 95% CI
ci_data <- all_data %>%
  group_by(scenario, season, buffer) %>%
  summarise(mean = mean(change_green),
            se   = sd(change_green) / sqrt(n()),
            .groups = "drop") %>%
  mutate(lo = mean - 1.96 * se, hi = mean + 1.96 * se)

fig6 <- ggplot(ci_data, aes(x = scenario, y = mean, color = scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.8) +
  scale_color_manual(values = c("Dormant 1000 m" = "#2c7bb6",
                                 "Growing 1000 m" = "#1a9641"),
                     guide = "none") +
  labs(title = "Tract-Level Mean ΔGHI with 95% CI",
       x = NULL, y = "Mean ΔGHI") +
  theme_lambert

ggsave(file.path(out_dir, "tract_fig6_mean_change_ci.png"),
       fig6, width = 7, height = 5, dpi = 300)


# ============================================================
# 9. DEMOGRAPHIC CHANGE: RACE & AGE
#
# Race: 6 individual categories (AIAN, Asian, Black, NHOPI,
#       Two or More, White) — NOT collapsed to White/non-White.
# Age:  4 groups (0-19, 20-29, 30-59, 60+).
# ============================================================

age  <- read.csv(here::here("data", "age", "age_10_20.csv"))
race <- read.csv(here::here("data", "race", "race_10_20.csv"))

# Normalize tract ID for joins
age$tract  <- as.numeric(gsub("[^0-9]", "", as.character(age$tract20l)))
race$tract <- as.numeric(gsub("[^0-9]", "", as.character(race$tract20l)))

# --- 9a. AGE: four groups (0-19, 20-29, 30-59, 60+) ---
age_bands <- c("age0_4","age5_9","age10_14","age15_19","age20_24","age25_29",
               "age30_34","age35_39","age40_44","age45_49","age50_54","age55_59",
               "age60_64","age65_69","age70_74","age75_79","age80_84","age85_up")

# Band indices for each group:
#   0-19:  bands 1-4  (age0_4 .. age15_19)
#   20-29: bands 5-6  (age20_24, age25_29)
#   30-59: bands 7-12 (age30_34 .. age55_59)
#   60+:   bands 13-18 (age60_64 .. age85_up)

age_group_pct <- function(df, yr) {
  cols <- paste0(age_bands, "_", yr)
  A <- as.matrix(df[, cols])
  total <- rowSums(A)
  total_safe <- ifelse(total == 0, NA_real_, total)
  data.frame(
    total      = total,
    pct_0_19   = rowSums(A[, 1:4])   / total_safe * 100,
    pct_20_29  = rowSums(A[, 5:6])   / total_safe * 100,
    pct_30_59  = rowSums(A[, 7:12])  / total_safe * 100,
    pct_60p    = rowSums(A[, 13:18]) / total_safe * 100
  )
}

a10 <- age_group_pct(age, 2010)
a20 <- age_group_pct(age, 2020)

# Delta-chi on the 4 age groups
age_delta_chi <- function(df) {
  A10 <- as.matrix(df[, paste0(age_bands, "_2010")])
  A20 <- as.matrix(df[, paste0(age_bands, "_2020")])
  # Aggregate to 4 groups before computing chi-square distance
  G10 <- cbind(rowSums(A10[, 1:4]),  rowSums(A10[, 5:6]),
               rowSums(A10[, 7:12]), rowSums(A10[, 13:18]))
  G20 <- cbind(rowSums(A20[, 1:4]),  rowSums(A20[, 5:6]),
               rowSums(A20[, 7:12]), rowSums(A20[, 13:18]))
  T10 <- rowSums(G10); T20 <- rowSums(G20)
  P10 <- G10 / ifelse(T10 == 0, NA_real_, T10)
  P20 <- G20 / ifelse(T20 == 0, NA_real_, T20)
  Pbar <- (P10 + P20) / 2
  num  <- (P20 - P10)^2
  ratio <- ifelse(Pbar > 0, num / Pbar, 0)
  rowSums(ratio, na.rm = TRUE)
}

age_out <- data.frame(
  tract            = age$tract,
  total_pop_2010   = a10$total,
  total_pop_2020   = a20$total,
  pct_0_19_2010    = a10$pct_0_19,
  pct_0_19_2020    = a20$pct_0_19,
  pct_20_29_2010   = a10$pct_20_29,
  pct_20_29_2020   = a20$pct_20_29,
  pct_30_59_2010   = a10$pct_30_59,
  pct_30_59_2020   = a20$pct_30_59,
  pct_60p_2010     = a10$pct_60p,
  pct_60p_2020     = a20$pct_60p,
  delta_chi_age    = age_delta_chi(age)
)
age_out$delta_pct_0_19  <- age_out$pct_0_19_2020  - age_out$pct_0_19_2010
age_out$delta_pct_20_29 <- age_out$pct_20_29_2020 - age_out$pct_20_29_2010
age_out$delta_pct_30_59 <- age_out$pct_30_59_2020 - age_out$pct_30_59_2010
age_out$delta_pct_60p   <- age_out$pct_60p_2020   - age_out$pct_60p_2010
age_out$delta_pop       <- age_out$total_pop_2020  - age_out$total_pop_2010

# --- 9b. RACE: 6 individual categories ---
race_cats <- c("AIAN", "ASIAN", "Black", "NHOPI", "Two", "White")
race_labels <- c("AIAN", "Asian", "Black", "NHOPI", "Two_or_More", "White")

R10 <- as.matrix(race[, paste0(race_cats, "2010")])
R20 <- as.matrix(race[, paste0(race_cats, "2020")])
T10 <- rowSums(R10); T20 <- rowSums(R20)
P10 <- R10 / ifelse(T10 == 0, NA_real_, T10)
P20 <- R20 / ifelse(T20 == 0, NA_real_, T20)
Pbar_race <- (P10 + P20) / 2
num_race  <- (P20 - P10)^2
ratio_race <- ifelse(Pbar_race > 0, num_race / Pbar_race, 0)

race_out <- data.frame(
  tract           = race$tract,
  total_race_2010 = T10,
  total_race_2020 = T20,
  delta_chi_race  = rowSums(ratio_race, na.rm = TRUE)
)

# Add pct columns for each race, both years
for (i in seq_along(race_cats)) {
  race_out[[paste0("pct_", race_labels[i], "_2010")]] <- P10[, paste0(race_cats[i], "2010")] * 100
  race_out[[paste0("pct_", race_labels[i], "_2020")]] <- P20[, paste0(race_cats[i], "2020")] * 100
}

# Compute deltas for each race
for (rl in race_labels) {
  race_out[[paste0("delta_pct_", rl)]] <-
    race_out[[paste0("pct_", rl, "_2020")]] - race_out[[paste0("pct_", rl, "_2010")]]
}

demo <- age_out %>%
  full_join(race_out, by = "tract")

# Merge demographics onto the analytical dataset
analysis_demo <- all_data %>%
  left_join(demo, by = "tract")

write.csv(demo,
          file.path(out_dir, "tract_demographics_with_delta_chi.csv"),
          row.names = FALSE)


# ============================================================
# 10. WHERE IS HOUSING BEING BUILT?
# ============================================================

# One row per tract (dormant scenario; HU/POP same across seasons)
tract_housing <- analysis_demo %>%
  filter(season == "Dormant_Season") %>%
  select(tract, area_ha, HU_2010, HU_2020, POP_2010, POP_2020,
         delta_HU, delta_POP, delta_HU_density,
         HU_density_2010, HU_density_2020,
         starts_with("pct_"), starts_with("delta_"),
         starts_with("total_pop"), starts_with("total_race"))

# --- TABLE 4: Housing growth summary ---
housing_summary <- tract_housing %>%
  summarise(
    n_tracts             = n(),
    total_HU_2010        = sum(HU_2010, na.rm = TRUE),
    total_HU_2020        = sum(HU_2020, na.rm = TRUE),
    total_delta_HU       = sum(delta_HU, na.rm = TRUE),
    pct_tracts_gained_HU = mean(delta_HU > 0, na.rm = TRUE) * 100,
    pct_tracts_lost_HU   = mean(delta_HU < 0, na.rm = TRUE) * 100,
    median_delta_HU      = median(delta_HU, na.rm = TRUE),
    mean_delta_HU        = mean(delta_HU, na.rm = TRUE),
    q25_delta_HU         = quantile(delta_HU, 0.25, na.rm = TRUE),
    q75_delta_HU         = quantile(delta_HU, 0.75, na.rm = TRUE),
    total_POP_2010       = sum(POP_2010, na.rm = TRUE),
    total_POP_2020       = sum(POP_2020, na.rm = TRUE),
    total_delta_POP      = sum(delta_POP, na.rm = TRUE)
  )

cat("\n================================================================\n")
cat("TABLE 4: Housing Growth Summary (All Tracts)\n")
cat("================================================================\n")
print(as.data.frame(housing_summary), digits = 4, row.names = FALSE)

write.csv(housing_summary,
          file.path(out_dir, "tract_table4_housing_summary.csv"),
          row.names = FALSE)

# --- Classify tracts by housing growth quartile ---
tract_housing <- tract_housing %>%
  mutate(
    hu_growth_q = cut(delta_HU,
                      breaks = quantile(delta_HU,
                                        probs = c(0, 0.25, 0.5, 0.75, 1),
                                        na.rm = TRUE),
                      labels = c("Q1 (least growth)",
                                 "Q2", "Q3",
                                 "Q4 (most growth)"),
                      include.lowest = TRUE)
  )

# --- TABLE 5: Demographics by housing growth quartile ---
demo_by_hu_growth <- tract_housing %>%
  group_by(hu_growth_q) %>%
  summarise(
    n_tracts = n(),
    mean_delta_HU   = mean(delta_HU, na.rm = TRUE),
    median_delta_HU = median(delta_HU, na.rm = TRUE),
    mean_delta_pop  = mean(delta_pop, na.rm = TRUE),
    # Race 2020
    across(all_of(paste0("pct_", race_labels, "_2020")),
           ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
    # Race 2010
    across(all_of(paste0("pct_", race_labels, "_2010")),
           ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
    # Race change
    across(all_of(paste0("delta_pct_", race_labels)),
           ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
    # Age 2020
    mean_pct_0_19_2020  = mean(pct_0_19_2020, na.rm = TRUE),
    mean_pct_20_29_2020 = mean(pct_20_29_2020, na.rm = TRUE),
    mean_pct_30_59_2020 = mean(pct_30_59_2020, na.rm = TRUE),
    mean_pct_60p_2020   = mean(pct_60p_2020, na.rm = TRUE),
    # Age 2010
    mean_pct_0_19_2010  = mean(pct_0_19_2010, na.rm = TRUE),
    mean_pct_20_29_2010 = mean(pct_20_29_2010, na.rm = TRUE),
    mean_pct_30_59_2010 = mean(pct_30_59_2010, na.rm = TRUE),
    mean_pct_60p_2010   = mean(pct_60p_2010, na.rm = TRUE),
    # Age change
    mean_delta_pct_0_19  = mean(delta_pct_0_19, na.rm = TRUE),
    mean_delta_pct_20_29 = mean(delta_pct_20_29, na.rm = TRUE),
    mean_delta_pct_30_59 = mean(delta_pct_30_59, na.rm = TRUE),
    mean_delta_pct_60p   = mean(delta_pct_60p, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n================================================================\n")
cat("TABLE 5: Demographics by Housing Growth Quartile\n")
cat("================================================================\n")
print(as.data.frame(demo_by_hu_growth), digits = 3, row.names = FALSE)

write.csv(demo_by_hu_growth,
          file.path(out_dir, "tract_table5_demo_by_hu_growth.csv"),
          row.names = FALSE)

# --- TABLE 6: Correlations — housing growth × demographics ---
race_pct_2010_cols <- paste0("pct_", race_labels, "_2010")
race_pct_2020_cols <- paste0("pct_", race_labels, "_2020")
race_delta_cols    <- paste0("delta_pct_", race_labels)

hu_demo_cors <- tract_housing %>%
  summarise(
    n = n(),
    # Race 2010
    across(all_of(race_pct_2010_cols),
           ~ cor(delta_HU, .x, use = "complete.obs"),
           .names = "cor_dHU_{.col}"),
    # Race 2020
    across(all_of(race_pct_2020_cols),
           ~ cor(delta_HU, .x, use = "complete.obs"),
           .names = "cor_dHU_{.col}"),
    # Race change
    across(all_of(race_delta_cols),
           ~ cor(delta_HU, .x, use = "complete.obs"),
           .names = "cor_dHU_{.col}"),
    # Age 2010
    cor_dHU_pct_0_19_2010  = cor(delta_HU, pct_0_19_2010, use = "complete.obs"),
    cor_dHU_pct_20_29_2010 = cor(delta_HU, pct_20_29_2010, use = "complete.obs"),
    cor_dHU_pct_30_59_2010 = cor(delta_HU, pct_30_59_2010, use = "complete.obs"),
    cor_dHU_pct_60p_2010   = cor(delta_HU, pct_60p_2010, use = "complete.obs"),
    # Age 2020
    cor_dHU_pct_0_19_2020  = cor(delta_HU, pct_0_19_2020, use = "complete.obs"),
    cor_dHU_pct_20_29_2020 = cor(delta_HU, pct_20_29_2020, use = "complete.obs"),
    cor_dHU_pct_30_59_2020 = cor(delta_HU, pct_30_59_2020, use = "complete.obs"),
    cor_dHU_pct_60p_2020   = cor(delta_HU, pct_60p_2020, use = "complete.obs"),
    # Age change
    cor_dHU_delta_pct_0_19  = cor(delta_HU, delta_pct_0_19, use = "complete.obs"),
    cor_dHU_delta_pct_20_29 = cor(delta_HU, delta_pct_20_29, use = "complete.obs"),
    cor_dHU_delta_pct_30_59 = cor(delta_HU, delta_pct_30_59, use = "complete.obs"),
    cor_dHU_delta_pct_60p   = cor(delta_HU, delta_pct_60p, use = "complete.obs"),
    # Chi-square distances
    cor_dHU_dchiAge  = cor(delta_HU, delta_chi_age, use = "complete.obs"),
    cor_dHU_dchiRace = cor(delta_HU, delta_chi_race, use = "complete.obs")
  )

cat("\n================================================================\n")
cat("TABLE 6: Correlations — Housing Growth × Demographics\n")
cat("================================================================\n")
print(as.data.frame(t(hu_demo_cors)), digits = 3)

write.csv(hu_demo_cors,
          file.path(out_dir, "tract_table6_hu_demo_correlations.csv"),
          row.names = FALSE)


# ============================================================
# 11. WHO IS EXPERIENCING CHANGES IN GREEN EXPOSURE?
#
# Stratify ΔGHI by:
#   - Each race category's share (terciles of % for each race)
#   - Age group composition (terciles of % 60+)
# ============================================================

# --- Classify tracts by % White tercile (2020) ---
# White is the largest group; terciles of % White are the most
# informative single-race stratifier for WA.
analysis_demo <- analysis_demo %>%
  mutate(
    white_tercile_2020 = cut(pct_White_2020,
                             breaks = quantile(pct_White_2020,
                                               probs = c(0, 1/3, 2/3, 1),
                                               na.rm = TRUE),
                             labels = c("Low White",
                                        "Mid White",
                                        "High White"),
                             include.lowest = TRUE),
    pct_60p_tercile_2020 = cut(pct_60p_2020,
                               breaks = quantile(pct_60p_2020,
                                                 probs = c(0, 1/3, 2/3, 1),
                                                 na.rm = TRUE),
                               labels = c("Low 60+",
                                          "Mid 60+",
                                          "High 60+"),
                               include.lowest = TRUE)
  )

# --- TABLE 7: ΔGHI by % White tercile (2020) ---
dghi_by_white <- analysis_demo %>%
  filter(!is.na(white_tercile_2020)) %>%
  group_by(scenario, white_tercile_2020) %>%
  summarise(
    n              = n(),
    mean_DGHI      = mean(change_green, na.rm = TRUE),
    sd_DGHI        = sd(change_green, na.rm = TRUE),
    median_DGHI    = median(change_green, na.rm = TRUE),
    pct_brightened  = mean(change_green > 0, na.rm = TRUE) * 100,
    mean_pct_White = mean(pct_White_2020, na.rm = TRUE),
    mean_NDVI_2010 = mean(NDVI_2010, na.rm = TRUE),
    mean_NDVI_2020 = mean(NDVI_2020, na.rm = TRUE),
    mean_delta_NDVI = mean(delta_NDVI, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n================================================================\n")
cat("TABLE 7: ΔGHI by % White Tercile (2020)\n")
cat("================================================================\n")
print(as.data.frame(dghi_by_white), digits = 3, row.names = FALSE)

write.csv(dghi_by_white,
          file.path(out_dir, "tract_table7_dghi_by_white_tercile.csv"),
          row.names = FALSE)

# --- TABLE 8: ΔGHI by % 60+ tercile (2020) ---
dghi_by_60p <- analysis_demo %>%
  filter(!is.na(pct_60p_tercile_2020)) %>%
  group_by(scenario, pct_60p_tercile_2020) %>%
  summarise(
    n              = n(),
    mean_DGHI      = mean(change_green, na.rm = TRUE),
    sd_DGHI        = sd(change_green, na.rm = TRUE),
    median_DGHI    = median(change_green, na.rm = TRUE),
    pct_brightened  = mean(change_green > 0, na.rm = TRUE) * 100,
    mean_pct_60p   = mean(pct_60p_2020, na.rm = TRUE),
    mean_NDVI_2010 = mean(NDVI_2010, na.rm = TRUE),
    mean_NDVI_2020 = mean(NDVI_2020, na.rm = TRUE),
    mean_delta_NDVI = mean(delta_NDVI, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n================================================================\n")
cat("TABLE 8: ΔGHI by % 60+ Tercile (2020)\n")
cat("================================================================\n")
print(as.data.frame(dghi_by_60p), digits = 3, row.names = FALSE)

write.csv(dghi_by_60p,
          file.path(out_dir, "tract_table8_dghi_by_60p_tercile.csv"),
          row.names = FALSE)

# --- TABLE 9: Bright/Dim transitions by % White tercile ---
transitions_by_white <- analysis_demo %>%
  filter(!is.na(white_tercile_2020)) %>%
  group_by(scenario, white_tercile_2020, transition_green) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(scenario, white_tercile_2020) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

cat("\n================================================================\n")
cat("TABLE 9: Bright/Dim Transitions by % White Tercile (2020)\n")
cat("================================================================\n")
print(as.data.frame(transitions_by_white %>%
                      select(scenario, white_tercile_2020, transition_green, n, pct) %>%
                      arrange(scenario, white_tercile_2020, transition_green)),
      digits = 2, row.names = FALSE)

write.csv(transitions_by_white,
          file.path(out_dir, "tract_table9_transitions_by_white.csv"),
          row.names = FALSE)

# --- TABLE 10: NDVI levels & housing by % White tercile ---
ndvi_by_white <- analysis_demo %>%
  filter(!is.na(white_tercile_2020)) %>%
  group_by(scenario, white_tercile_2020) %>%
  summarise(
    n               = n(),
    mean_NDVI_2010  = mean(NDVI_2010, na.rm = TRUE),
    mean_NDVI_2020  = mean(NDVI_2020, na.rm = TRUE),
    mean_delta_NDVI = mean(delta_NDVI, na.rm = TRUE),
    mean_HU_dens_2010 = mean(HU_density_2010, na.rm = TRUE),
    mean_HU_dens_2020 = mean(HU_density_2020, na.rm = TRUE),
    mean_delta_HU   = mean(delta_HU, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n================================================================\n")
cat("TABLE 10: NDVI Levels & Housing by % White Tercile\n")
cat("================================================================\n")
print(as.data.frame(ndvi_by_white), digits = 3, row.names = FALSE)

write.csv(ndvi_by_white,
          file.path(out_dir, "tract_table10_ndvi_by_white.csv"),
          row.names = FALSE)


# --- TABLE 11: ΔGHI correlations with all demographic variables ---
cor_by_scenario <- analysis_demo %>%
  group_by(scenario) %>%
  summarise(
    n = sum(!is.na(change_green)),
    # Race level 2020
    across(all_of(paste0("pct_", race_labels, "_2020")),
           ~ cor(change_green, .x, use = "complete.obs"),
           .names = "cor_dGHI_{.col}"),
    # Race change
    across(all_of(paste0("delta_pct_", race_labels)),
           ~ cor(change_green, .x, use = "complete.obs"),
           .names = "cor_dGHI_{.col}"),
    # Age level 2020
    cor_dGHI_pct_0_19_2020  = cor(change_green, pct_0_19_2020, use = "complete.obs"),
    cor_dGHI_pct_20_29_2020 = cor(change_green, pct_20_29_2020, use = "complete.obs"),
    cor_dGHI_pct_30_59_2020 = cor(change_green, pct_30_59_2020, use = "complete.obs"),
    cor_dGHI_pct_60p_2020   = cor(change_green, pct_60p_2020, use = "complete.obs"),
    # Age change
    cor_dGHI_delta_pct_0_19  = cor(change_green, delta_pct_0_19, use = "complete.obs"),
    cor_dGHI_delta_pct_20_29 = cor(change_green, delta_pct_20_29, use = "complete.obs"),
    cor_dGHI_delta_pct_30_59 = cor(change_green, delta_pct_30_59, use = "complete.obs"),
    cor_dGHI_delta_pct_60p   = cor(change_green, delta_pct_60p, use = "complete.obs"),
    # Chi-square distances
    cor_dGHI_dchi_age  = cor(change_green, delta_chi_age, use = "complete.obs"),
    cor_dGHI_dchi_race = cor(change_green, delta_chi_race, use = "complete.obs"),
    # Housing growth
    cor_dGHI_dHU  = cor(change_green, delta_HU, use = "complete.obs"),
    cor_dGHI_dPOP = cor(change_green, delta_POP, use = "complete.obs"),
    .groups = "drop"
  )

cat("\n================================================================\n")
cat("TABLE 11: ΔGHI Correlations — Demographics & Housing Growth\n")
cat("================================================================\n")
print(as.data.frame(cor_by_scenario), digits = 3, row.names = FALSE)

write.csv(cor_by_scenario,
          file.path(out_dir, "tract_table11_dghi_correlations.csv"),
          row.names = FALSE)

# --- TABLE 12: Demographic summary (all tracts with demo data) ---
demo_summary <- demo %>%
  summarise(
    n_tracts = n(),
    # Age groups
    mean_pct_0_19_2010  = mean(pct_0_19_2010, na.rm = TRUE),
    mean_pct_0_19_2020  = mean(pct_0_19_2020, na.rm = TRUE),
    mean_pct_20_29_2010 = mean(pct_20_29_2010, na.rm = TRUE),
    mean_pct_20_29_2020 = mean(pct_20_29_2020, na.rm = TRUE),
    mean_pct_30_59_2010 = mean(pct_30_59_2010, na.rm = TRUE),
    mean_pct_30_59_2020 = mean(pct_30_59_2020, na.rm = TRUE),
    mean_pct_60p_2010   = mean(pct_60p_2010, na.rm = TRUE),
    mean_pct_60p_2020   = mean(pct_60p_2020, na.rm = TRUE),
    median_delta_chi_age  = median(delta_chi_age, na.rm = TRUE),
    median_delta_chi_race = median(delta_chi_race, na.rm = TRUE),
    # Race
    across(all_of(paste0("pct_", race_labels, "_2010")),
           ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}"),
    across(all_of(paste0("pct_", race_labels, "_2020")),
           ~ mean(.x, na.rm = TRUE), .names = "mean_{.col}")
  )

cat("\n================================================================\n")
cat("TABLE 12: Demographic Summary (All WA Tracts)\n")
cat("================================================================\n")
print(as.data.frame(demo_summary), digits = 4, row.names = FALSE)

write.csv(demo_summary,
          file.path(out_dir, "tract_table12_demographic_summary.csv"),
          row.names = FALSE)


# ============================================================
# 12. DEMOGRAPHIC FIGURES
# ============================================================

# FIGURE 7: ΔGHI by % White tercile (boxplot)
fig7 <- analysis_demo %>%
  filter(!is.na(white_tercile_2020)) %>%
  ggplot(aes(x = white_tercile_2020, y = change_green,
             fill = white_tercile_2020)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.4, width = 0.6) +
  scale_fill_manual(values = c("Low White"  = "#e34a33",
                                "Mid White"  = "#fdbb84",
                                "High White" = "#fee8c8"),
                    guide = "none") +
  facet_wrap(~ scenario) +
  labs(title    = "Greenness Change by Tract % White (2020)",
       subtitle = "Who is experiencing dimming vs. brightening?",
       x = "Tract % White tercile (2020)",
       y = "ΔGHI (pseudo-residual change)") +
  theme_lambert

ggsave(file.path(out_dir, "tract_fig7_dghi_by_white.png"),
       fig7, width = 10, height = 5, dpi = 300)

# FIGURE 8: ΔGHI by % 60+ tercile (boxplot)
fig8 <- analysis_demo %>%
  filter(!is.na(pct_60p_tercile_2020)) %>%
  ggplot(aes(x = pct_60p_tercile_2020, y = change_green,
             fill = pct_60p_tercile_2020)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.4, width = 0.6) +
  scale_fill_manual(values = c("Low 60+"  = "#deebf7",
                                "Mid 60+"  = "#9ecae1",
                                "High 60+" = "#3182bd"),
                    guide = "none") +
  facet_wrap(~ scenario) +
  labs(title    = "Greenness Change by Tract % Aged 60+ (2020)",
       subtitle = "Who is experiencing dimming vs. brightening?",
       x = "Tract % aged 60+ tercile (2020)",
       y = "ΔGHI (pseudo-residual change)") +
  theme_lambert

ggsave(file.path(out_dir, "tract_fig8_dghi_by_60p.png"),
       fig8, width = 10, height = 5, dpi = 300)

# FIGURE 9: Housing growth vs. ΔGHI, colored by % White tercile
fig9 <- analysis_demo %>%
  filter(!is.na(white_tercile_2020)) %>%
  ggplot(aes(x = delta_HU, y = change_green,
             color = white_tercile_2020)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.4, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  scale_color_manual(values = c("Low White"  = "#e34a33",
                                 "Mid White"  = "#fdbb84",
                                 "High White" = "#fee8c8"),
                     name = "% White tercile (2020)") +
  facet_wrap(~ scenario) +
  labs(title    = "Housing Growth vs. Greenness Change, by % White",
       subtitle = "Do low-White tracts that gain housing also lose greenness?",
       x = "ΔHU (housing units gained, 2010-2020)",
       y = "ΔGHI") +
  theme_lambert

ggsave(file.path(out_dir, "tract_fig9_hu_growth_dghi_by_white.png"),
       fig9, width = 10, height = 5, dpi = 300)

# FIGURE 10: Baseline NDVI (2010) by % White tercile
fig10 <- analysis_demo %>%
  filter(!is.na(white_tercile_2020)) %>%
  ggplot(aes(x = white_tercile_2020, y = NDVI_2010,
             fill = white_tercile_2020)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.4, width = 0.6) +
  scale_fill_manual(values = c("Low White"  = "#e34a33",
                                "Mid White"  = "#fdbb84",
                                "High White" = "#fee8c8"),
                    guide = "none") +
  facet_wrap(~ scenario) +
  labs(title    = "Baseline Greenness (2010) by Tract % White",
       subtitle = "Do low-White tracts start with less greenness?",
       x = "Tract % White tercile (2020)",
       y = "NDVI, 2010") +
  theme_lambert

ggsave(file.path(out_dir, "tract_fig10_baseline_ndvi_by_white.png"),
       fig10, width = 10, height = 5, dpi = 300)

# FIGURE 11: Race composition of housing growth quartiles (stacked bar)
# Build long-format race composition for 2010 and 2020 by HU growth quartile
race_by_hu_long <- demo_by_hu_growth %>%
  select(hu_growth_q,
         starts_with("mean_pct_") & matches("(AIAN|Asian|Black|NHOPI|Two_or_More|White)")) %>%
  pivot_longer(-hu_growth_q, names_to = "variable", values_to = "value") %>%
  mutate(
    year = ifelse(grepl("2010", variable), "2010", "2020"),
    race = gsub("mean_pct_(.+)_(2010|2020)", "\\1", variable),
    race = gsub("_", " ", race)
  )

fig11 <- ggplot(race_by_hu_long,
                aes(x = hu_growth_q, y = value, fill = race)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_brewer(palette = "Set2", name = "Race") +
  facet_wrap(~ year) +
  labs(title    = "Race Composition of Housing Growth Quartiles",
       subtitle = "Who lives in tracts that gained the most housing?",
       x = "Housing growth quartile (ΔHU, 2010-2020)",
       y = "Mean % of tract population") +
  theme_lambert +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

ggsave(file.path(out_dir, "tract_fig11_race_by_hu_growth.png"),
       fig11, width = 10, height = 5, dpi = 300)

# FIGURE 12: Age composition of housing growth quartiles (stacked bar)
age_by_hu_long <- demo_by_hu_growth %>%
  select(hu_growth_q,
         matches("mean_pct_(0_19|20_29|30_59|60p)_(2010|2020)")) %>%
  pivot_longer(-hu_growth_q, names_to = "variable", values_to = "value") %>%
  mutate(
    year = ifelse(grepl("2010", variable), "2010", "2020"),
    age_group = case_when(
      grepl("0_19", variable)  ~ "0-19",
      grepl("20_29", variable) ~ "20-29",
      grepl("30_59", variable) ~ "30-59",
      grepl("60p", variable)   ~ "60+"
    ),
    age_group = factor(age_group, levels = c("0-19", "20-29", "30-59", "60+"))
  )

fig12 <- ggplot(age_by_hu_long,
                aes(x = hu_growth_q, y = value, fill = age_group)) +
  geom_col(position = "stack", width = 0.7) +
  scale_fill_brewer(palette = "Blues", name = "Age group") +
  facet_wrap(~ year) +
  labs(title    = "Age Composition of Housing Growth Quartiles",
       subtitle = "Who lives in tracts that gained the most housing?",
       x = "Housing growth quartile (ΔHU, 2010-2020)",
       y = "Mean % of tract population") +
  theme_lambert +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

ggsave(file.path(out_dir, "tract_fig12_age_by_hu_growth.png"),
       fig12, width = 10, height = 5, dpi = 300)

# FIGURE 13: Δ race % vs Δχ²_race, faceted by scenario × race
race_delta_long <- analysis_demo %>%
  select(tract, scenario, change_green, delta_chi_race,
         all_of(paste0("delta_pct_", race_labels))) %>%
  pivot_longer(cols = starts_with("delta_pct_"),
               names_to = "y_var", values_to = "y_value") %>%
  mutate(y_var = gsub("delta_pct_", "Δ % ", y_var),
         y_var = gsub("_", " ", y_var))

fig13 <- ggplot(race_delta_long,
               aes(x = delta_chi_race, y = y_value)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = change_green), alpha = 0.5, size = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
  scale_color_gradient2(low = "#878787", mid = "grey90", high = "#fdae61",
                        midpoint = 0, name = "ΔGHI") +
  facet_grid(y_var ~ scenario, scales = "free_y") +
  labs(title    = "Race Composition Change vs. Δχ² (Race)",
       x = expression(Delta ~ chi[Race]^2),
       y = "Δ (race %)") +
  theme_lambert +
  theme(strip.text.y = element_text(size = 8))

ggsave(file.path(out_dir, "tract_fig13_race_vs_delta_chi.png"),
       fig13, width = 10, height = 12, dpi = 300)

# FIGURE 14: Δ age group % vs Δχ²_age, faceted by scenario × age group
age_delta_long <- analysis_demo %>%
  select(tract, scenario, change_green, delta_chi_age,
         delta_pct_0_19, delta_pct_20_29, delta_pct_30_59, delta_pct_60p) %>%
  pivot_longer(cols = starts_with("delta_pct_"),
               names_to = "y_var", values_to = "y_value") %>%
  mutate(y_var = recode(y_var,
                        delta_pct_0_19  = "Δ % 0-19",
                        delta_pct_20_29 = "Δ % 20-29",
                        delta_pct_30_59 = "Δ % 30-59",
                        delta_pct_60p   = "Δ % 60+"))

fig14 <- ggplot(age_delta_long,
               aes(x = delta_chi_age, y = y_value)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(color = change_green), alpha = 0.5, size = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.5) +
  scale_color_gradient2(low = "#878787", mid = "grey90", high = "#fdae61",
                        midpoint = 0, name = "ΔGHI") +
  facet_grid(y_var ~ scenario, scales = "free_y") +
  labs(title    = "Age Composition Change vs. Δχ² (Age)",
       x = expression(Delta ~ chi[Age]^2),
       y = "Δ (age group %)") +
  theme_lambert

ggsave(file.path(out_dir, "tract_fig14_age_vs_delta_chi.png"),
       fig14, width = 10, height = 8, dpi = 300)

# FIGURE 15: ΔGHI vs Δχ² (greenness change vs demographic shift)
chi_long <- analysis_demo %>%
  select(tract, scenario, change_green, delta_chi_age, delta_chi_race) %>%
  pivot_longer(cols = c(delta_chi_age, delta_chi_race),
               names_to = "chi_var", values_to = "chi_value") %>%
  mutate(chi_var = recode(chi_var,
                          delta_chi_age  = "Δχ² (Age)",
                          delta_chi_race = "Δχ² (Race)"))

fig15 <- ggplot(chi_long,
               aes(x = chi_value, y = change_green)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(alpha = 0.4, size = 0.7, color = "#2c7bb6") +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
  facet_grid(scenario ~ chi_var, scales = "free_x") +
  labs(title    = "ΔGHI vs. Demographic Δχ²",
       subtitle = "Does greenness change track demographic reshuffling?",
       x = "Δχ²", y = expression(Delta ~ "GHI")) +
  theme_lambert

ggsave(file.path(out_dir, "tract_fig15_dghi_vs_delta_chi.png"),
       fig15, width = 10, height = 6, dpi = 300)


# ============================================================
# SAVE FULL ANALYTICAL DATASET
# ============================================================
write.csv(analysis_demo,
          file.path(out_dir, "tract_pseudo_residual_analysis.csv"),
          row.names = FALSE)

cat("\n================================================================\n")
cat("All outputs saved to:", out_dir, "\n")
cat("  Tables: tract_table1-tract_table12 CSVs\n")
cat("  Figures: tract_fig1-tract_fig15 PNGs\n")
cat("  Demographics: tract_demographics_with_delta_chi.csv\n")
cat("  Data: tract_pseudo_residual_analysis.csv\n")
cat("================================================================\n")
cat("ANALYSIS COMPLETE\n")
