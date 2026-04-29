library(here)
library(tidyverse)
library(ggplot2)

here::i_am("code/1_2_block_ndvi.R")

# ============================================================
# PSEUDO-RESIDUAL CHANGE ANALYSIS, 2010 -> 2020
# Greenness (NDVI) only, relative to housing density.
#
# Revised scenario set based on block_ndvi_scenario_comparison;
# Final scenario set (4, down from 6 in the original):
#   - Dormant_Season x 1000 m     (representative dormant buffer)
#   - Growing_Season x 500 m
#   - Growing_Season x 1000 m
#   - Growing_Season x 2000 m
# ============================================================


out_dir <- here::here("output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------
# 1. READ DATA
# -----------------------------
ndvi_block2 <- read.csv(here::here("data", "block_final.csv"))

# -----------------------------
# 2. BASIC CLEANING
# -----------------------------
df <- ndvi_block2 %>%
  mutate(
    across(all_of(c("NDVI", "HU_density")), as.numeric),
    year   = as.integer(year),
    buffer = as.integer(buffer)
  ) %>%
  filter(year %in% c(2010, 2020)) %>%
  filter(!is.na(NDVI), !is.na(HU_density)) %>%
  mutate(log_HU_density = log(pmax(HU_density, 1e-10)))

# -----------------------------
# 3. PIVOT TO WIDE
# -----------------------------
wide <- df %>%
  select(GEOID20, season, buffer, year, log_HU_density, NDVI) %>%
  pivot_wider(names_from  = year,
              values_from = c(log_HU_density, NDVI),
              names_sep   = "_") %>%
  drop_na() %>%
  mutate(delta_log_HU = log_HU_density_2020 - log_HU_density_2010,
         delta_NDVI   = NDVI_2020 - NDVI_2010)

# -----------------------------
# 4. RESTRICT TO UPDATED SCENARIO SET
# -----------------------------
scenarios <- tibble::tribble(
  ~season,           ~buffer,
  "Dormant_Season",     1000,
  "Growing_Season",      500,
  "Growing_Season",     1000,
  "Growing_Season",     2000
)

wide <- wide %>% semi_join(scenarios, by = c("season", "buffer"))

cat("Retained rows:", nrow(wide), "\n")
cat("Unique blocks:", n_distinct(wide$GEOID20), "\n\n")

# ============================================================
# 5. FIT 2010 BASELINE REGRESSIONS (one per scenario)
# ============================================================
results <- scenarios %>%
  mutate(data = map2(season, buffer, ~ wide %>%
                       filter(season == .x, buffer == .y))) %>%
  mutate(mod_green = map(data, ~ lm(NDVI_2010 ~ log_HU_density_2010, data = .x)))

# -------------------------------------------------------
# TABLE 1: 2010 Regression coefficients
# -------------------------------------------------------
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
cat("TABLE 1: 2010 OLS Regression Coefficients (Greenness)\n")
cat("  NDVI_2010 ~ log(HU_density_2010)\n")
cat("================================================================\n")
print(as.data.frame(reg_table), digits = 4, row.names = FALSE)

# write.csv(reg_table,
#           file.path(out_dir, "table1_regression_coefficients.csv"),
#           row.names = FALSE)


# ============================================================
# 6. COMPUTE RESIDUALS (2010) AND PSEUDO-RESIDUALS (2020)
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
    )
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
cat("TABLE 2: Greenness Change Summary by Scenario\n")
cat("================================================================\n")
print(as.data.frame(change_summary), digits = 4, row.names = FALSE)

# write.csv(change_summary,
#           file.path(out_dir, "table2_change_summary.csv"),
#           row.names = FALSE)


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

# write.csv(transition_green,
#           file.path(out_dir, "table3_transitions.csv"),
#           row.names = FALSE)


# ============================================================
# 7. STATISTICAL COMPARISONS
# ============================================================

# --- 7a. Paired t-test: Dormant (1000m) vs Growing (1000m) ---
cat("\n================================================================\n")
cat("PAIRED t-TEST: Dormant_Season_1000 vs Growing_Season_1000\n")
cat("  (tests the only within-block cross-season comparison retained)\n")
cat("================================================================\n")

paired <- all_data %>%
  filter(buffer == 1000) %>%
  select(GEOID20, season, change_green) %>%
  pivot_wider(names_from = season, values_from = change_green) %>%
  drop_na()

print(t.test(paired$Dormant_Season, paired$Growing_Season, paired = TRUE))

# --- 7b. ANOVA: buffer effect within growing season ---
cat("\n================================================================\n")
cat("ANOVA: change_green ~ buffer (Growing Season only)\n")
cat("================================================================\n")

grow <- all_data %>% filter(season == "Growing_Season")
grow$buffer_f <- factor(grow$buffer)
aov_grow <- aov(change_green ~ buffer_f, data = grow)
print(summary(aov_grow))

cat("\n--- Tukey HSD pairwise buffer comparisons (Growing Season) ---\n")
print(TukeyHSD(aov_grow))

# --- 7c. Cross-buffer correlation within growing season ---
cat("\n================================================================\n")
cat("Cross-buffer ΔGHI correlation (Growing Season)\n")
cat("================================================================\n")

cross_buf <- grow %>%
  select(GEOID20, buffer, change_green) %>%
  pivot_wider(names_from = buffer, values_from = change_green,
              names_prefix = "b") %>%
  drop_na()

print(round(cor(cross_buf %>% select(-GEOID20)), 4))


# ============================================================
# 8. VISUALIZATIONS
# ============================================================
theme_lambert <- theme_minimal(base_size = 11) +
  theme(strip.text       = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        legend.position  = "bottom")

# Scenario label for cleaner faceting (1 dormant + 3 growing)
all_data <- all_data %>%
  mutate(scenario = paste(
    ifelse(season == "Dormant_Season", "Dormant", "Growing"),
    paste0(buffer, " m")),
    scenario = factor(scenario,
                      levels = c("Dormant 1000 m",
                                 "Growing 500 m",
                                 "Growing 1000 m",
                                 "Growing 2000 m")))

# -------------------------------------------------------
# FIGURE 1: 2010 baseline — NDVI vs log(HU_density)
# -------------------------------------------------------
fig1 <- ggplot(all_data,
               aes(x = log_HU_density_2010, y = NDVI_2010)) +
  geom_point(alpha = 0.05, size = 0.3, color = "grey40") +
  geom_smooth(method = "lm", se = TRUE, color = "#2c7bb6", linewidth = 0.8) +
  facet_wrap(~ scenario, nrow = 1) +
  labs(title    = "2010 Baseline: Greenness vs. Housing Density",
       subtitle = "OLS regression of NDVI on log(HU density/ha), 4 retained scenarios",
       x = "log(Housing Unit Density per ha), 2010",
       y = "NDVI, 2010") +
  theme_lambert

# ggsave(file.path(out_dir, "fig1_baseline_green_regression.png"),
#        fig1, width = 12, height = 4, dpi = 300)

# -------------------------------------------------------
# FIGURE 2: ΔGHI vs housing density change
# -------------------------------------------------------
fig2 <- ggplot(all_data,
               aes(x = delta_log_HU, y = change_green)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_point(aes(color = direction_green), alpha = 0.15, size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7) +
  scale_color_manual(
    values = c("Brightened" = "#fdae61", "Dimmed" = "#878787", "No change" = "grey50"),
    name = "Direction") +
  facet_wrap(~ scenario, nrow = 1) +
  labs(title    = "Change in Greenness Index vs. Housing Density Change",
       # subtitle = "Positive = brightening (greener than expected); Negative = dimming",
       x = expression(Delta ~ "log(HU density/ha), 2010" %->% "2020"),
       y = expression(Delta ~ "Greenness Index")) +
  theme_lambert

# ggsave(file.path(out_dir, "fig2_change_green_vs_housing.png"),
#        fig2, width = 12, height = 4, dpi = 300)

# -------------------------------------------------------
# FIGURE 3: ΔGHI vs raw NDVI change
# -------------------------------------------------------
fig3 <- ggplot(all_data,
               aes(x = delta_NDVI, y = change_green)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_point(aes(color = direction_green), alpha = 0.15, size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7) +
  scale_color_manual(
    values = c("Brightened" = "#fdae61", "Dimmed" = "#878787", "No change" = "grey50"),
    name = "Direction") +
  facet_wrap(~ scenario, nrow = 1) +
  labs(title    = "Change in Greenness Index vs. Raw NDVI Change",
       # subtitle = "Analogous to Lambert Fig 3b: index change as a function of greenness loss/gain",
       x = expression(Delta ~ "NDVI, 2010" %->% "2020"),
       y = expression(Delta ~ "Greenness Index")) +
  theme_lambert

# ggsave(file.path(out_dir, "fig3_change_green_vs_ndvi.png"),
#        fig3, width = 12, height = 4, dpi = 300)

# -------------------------------------------------------
# FIGURE 4: Arrow plot — NDVI vs HU density trajectories
# -------------------------------------------------------
set.seed(42)
arrow_sample <- all_data %>%
  filter(abs(delta_NDVI) > 0.01 | abs(delta_log_HU) > 0.05) %>%
  group_by(scenario) %>%
  slice_sample(n = min(500, n())) %>%
  ungroup()

fig4 <- ggplot(arrow_sample) +
  geom_segment(
    aes(x = log_HU_density_2010, xend = log_HU_density_2020,
        y = NDVI_2010, yend = NDVI_2020,
        color = direction_green),
    arrow = arrow(length = unit(0.06, "inches"), type = "closed"),
    alpha = 0.35, linewidth = 0.3) +
  scale_color_manual(
    values = c("Brightened" = "#fdae61", "Dimmed" = "#878787", "No change" = "grey50"),
    name = "Index Direction") +
  facet_wrap(~ scenario, nrow = 1) +
  labs(title    = "Trajectories: Greenness vs. Housing Density, 2010 to 2020",
       subtitle = "Arrows show block-level change (subsample); color = index brightening/dimming",
       x = "log(Housing Unit Density per ha)",
       y = "NDVI") +
  theme_lambert

# ggsave(file.path(out_dir, "fig4_arrow_green_trajectories.png"),
#        fig4, width = 12, height = 4, dpi = 300)

# -------------------------------------------------------
# FIGURE 5: Boxplot of ΔGHI by scenario
# -------------------------------------------------------
fig5 <- ggplot(all_data, aes(x = scenario, y = change_green, fill = scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2, width = 0.6) +
  scale_fill_manual(values = c("Dormant 1000 m" = "#abd9e9",
                                "Growing 500 m"  = "#c7e9c0",
                                "Growing 1000 m" = "#66c164",
                                "Growing 2000 m" = "#1a9641"),
                    guide = "none") +
  labs(title    = "Distribution of ΔGHI by Scenario",
       subtitle = "Positive = brightening; Negative = dimming (relative to 2010 baseline)",
       x = NULL,
       y = "ΔGHI (pseudo-residual change)") +
  theme_lambert

# ggsave(file.path(out_dir, "fig5_boxplot_change.png"),
#        fig5, width = 8, height = 5, dpi = 300)

# -------------------------------------------------------
# FIGURE 6: Mean ΔGHI with 95% CI
# -------------------------------------------------------
ci_data <- all_data %>%
  group_by(scenario, season, buffer) %>%
  summarise(
    mean = mean(change_green),
    se   = sd(change_green) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(lo = mean - 1.96 * se,
         hi = mean + 1.96 * se)

fig6 <- ggplot(ci_data, aes(x = scenario, y = mean, color = scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = 0.6) +
  scale_color_manual(values = c("Dormant 1000 m" = "#2c7bb6",
                                 "Growing 500 m"  = "#74c476",
                                 "Growing 1000 m" = "#31a354",
                                 "Growing 2000 m" = "#006d2c"),
                     guide = "none") +
  labs(title = "Mean ΔGHI with 95% CI by Scenario",
       x = NULL, y = "Mean ΔGHI") +
  theme_lambert

# ggsave(file.path(out_dir, "fig6_mean_change_ci.png"),
#        fig6, width = 8, height = 5, dpi = 300)


# ============================================================
# SAVE FULL ANALYTICAL DATASET
# ============================================================
# write.csv(all_data,
#           file.path(out_dir, "block_pseudo_residual_analysis.csv"),
#           row.names = FALSE)

cat("\n================================================================\n")
cat("All outputs saved to:", out_dir, "\n")
cat("  Tables: table1, table2, table3 CSVs\n")
cat("  Figures: fig1-fig6 PNGs\n")
cat("  Data: block_pseudo_residual_analysis.csv\n")
cat("================================================================\n")
cat("ANALYSIS COMPLETE\n")
