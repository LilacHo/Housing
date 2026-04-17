library(here)
library(tidyverse)

here::i_am("code/1_block_scenario_comparison.R")

# ============================================================
# SCENARIO COMPARISON: Are the 6 season x buffer combos
# meaningfully different, or can some be collapsed?
# ============================================================

out_dir <- here::here("output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --- Load and prepare data (same as 1_block.R) ---
ndvi_block2 <- read.csv(here::here("data", "block_final.csv"))

df <- ndvi_block2 %>%
  mutate(
    across(all_of(c("PCT_BUILT", "NDVI", "POP_density", "HU_density")), as.numeric),
    year   = as.integer(year),
    buffer = as.integer(buffer)
  ) %>%
  filter(year %in% c(2010, 2020)) %>%
  filter(!is.na(PCT_BUILT), !is.na(NDVI),
         !is.na(POP_density), !is.na(HU_density)) %>%
  mutate(log_HU_density = log(pmax(HU_density, 1e-10)))

wide <- df %>%
  select(GEOID20, season, buffer, year, log_HU_density, NDVI, PCT_BUILT) %>%
  pivot_wider(names_from = year,
              values_from = c(log_HU_density, NDVI, PCT_BUILT),
              names_sep = "_") %>%
  drop_na()

cat("Total block-season-buffer observations:", nrow(wide), "\n")
cat("Unique blocks:", n_distinct(wide$GEOID20), "\n\n")


# ============================================================
# TEST 1: Are dormant and growing season NDVI linearly related?
#         If NDVI_dormant = a + b * NDVI_growing with high R²,
#         the two seasons carry redundant information.
# ============================================================

cat("============================================================\n")
cat("TEST 1: Cross-season NDVI correlation (within block x buffer)\n")
cat("============================================================\n\n")

season_wide <- wide %>%
  select(GEOID20, buffer, season, NDVI_2010, NDVI_2020) %>%
  pivot_wider(names_from = season,
              values_from = c(NDVI_2010, NDVI_2020),
              names_sep = ".") %>%
  drop_na()

for (b in sort(unique(season_wide$buffer))) {
  sub <- season_wide %>% filter(buffer == b)
  
  # 2010 cross-season
  r10 <- cor(sub$NDVI_2010.Dormant_Season, sub$NDVI_2010.Growing_Season)
  m10 <- lm(NDVI_2010.Dormant_Season ~ NDVI_2010.Growing_Season, data = sub)
  s10 <- summary(m10)
  
  # 2020 cross-season
  r20 <- cor(sub$NDVI_2020.Dormant_Season, sub$NDVI_2020.Growing_Season)
  m20 <- lm(NDVI_2020.Dormant_Season ~ NDVI_2020.Growing_Season, data = sub)
  s20 <- summary(m20)
  
  cat(sprintf("Buffer %d m:\n", b))
  cat(sprintf("  2010: r = %.4f, R² = %.4f, slope = %.4f, intercept = %.4f\n",
              r10, s10$r.squared, coef(m10)[2], coef(m10)[1]))
  cat(sprintf("  2020: r = %.4f, R² = %.4f, slope = %.4f, intercept = %.4f\n",
              r20, s20$r.squared, coef(m20)[2], coef(m20)[1]))
}

# Same for PCT_BUILT
cat("\n--- Cross-season PCT_BUILT correlation ---\n")
season_wide_imp <- wide %>%
  select(GEOID20, buffer, season, PCT_BUILT_2010, PCT_BUILT_2020) %>%
  pivot_wider(names_from = season,
              values_from = c(PCT_BUILT_2010, PCT_BUILT_2020),
              names_sep = ".") %>%
  drop_na()

for (b in sort(unique(season_wide_imp$buffer))) {
  sub <- season_wide_imp %>% filter(buffer == b)
  r10 <- cor(sub$PCT_BUILT_2010.Dormant_Season, sub$PCT_BUILT_2010.Growing_Season)
  r20 <- cor(sub$PCT_BUILT_2020.Dormant_Season, sub$PCT_BUILT_2020.Growing_Season)
  cat(sprintf("Buffer %d m: r_2010 = %.4f, r_2020 = %.4f\n", b, r10, r20))
}


# ============================================================
# TEST 2: Are the 2010 regression coefficients meaningfully
#         different across scenarios?
#         Compare slopes and intercepts with confidence intervals.
# ============================================================

cat("\n\n============================================================\n")
cat("TEST 2: Regression coefficient comparison across scenarios\n")
cat("============================================================\n\n")

combos <- expand.grid(
  season = c("Dormant_Season", "Growing_Season"),
  buffer = c(500, 1000, 2000),
  stringsAsFactors = FALSE
) %>% arrange(season, buffer)

reg_results <- pmap_dfr(combos, function(season, buffer) {
  sub <- wide %>% filter(season == !!season, buffer == !!buffer)
  mod_g <- lm(NDVI_2010 ~ log_HU_density_2010, data = sub)
  mod_i <- lm(PCT_BUILT_2010 ~ log_HU_density_2010, data = sub)
  tibble(
    season      = season,
    buffer      = buffer,
    n           = nrow(sub),
    g_intercept = coef(mod_g)[1],
    g_slope     = coef(mod_g)[2],
    g_slope_se  = summary(mod_g)$coefficients[2, 2],
    g_R2        = summary(mod_g)$r.squared,
    i_intercept = coef(mod_i)[1],
    i_slope     = coef(mod_i)[2],
    i_slope_se  = summary(mod_i)$coefficients[2, 2],
    i_R2        = summary(mod_i)$r.squared
  )
})

cat("--- Greenness (NDVI) regression coefficients ---\n")
reg_results %>%
  select(season, buffer, n, g_intercept, g_slope, g_slope_se, g_R2) %>%
  mutate(across(where(is.numeric) & !matches("^n$"), ~ round(.x, 5))) %>%
  print(n = Inf, width = Inf)

cat("\n--- Imperviousness (PCT_BUILT) regression coefficients ---\n")
reg_results %>%
  select(season, buffer, n, i_intercept, i_slope, i_slope_se, i_R2) %>%
  mutate(across(where(is.numeric) & !matches("^n$"), ~ round(.x, 5))) %>%
  print(n = Inf, width = Inf)

# Formal test: does adding season and/or buffer improve the model?
cat("\n--- Nested model comparison: does season matter for NDVI? ---\n")
for (b in c(500, 1000, 2000)) {
  sub <- wide %>% filter(buffer == b)
  m_pooled   <- lm(NDVI_2010 ~ log_HU_density_2010, data = sub)
  m_season   <- lm(NDVI_2010 ~ log_HU_density_2010 * season, data = sub)
  cat(sprintf("\nBuffer %d m:\n", b))
  print(anova(m_pooled, m_season))
}

cat("\n--- Nested model comparison: does buffer matter for NDVI? ---\n")
for (s in c("Dormant_Season", "Growing_Season")) {
  sub <- wide %>% filter(season == s)
  sub$buffer_f <- factor(sub$buffer)
  m_pooled <- lm(NDVI_2010 ~ log_HU_density_2010, data = sub)
  m_buffer <- lm(NDVI_2010 ~ log_HU_density_2010 * buffer_f, data = sub)
  cat(sprintf("\n%s:\n", s))
  print(anova(m_pooled, m_buffer))
}


# ============================================================
# TEST 3: Are the CHANGE variables correlated across scenarios?
#         If change_green for dormant ≈ change_green for growing,
#         the two seasons tell the same story about change.
# ============================================================

cat("\n\n============================================================\n")
cat("TEST 3: Cross-scenario correlation of CHANGE variables\n")
cat("============================================================\n\n")

# First, compute change variables per scenario (replicating 1_block.R logic)
compute_change <- function(d, s, b) {
  sub <- d %>% filter(season == s, buffer == b)
  mg <- lm(NDVI_2010 ~ log_HU_density_2010, data = sub)
  mi <- lm(PCT_BUILT_2010 ~ log_HU_density_2010, data = sub)
  sub %>% mutate(
    resid_g_2010 = NDVI_2010 - predict(mg),
    resid_i_2010 = PCT_BUILT_2010 - predict(mi),
    resid_g_2020 = NDVI_2020 - predict(mg, newdata = data.frame(
      log_HU_density_2010 = log_HU_density_2020)),
    resid_i_2020 = PCT_BUILT_2020 - predict(mi, newdata = data.frame(
      log_HU_density_2010 = log_HU_density_2020)),
    change_green = resid_g_2020 - resid_g_2010,
    change_imp   = resid_i_2020 - resid_i_2010
  )
}

all_changes <- map2_dfr(
  rep(c("Dormant_Season", "Growing_Season"), each = 3),
  rep(c(500, 1000, 2000), 2),
  ~ compute_change(wide, .x, .y)
)

# Pivot to get one column per scenario for each block
change_wide_g <- all_changes %>%
  mutate(scenario = paste(season, buffer, sep = "_")) %>%
  select(GEOID20, scenario, change_green) %>%
  pivot_wider(names_from = scenario, values_from = change_green,
              names_prefix = "g_") %>%
  drop_na()

change_wide_i <- all_changes %>%
  mutate(scenario = paste(season, buffer, sep = "_")) %>%
  select(GEOID20, scenario, change_imp) %>%
  pivot_wider(names_from = scenario, values_from = change_imp,
              names_prefix = "i_") %>%
  drop_na()

cat("--- Greenness change correlation matrix ---\n")
cor_g <- cor(change_wide_g %>% select(-GEOID20))
print(round(cor_g, 4))

cat("\n--- Imperviousness change correlation matrix ---\n")
cor_i <- cor(change_wide_i %>% select(-GEOID20))
print(round(cor_i, 4))


# ============================================================
# TEST 4: Visualize cross-scenario relationships
# ============================================================

library(ggplot2)

theme_diag <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none")

# --- 4a: Dormant vs Growing NDVI change, by buffer ---
cat("\n\n============================================================\n")
cat("TEST 4: Cross-season scatter plots of change variables\n")
cat("============================================================\n\n")

change_cross <- all_changes %>%
  select(GEOID20, buffer, season, change_green, change_imp) %>%
  pivot_wider(names_from = season,
              values_from = c(change_green, change_imp),
              names_sep = ".") %>%
  drop_na()

for (b in c(500, 1000, 2000)) {
  sub <- change_cross %>% filter(buffer == b)
  r_g <- cor(sub$change_green.Dormant_Season, sub$change_green.Growing_Season)
  r_i <- cor(sub$change_imp.Dormant_Season, sub$change_imp.Growing_Season)
  cat(sprintf("Buffer %d m: r(change_green) = %.4f, r(change_imp) = %.4f\n",
              b, r_g, r_i))
}

fig_cross_green <- ggplot(change_cross,
       aes(x = change_green.Dormant_Season,
           y = change_green.Growing_Season)) +
  geom_point(alpha = 0.05, size = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", se = TRUE, color = "#2c7bb6", linewidth = 0.7) +
  facet_wrap(~ buffer, labeller = labeller(
    buffer = c("500" = "500 m", "1000" = "1000 m", "2000" = "2000 m"))) +
  labs(
    title = "Greenness Index Change: Dormant vs. Growing Season",
    subtitle = "Red dashed = 1:1 line; Blue = OLS fit",
    x = "Dormant Season Change",
    y = "Growing Season Change"
  ) +
  coord_equal() +
  theme_diag

ggsave(file.path(out_dir, "diag_cross_season_green_change.png"),
       fig_cross_green, width = 10, height = 4, dpi = 300)

fig_cross_imp <- ggplot(change_cross,
       aes(x = change_imp.Dormant_Season,
           y = change_imp.Growing_Season)) +
  geom_point(alpha = 0.05, size = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", se = TRUE, color = "#d7191c", linewidth = 0.7) +
  facet_wrap(~ buffer, labeller = labeller(
    buffer = c("500" = "500 m", "1000" = "1000 m", "2000" = "2000 m"))) +
  labs(
    title = "Imperviousness Index Change: Dormant vs. Growing Season",
    subtitle = "Red dashed = 1:1 line; Blue = OLS fit",
    x = "Dormant Season Change",
    y = "Growing Season Change"
  ) +
  coord_equal() +
  theme_diag

ggsave(file.path(out_dir, "diag_cross_season_imp_change.png"),
       fig_cross_imp, width = 10, height = 4, dpi = 300)


# --- 4b: Cross-buffer scatter plots (within same season) ---

change_cross_buf <- all_changes %>%
  select(GEOID20, buffer, season, change_green, change_imp) %>%
  pivot_wider(names_from = buffer,
              values_from = c(change_green, change_imp),
              names_sep = ".") %>%
  drop_na()

cat("\n--- Cross-buffer correlations for greenness change ---\n")
for (s in c("Dormant_Season", "Growing_Season")) {
  sub <- change_cross_buf %>% filter(season == s)
  cat(sprintf("%s:\n", s))
  cat(sprintf("  500 vs 1000: r = %.4f\n",
              cor(sub$change_green.500, sub$change_green.1000)))
  cat(sprintf("  500 vs 2000: r = %.4f\n",
              cor(sub$change_green.500, sub$change_green.2000)))
  cat(sprintf("  1000 vs 2000: r = %.4f\n",
              cor(sub$change_green.1000, sub$change_green.2000)))
}

cat("\n--- Cross-buffer correlations for imperviousness change ---\n")
for (s in c("Dormant_Season", "Growing_Season")) {
  sub <- change_cross_buf %>% filter(season == s)
  cat(sprintf("%s:\n", s))
  cat(sprintf("  500 vs 1000: r = %.4f\n",
              cor(sub$change_imp.500, sub$change_imp.1000)))
  cat(sprintf("  500 vs 2000: r = %.4f\n",
              cor(sub$change_imp.500, sub$change_imp.2000)))
  cat(sprintf("  1000 vs 2000: r = %.4f\n",
              cor(sub$change_imp.1000, sub$change_imp.2000)))
}

# Cross-buffer scatter for greenness
fig_cross_buf_green <- change_cross_buf %>%
  select(GEOID20, season, change_green.500, change_green.1000, change_green.2000) %>%
  pivot_longer(cols = c(change_green.1000, change_green.2000),
               names_to = "comparison", values_to = "y") %>%
  mutate(comparison = recode(comparison,
    change_green.1000 = "500 m vs 1000 m",
    change_green.2000 = "500 m vs 2000 m")) %>%
  ggplot(aes(x = change_green.500, y = y)) +
  geom_point(alpha = 0.03, size = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", se = TRUE, color = "#2c7bb6", linewidth = 0.7) +
  facet_grid(season ~ comparison) +
  labs(
    title = "Greenness Change: Cross-Buffer Comparison",
    x = "Change at 500 m buffer",
    y = "Change at comparison buffer"
  ) +
  coord_equal() +
  theme_diag

ggsave(file.path(out_dir, "diag_cross_buffer_green_change.png"),
       fig_cross_buf_green, width = 8, height = 6, dpi = 300)


# ============================================================
# TEST 5: Information-theoretic comparison
#         Does a model with season/buffer interactions fit
#         meaningfully better than a pooled model?
# ============================================================

cat("\n\n============================================================\n")
cat("TEST 5: AIC/BIC comparison — pooled vs stratified models\n")
cat("============================================================\n\n")

# For NDVI
all_for_aic <- wide %>%
  mutate(buffer_f = factor(buffer))

cat("--- NDVI_2010 ~ log_HU_density_2010 ---\n")
m1 <- lm(NDVI_2010 ~ log_HU_density_2010, data = all_for_aic)
m2 <- lm(NDVI_2010 ~ log_HU_density_2010 + season, data = all_for_aic)
m3 <- lm(NDVI_2010 ~ log_HU_density_2010 + buffer_f, data = all_for_aic)
m4 <- lm(NDVI_2010 ~ log_HU_density_2010 + season + buffer_f, data = all_for_aic)
m5 <- lm(NDVI_2010 ~ log_HU_density_2010 * season, data = all_for_aic)
m6 <- lm(NDVI_2010 ~ log_HU_density_2010 * buffer_f, data = all_for_aic)
m7 <- lm(NDVI_2010 ~ log_HU_density_2010 * season * buffer_f, data = all_for_aic)

aic_table_ndvi <- tibble(
  model = c("pooled", "+season", "+buffer", "+season+buffer",
            "*season", "*buffer", "*season*buffer"),
  AIC   = c(AIC(m1), AIC(m2), AIC(m3), AIC(m4), AIC(m5), AIC(m6), AIC(m7)),
  BIC   = c(BIC(m1), BIC(m2), BIC(m3), BIC(m4), BIC(m5), BIC(m6), BIC(m7)),
  R2    = c(summary(m1)$r.squared, summary(m2)$r.squared,
            summary(m3)$r.squared, summary(m4)$r.squared,
            summary(m5)$r.squared, summary(m6)$r.squared,
            summary(m7)$r.squared)
) %>%
  mutate(delta_AIC = AIC - min(AIC),
         delta_BIC = BIC - min(BIC))

cat("NDVI model comparison:\n")
print(as.data.frame(aic_table_ndvi), digits = 2, row.names = FALSE)

# For PCT_BUILT
cat("\n--- PCT_BUILT_2010 ~ log_HU_density_2010 ---\n")
p1 <- lm(PCT_BUILT_2010 ~ log_HU_density_2010, data = all_for_aic)
p2 <- lm(PCT_BUILT_2010 ~ log_HU_density_2010 + season, data = all_for_aic)
p3 <- lm(PCT_BUILT_2010 ~ log_HU_density_2010 + buffer_f, data = all_for_aic)
p4 <- lm(PCT_BUILT_2010 ~ log_HU_density_2010 + season + buffer_f, data = all_for_aic)
p5 <- lm(PCT_BUILT_2010 ~ log_HU_density_2010 * season, data = all_for_aic)
p6 <- lm(PCT_BUILT_2010 ~ log_HU_density_2010 * buffer_f, data = all_for_aic)
p7 <- lm(PCT_BUILT_2010 ~ log_HU_density_2010 * season * buffer_f, data = all_for_aic)

aic_table_imp <- tibble(
  model = c("pooled", "+season", "+buffer", "+season+buffer",
            "*season", "*buffer", "*season*buffer"),
  AIC   = c(AIC(p1), AIC(p2), AIC(p3), AIC(p4), AIC(p5), AIC(p6), AIC(p7)),
  BIC   = c(BIC(p1), BIC(p2), BIC(p3), BIC(p4), BIC(p5), BIC(p6), BIC(p7)),
  R2    = c(summary(p1)$r.squared, summary(p2)$r.squared,
            summary(p3)$r.squared, summary(p4)$r.squared,
            summary(p5)$r.squared, summary(p6)$r.squared,
            summary(p7)$r.squared)
) %>%
  mutate(delta_AIC = AIC - min(AIC),
         delta_BIC = BIC - min(BIC))

cat("PCT_BUILT model comparison:\n")
print(as.data.frame(aic_table_imp), digits = 2, row.names = FALSE)

write.csv(aic_table_ndvi,
          file.path(out_dir, "diag_aic_ndvi.csv"), row.names = FALSE)
write.csv(aic_table_imp,
          file.path(out_dir, "diag_aic_imp.csv"), row.names = FALSE)


# ============================================================
# TEST 6: Practical significance — how different are the
#         change variable means across scenarios?
#         Cohen's d for season effect within each buffer.
# ============================================================

cat("\n\n============================================================\n")
cat("TEST 6: Effect sizes (Cohen's d) for season differences\n")
cat("============================================================\n\n")

cohens_d <- function(x, y) {
  nx <- length(x); ny <- length(y)
  sp <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  (mean(x) - mean(y)) / sp
}

for (b in c(500, 1000, 2000)) {
  sub_d <- all_changes %>% filter(buffer == b, season == "Dormant_Season")
  sub_g <- all_changes %>% filter(buffer == b, season == "Growing_Season")
  
  d_green <- cohens_d(sub_d$change_green, sub_g$change_green)
  d_imp   <- cohens_d(sub_d$change_imp, sub_g$change_imp)
  
  cat(sprintf("Buffer %d m:\n", b))
  cat(sprintf("  Greenness change: Cohen's d = %.4f  (%s)\n",
              d_green,
              ifelse(abs(d_green) < 0.2, "negligible",
              ifelse(abs(d_green) < 0.5, "small",
              ifelse(abs(d_green) < 0.8, "medium", "large")))))
  cat(sprintf("  Imperviousness change: Cohen's d = %.4f  (%s)\n",
              d_imp,
              ifelse(abs(d_imp) < 0.2, "negligible",
              ifelse(abs(d_imp) < 0.5, "small",
              ifelse(abs(d_imp) < 0.8, "medium", "large")))))
}

# Cross-buffer effect sizes within each season
cat("\n--- Effect sizes for buffer differences (500 vs 2000) ---\n")
for (s in c("Dormant_Season", "Growing_Season")) {
  sub_500  <- all_changes %>% filter(season == s, buffer == 500)
  sub_2000 <- all_changes %>% filter(season == s, buffer == 2000)
  
  d_green <- cohens_d(sub_500$change_green, sub_2000$change_green)
  d_imp   <- cohens_d(sub_500$change_imp, sub_2000$change_imp)
  
  cat(sprintf("%s (500 vs 2000 m):\n", s))
  cat(sprintf("  Greenness: Cohen's d = %.4f  (%s)\n",
              d_green,
              ifelse(abs(d_green) < 0.2, "negligible",
              ifelse(abs(d_green) < 0.5, "small",
              ifelse(abs(d_green) < 0.8, "medium", "large")))))
  cat(sprintf("  Imperviousness: Cohen's d = %.4f  (%s)\n",
              d_imp,
              ifelse(abs(d_imp) < 0.2, "negligible",
              ifelse(abs(d_imp) < 0.5, "small",
              ifelse(abs(d_imp) < 0.8, "medium", "large")))))
}


# ============================================================
# TEST 7: Summary visualization — coefficient comparison plot
# ============================================================

coef_plot_data <- reg_results %>%
  select(season, buffer, g_slope, g_slope_se, g_R2,
         i_slope, i_slope_se, i_R2) %>%
  pivot_longer(
    cols = c(g_slope, g_slope_se, g_R2, i_slope, i_slope_se, i_R2),
    names_to = c("variable", ".value"),
    names_pattern = "(g|i)_(slope|slope_se|R2)"
  ) %>%
  mutate(
    variable = recode(variable, g = "Greenness (NDVI)", i = "Imperviousness (PCT_BUILT)"),
    lo = slope - 1.96 * slope_se,
    hi = slope + 1.96 * slope_se,
    buffer_label = paste0(buffer, " m")
  )

fig_coef <- ggplot(coef_plot_data,
       aes(x = buffer_label, y = slope, color = season)) +
  geom_pointrange(aes(ymin = lo, ymax = hi),
                  position = position_dodge(width = 0.4), size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ variable, scales = "free_y") +
  scale_color_manual(values = c("Dormant_Season" = "#2c7bb6",
                                 "Growing_Season" = "#1a9641"),
                     name = "Season",
                     labels = c("Dormant", "Growing")) +
  labs(
    title = "2010 Regression Slopes (± 95% CI) Across Scenarios",
    subtitle = "Overlap indicates scenarios may be collapsible",
    x = "Buffer", y = "Slope of log(HU_density)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "diag_coefficient_comparison.png"),
       fig_coef, width = 8, height = 5, dpi = 300)

# R² comparison
fig_r2 <- coef_plot_data %>%
  ggplot(aes(x = buffer_label, y = R2, fill = season)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  facet_wrap(~ variable) +
  scale_fill_manual(values = c("Dormant_Season" = "#abd9e9",
                                "Growing_Season" = "#66c164"),
                    name = "Season",
                    labels = c("Dormant", "Growing")) +
  labs(
    title = "R² of 2010 Baseline Regressions Across Scenarios",
    x = "Buffer", y = "R²"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "diag_r2_comparison.png"),
       fig_r2, width = 8, height = 5, dpi = 300)

cat("\n\n============================================================\n")
cat("SCENARIO COMPARISON COMPLETE\n")
cat("Outputs saved to:", out_dir, "\n")
cat("============================================================\n")
