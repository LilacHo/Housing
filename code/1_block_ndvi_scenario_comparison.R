library(here)
library(tidyverse)

here::i_am("code/1_block_ndvi_scenario_comparison.R")

# ============================================================
# SCENARIO COMPARISON (GREENNESS ONLY)
# Are the 6 season x buffer combos meaningfully different,
# or can some be collapsed?
#
# Focuses exclusively on NDVI; PCT_BUILT is ignored.
# ============================================================

out_dir <- here::here("output")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# -----------------------------
# 1. Load and prepare data
# -----------------------------
ndvi_block2 <- read.csv(here::here("data", "block_final.csv"))

df <- ndvi_block2 %>%
  mutate(
    across(all_of(c("NDVI", "HU_density")), as.numeric),
    year   = as.integer(year),
    buffer = as.integer(buffer)
  ) %>%
  filter(year %in% c(2010, 2020)) %>%
  filter(!is.na(NDVI), !is.na(HU_density)) %>%
  mutate(log_HU_density = log(pmax(HU_density, 1e-10)))

wide <- df %>%
  select(GEOID20, season, buffer, year, log_HU_density, NDVI) %>%
  pivot_wider(names_from  = year,
              values_from = c(log_HU_density, NDVI),
              names_sep   = "_") %>%
  drop_na() %>%
  mutate(delta_log_HU = log_HU_density_2020 - log_HU_density_2010,
         delta_NDVI   = NDVI_2020 - NDVI_2010)

cat("Total block-season-buffer observations:", nrow(wide), "\n")
cat("Unique blocks:", n_distinct(wide$GEOID20), "\n\n")


# ============================================================
# TEST 1: Cross-season NDVI correlation (within block x buffer)
# ============================================================
cat("============================================================\n")
cat("TEST 1: Cross-season NDVI correlation (within block x buffer)\n")
cat("============================================================\n\n")

season_wide <- wide %>%
  select(GEOID20, buffer, season, NDVI_2010, NDVI_2020) %>%
  pivot_wider(names_from  = season,
              values_from = c(NDVI_2010, NDVI_2020),
              names_sep   = ".") %>%
  drop_na()

for (b in sort(unique(season_wide$buffer))) {
  sub <- season_wide %>% filter(buffer == b)
  r10 <- cor(sub$NDVI_2010.Dormant_Season, sub$NDVI_2010.Growing_Season)
  r20 <- cor(sub$NDVI_2020.Dormant_Season, sub$NDVI_2020.Growing_Season)
  cat(sprintf("Buffer %d m: r_2010 = %.4f, r_2020 = %.4f\n", b, r10, r20))
}


# ============================================================
# TEST 2: Regression coefficient comparison across scenarios
# ============================================================
cat("\n\n============================================================\n")
cat("TEST 2: 2010 regression coefficients by season x buffer\n")
cat("============================================================\n\n")

combos <- expand.grid(
  season = c("Dormant_Season", "Growing_Season"),
  buffer = c(500, 1000, 2000),
  stringsAsFactors = FALSE
) %>% arrange(season, buffer)

reg_results <- pmap_dfr(combos, function(season, buffer) {
  sub <- wide %>% filter(season == !!season, buffer == !!buffer)
  m <- lm(NDVI_2010 ~ log_HU_density_2010, data = sub)
  tibble(
    season    = season,
    buffer    = buffer,
    n         = nrow(sub),
    intercept = coef(m)[1],
    slope     = coef(m)[2],
    slope_se  = summary(m)$coefficients[2, 2],
    slope_p   = summary(m)$coefficients[2, 4],
    R2        = summary(m)$r.squared
  )
})

cat("--- Greenness (NDVI) regression coefficients ---\n")
print(reg_results %>% mutate(across(where(is.numeric) & !matches("^n$"),
                                    ~ round(.x, 5))),
      n = Inf, width = Inf)

# Nested model comparisons
cat("\n--- Nested model: does season matter for NDVI? ---\n")
for (b in c(500, 1000, 2000)) {
  sub <- wide %>% filter(buffer == b)
  m_pooled <- lm(NDVI_2010 ~ log_HU_density_2010, data = sub)
  m_season <- lm(NDVI_2010 ~ log_HU_density_2010 * season, data = sub)
  cat(sprintf("\nBuffer %d m:\n", b))
  print(anova(m_pooled, m_season))
}

cat("\n--- Nested model: does buffer matter for NDVI? ---\n")
for (s in c("Dormant_Season", "Growing_Season")) {
  sub <- wide %>% filter(season == s)
  sub$buffer_f <- factor(sub$buffer)
  m_pooled <- lm(NDVI_2010 ~ log_HU_density_2010, data = sub)
  m_buffer <- lm(NDVI_2010 ~ log_HU_density_2010 * buffer_f, data = sub)
  cat(sprintf("\n%s:\n", s))
  print(anova(m_pooled, m_buffer))
}


# ============================================================
# TEST 3: Compute change variables per scenario, then correlate
# ============================================================
cat("\n\n============================================================\n")
cat("TEST 3: Cross-scenario correlation of greenness CHANGE\n")
cat("============================================================\n\n")

compute_change <- function(d, s, b) {
  sub <- d %>% filter(season == s, buffer == b)
  m <- lm(NDVI_2010 ~ log_HU_density_2010, data = sub)
  sub %>% mutate(
    resid_g_2010 = NDVI_2010 - predict(m),
    resid_g_2020 = NDVI_2020 - predict(m, newdata = data.frame(
      log_HU_density_2010 = log_HU_density_2020)),
    change_green = resid_g_2020 - resid_g_2010
  )
}

all_changes <- map2_dfr(
  rep(c("Dormant_Season", "Growing_Season"), each = 3),
  rep(c(500, 1000, 2000), 2),
  ~ compute_change(wide, .x, .y)
)

change_wide <- all_changes %>%
  mutate(scenario = paste(season, buffer, sep = "_")) %>%
  select(GEOID20, scenario, change_green) %>%
  pivot_wider(names_from = scenario, values_from = change_green,
              names_prefix = "g_") %>%
  drop_na()

cat("Greenness change correlation matrix:\n")
cor_g <- cor(change_wide %>% select(-GEOID20))
print(round(cor_g, 4))


# ============================================================
# TEST 4: Cross-season and cross-buffer scatter plots
# ============================================================
cat("\n\n============================================================\n")
cat("TEST 4: Cross-scenario scatter plots of change variables\n")
cat("============================================================\n\n")

change_cross <- all_changes %>%
  select(GEOID20, buffer, season, change_green) %>%
  pivot_wider(names_from  = season,
              values_from = change_green,
              names_sep   = ".") %>%
  drop_na()

for (b in c(500, 1000, 2000)) {
  sub <- change_cross %>% filter(buffer == b)
  r_g <- cor(sub$Dormant_Season, sub$Growing_Season)
  cat(sprintf("Buffer %d m: r(change_green, Dormant vs Growing) = %.4f\n",
              b, r_g))
}

theme_diag <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

fig_cross_green <- ggplot(change_cross,
       aes(x = Dormant_Season, y = Growing_Season)) +
  geom_point(alpha = 0.05, size = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", se = TRUE, color = "#2c7bb6", linewidth = 0.7) +
  facet_wrap(~ buffer, labeller = labeller(
    buffer = c("500" = "500 m", "1000" = "1000 m", "2000" = "2000 m"))) +
  labs(
    title    = "Greenness Index Change: Dormant vs. Growing Season",
    subtitle = "Red dashed = 1:1 line; Blue = OLS fit",
    x = "Dormant Season Change", y = "Growing Season Change"
  ) +
  coord_equal() + theme_diag

ggsave(file.path(out_dir, "diag_cross_season_green_change.png"),
       fig_cross_green, width = 10, height = 4, dpi = 300)

# Cross-buffer correlations within each season
change_cross_buf <- all_changes %>%
  select(GEOID20, buffer, season, change_green) %>%
  pivot_wider(names_from = buffer, values_from = change_green,
              names_prefix = "b") %>%
  drop_na()

cat("\n--- Cross-buffer correlations for greenness change ---\n")
for (s in c("Dormant_Season", "Growing_Season")) {
  sub <- change_cross_buf %>% filter(season == s)
  cat(sprintf("%s:\n", s))
  cat(sprintf("  500 vs 1000: r = %.4f\n",  cor(sub$b500,  sub$b1000)))
  cat(sprintf("  500 vs 2000: r = %.4f\n",  cor(sub$b500,  sub$b2000)))
  cat(sprintf("  1000 vs 2000: r = %.4f\n", cor(sub$b1000, sub$b2000)))
}


# ============================================================
# TEST 5: AIC/BIC model comparison
# ============================================================
cat("\n\n============================================================\n")
cat("TEST 5: AIC/BIC comparison — pooled vs stratified NDVI models\n")
cat("============================================================\n\n")

all_for_aic <- wide %>% mutate(buffer_f = factor(buffer))

m1 <- lm(NDVI_2010 ~ log_HU_density_2010, data = all_for_aic)
m2 <- lm(NDVI_2010 ~ log_HU_density_2010 + season, data = all_for_aic)
m3 <- lm(NDVI_2010 ~ log_HU_density_2010 + buffer_f, data = all_for_aic)
m4 <- lm(NDVI_2010 ~ log_HU_density_2010 + season + buffer_f, data = all_for_aic)
m5 <- lm(NDVI_2010 ~ log_HU_density_2010 * season, data = all_for_aic)
m6 <- lm(NDVI_2010 ~ log_HU_density_2010 * buffer_f, data = all_for_aic)
m7 <- lm(NDVI_2010 ~ log_HU_density_2010 * season * buffer_f, data = all_for_aic)

aic_table <- tibble(
  model = c("pooled", "+season", "+buffer", "+season+buffer",
            "*season", "*buffer", "*season*buffer"),
  AIC   = c(AIC(m1), AIC(m2), AIC(m3), AIC(m4), AIC(m5), AIC(m6), AIC(m7)),
  BIC   = c(BIC(m1), BIC(m2), BIC(m3), BIC(m4), BIC(m5), BIC(m6), BIC(m7)),
  R2    = c(summary(m1)$r.squared, summary(m2)$r.squared,
            summary(m3)$r.squared, summary(m4)$r.squared,
            summary(m5)$r.squared, summary(m6)$r.squared,
            summary(m7)$r.squared)
) %>% mutate(delta_AIC = AIC - min(AIC),
             delta_BIC = BIC - min(BIC))

print(as.data.frame(aic_table), digits = 2, row.names = FALSE)
write.csv(aic_table, file.path(out_dir, "scenario_aic_ndvi.csv"),
          row.names = FALSE)


# ============================================================
# TEST 6: Cohen's d effect sizes for change variables
# ============================================================
cat("\n\n============================================================\n")
cat("TEST 6: Effect sizes (Cohen's d) for greenness change\n")
cat("============================================================\n\n")

cohens_d <- function(x, y) {
  nx <- length(x); ny <- length(y)
  sp <- sqrt(((nx - 1) * var(x) + (ny - 1) * var(y)) / (nx + ny - 2))
  (mean(x) - mean(y)) / sp
}

label_d <- function(d) {
  ifelse(abs(d) < 0.2, "negligible",
  ifelse(abs(d) < 0.5, "small",
  ifelse(abs(d) < 0.8, "medium", "large")))
}

cat("--- Season effect within each buffer ---\n")
for (b in c(500, 1000, 2000)) {
  x <- all_changes %>% filter(buffer == b, season == "Dormant_Season") %>% pull(change_green)
  y <- all_changes %>% filter(buffer == b, season == "Growing_Season") %>% pull(change_green)
  d <- cohens_d(x, y)
  cat(sprintf("Buffer %d m: Cohen's d = %.4f (%s)\n", b, d, label_d(d)))
}

cat("\n--- Buffer effect (500 vs 2000) within each season ---\n")
for (s in c("Dormant_Season", "Growing_Season")) {
  x <- all_changes %>% filter(season == s, buffer == 500)  %>% pull(change_green)
  y <- all_changes %>% filter(season == s, buffer == 2000) %>% pull(change_green)
  d <- cohens_d(x, y)
  cat(sprintf("%s: Cohen's d = %.4f (%s)\n", s, d, label_d(d)))
}


# ============================================================
# TEST 7: Coefficient comparison plot
# ============================================================
coef_plot_data <- reg_results %>%
  mutate(lo = slope - 1.96 * slope_se,
         hi = slope + 1.96 * slope_se,
         buffer_label = paste0(buffer, " m"))

fig_coef <- ggplot(coef_plot_data,
       aes(x = buffer_label, y = slope, color = season)) +
  geom_pointrange(aes(ymin = lo, ymax = hi),
                  position = position_dodge(width = 0.4), size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_color_manual(values = c("Dormant_Season" = "#2c7bb6",
                                 "Growing_Season" = "#1a9641"),
                     name = "Season", labels = c("Dormant", "Growing")) +
  labs(
    title    = "2010 NDVI Regression Slopes (± 95% CI)",
    subtitle = "Overlap indicates scenarios may be collapsible",
    x = "Buffer", y = "Slope of log(HU_density)"
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "scenario_coefficient_comparison_green.png"),
       fig_coef, width = 8, height = 5, dpi = 300)

fig_r2 <- coef_plot_data %>%
  ggplot(aes(x = buffer_label, y = R2, fill = season)) +
  geom_col(position = position_dodge(0.7), width = 0.6) +
  scale_fill_manual(values = c("Dormant_Season" = "#abd9e9",
                                "Growing_Season" = "#66c164"),
                    name = "Season", labels = c("Dormant", "Growing")) +
  labs(title = "R² of 2010 NDVI Baseline Regressions",
       x = "Buffer", y = "R²") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

ggsave(file.path(out_dir, "scenario_r2_comparison_green.png"),
       fig_r2, width = 8, height = 5, dpi = 300)

cat("\n\n============================================================\n")
cat("SCENARIO COMPARISON (GREENNESS) COMPLETE\n")
cat("Outputs saved to:", out_dir, "\n")
cat("============================================================\n")
