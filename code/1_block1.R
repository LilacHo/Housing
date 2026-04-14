library(here)
library(tidyverse)

here::i_am("code/1_block(1).R")

# ============================================================
# PSEUDO-RESIDUAL CHANGE ANALYSIS, 2010 -> 2020
# Greenness (NDVI) and Imperviousness (PCT_BUILT)
# relative to housing density (HU_density)
#
# Modeled after Lambert et al. (2025) Conservation Science
# and Practice: "Building the neighborhood for the trees"
# ============================================================

# -----------------------------
# 0. Packages
# -----------------------------
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(stringr)
library(ggplot2)

# Output directory for figures and tables
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
    across(all_of(c("PCT_BUILT", "NDVI", "POP_density", "HU_density")), as.numeric),
    year   = as.integer(year),
    buffer = as.integer(buffer)
  ) %>%
  filter(year %in% c(2010, 2020)) %>%
  filter(!is.na(PCT_BUILT), !is.na(NDVI),
         !is.na(POP_density), !is.na(HU_density)) %>%
  mutate(
    log_HU_density = log(pmax(HU_density, 1e-10))
  )

id_vars <- c("GEOID20", "season", "buffer")

# -----------------------------
# 3. PIVOT TO WIDE: one row per GEOID20 x season x buffer
# -----------------------------

wide <- df %>%
  select(all_of(id_vars), year,
         log_HU_density, NDVI, PCT_BUILT) %>%
  pivot_wider(
    names_from  = year,
    values_from = c(log_HU_density, NDVI, PCT_BUILT),
    names_sep   = "_"
  ) %>%
  drop_na()

# Derived change columns (raw, before residual analysis)
wide <- wide %>%
  mutate(
    delta_log_HU  = log_HU_density_2020 - log_HU_density_2010,
    delta_NDVI    = NDVI_2020 - NDVI_2010,
    delta_IMP     = PCT_BUILT_2020 - PCT_BUILT_2010
  )

# ============================================================
# 4. FIT 2010 REGRESSIONS (by season x buffer)
# ============================================================

combos <- wide %>%
  distinct(season, buffer) %>%
  arrange(season, buffer)

results <- combos %>%
  mutate(data = map2(season, buffer, ~ wide %>%
                       filter(season == .x, buffer == .y)))

results <- results %>%
  mutate(
    mod_green = map(data, ~ lm(NDVI_2010 ~ log_HU_density_2010, data = .x)),
    mod_imp   = map(data, ~ lm(PCT_BUILT_2010 ~ log_HU_density_2010, data = .x))
  )

# -------------------------------------------------------
# TABLE 1: 2010 Regression coefficients (Lambert-style)
# -------------------------------------------------------
reg_table <- results %>%
  mutate(
    tidy_green = map(mod_green, ~ {
      s <- summary(.x)
      tibble(
        green_intercept = coef(.x)[1],
        green_slope     = coef(.x)[2],
        green_slope_se  = coef(s)[2, 2],
        green_slope_p   = coef(s)[2, 4],
        green_R2        = s$r.squared,
        green_n         = nobs(.x)
      )
    }),
    tidy_imp = map(mod_imp, ~ {
      s <- summary(.x)
      tibble(
        imp_intercept = coef(.x)[1],
        imp_slope     = coef(.x)[2],
        imp_slope_se  = coef(s)[2, 2],
        imp_slope_p   = coef(s)[2, 4],
        imp_R2        = s$r.squared,
        imp_n         = nobs(.x)
      )
    })
  ) %>%
  select(season, buffer, tidy_green, tidy_imp) %>%
  unnest(c(tidy_green, tidy_imp))

cat("\n================================================================\n")
cat("TABLE 1: 2010 OLS Regression Coefficients by Season and Buffer\n")
cat("  Response ~ log(HU_density)\n")
cat("================================================================\n")
print(as.data.frame(reg_table), digits = 4, row.names = FALSE)

write.csv(reg_table,
          file.path(out_dir, "table1_regression_coefficients.csv"),
          row.names = FALSE)


# ============================================================
# 5. COMPUTE RESIDUALS (2010) AND PSEUDO-RESIDUALS (2020)
#    Following Lambert et al. (2025) methodology:
#    THI_2010 = y_2010 - (beta_2006 * log(x_2010) + alpha_2006)
#    THI_2020 = y_2020 - (beta_2006 * log(x_2020) + alpha_2006)
#    delta_THI = THI_2020 - THI_2010
# ============================================================

results <- results %>%
  mutate(
    data = pmap(list(data, mod_green, mod_imp), function(d, mg, mi) {
      d %>%
        mutate(
          # 2010 residuals (baseline index)
          resid_green_2010 = NDVI_2010 - predict(mg, newdata = d),
          resid_imp_2010   = PCT_BUILT_2010 - predict(mi, newdata = d),

          # 2020 predicted using 2010 model coefficients + 2020 housing
          pred_green_2020  = predict(mg, newdata = data.frame(
            log_HU_density_2010 = d$log_HU_density_2020)),
          pred_imp_2020    = predict(mi, newdata = data.frame(
            log_HU_density_2010 = d$log_HU_density_2020)),

          # 2020 pseudo-residuals
          resid_green_2020 = NDVI_2020 - pred_green_2020,
          resid_imp_2020   = PCT_BUILT_2020 - pred_imp_2020,

          # Change: brightening (+) / dimming (-)
          change_green = resid_green_2020 - resid_green_2010,
          change_imp   = resid_imp_2020   - resid_imp_2010
        )
    })
  )

all_data <- results %>%
  select(data) %>%
  unnest(data)

# Classify bright/dim states (Lambert-style)
all_data <- all_data %>%
  mutate(
    state_green_2010 = ifelse(resid_green_2010 >= 0, "Bright", "Dim"),
    state_green_2020 = ifelse(resid_green_2020 >= 0, "Bright", "Dim"),
    state_imp_2010   = ifelse(resid_imp_2010 >= 0, "Bright", "Dim"),
    state_imp_2020   = ifelse(resid_imp_2020 >= 0, "Bright", "Dim"),
    transition_green = paste0(state_green_2010, " -> ", state_green_2020),
    transition_imp   = paste0(state_imp_2010, " -> ", state_imp_2020),
    direction_green  = case_when(
      change_green > 0 ~ "Brightened",
      change_green < 0 ~ "Dimmed",
      TRUE ~ "No change"
    ),
    direction_imp = case_when(
      change_imp > 0 ~ "Brightened",
      change_imp < 0 ~ "Dimmed",
      TRUE ~ "No change"
    )
  )


# ============================================================
# TABLE 2: Summary of change variables by season x buffer
#          (analogous to Lambert's summary statistics)
# ============================================================

change_summary_green <- all_data %>%
  group_by(season, buffer) %>%
  summarise(
    n = n(),
    mean_change_green   = mean(change_green),
    sd_change_green     = sd(change_green),
    median_change_green = median(change_green),
    pct_brightened_green = mean(change_green > 0) * 100,
    pct_dimmed_green     = mean(change_green < 0) * 100,
    .groups = "drop"
  )

change_summary_imp <- all_data %>%
  group_by(season, buffer) %>%
  summarise(
    n = n(),
    mean_change_imp     = mean(change_imp),
    sd_change_imp       = sd(change_imp),
    median_change_imp   = median(change_imp),
    pct_brightened_imp   = mean(change_imp > 0) * 100,
    pct_dimmed_imp       = mean(change_imp < 0) * 100,
    .groups = "drop"
  )

cat("\n================================================================\n")
cat("TABLE 2: Change Variable Summaries by Season and Buffer\n")
cat("================================================================\n")
print(as.data.frame(change_summary_green), digits = 4, row.names = FALSE)
print(as.data.frame(change_summary_imp), digits = 4, row.names = FALSE)

write.csv(change_summary,
          file.path(out_dir, "table2_change_summaries.csv"),
          row.names = FALSE)

# ============================================================
# TABLE 3: Bright/Dim transition counts (Lambert Table 1 style)
# ============================================================

transition_green <- all_data %>%
  group_by(season, buffer, transition_green) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = transition_green, values_from = n, values_fill = 0)

transition_imp <- all_data %>%
  group_by(season, buffer, transition_imp) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = transition_imp, values_from = n, values_fill = 0)

cat("\n================================================================\n")
cat("TABLE 3a: Greenness Bright/Dim Transitions (2010 -> 2020)\n")
cat("================================================================\n")
print(as.data.frame(transition_green), row.names = FALSE)

cat("\n================================================================\n")
cat("TABLE 3b: Imperviousness Bright/Dim Transitions (2010 -> 2020)\n")
cat("================================================================\n")
print(as.data.frame(transition_imp), row.names = FALSE)

write.csv(transition_green,
          file.path(out_dir, "table3a_green_transitions.csv"),
          row.names = FALSE)
write.csv(transition_imp,
          file.path(out_dir, "table3b_imp_transitions.csv"),
          row.names = FALSE)


# ============================================================
# 6. STATISTICAL COMPARISONS
# ============================================================

# --- 6a. Two-way ANOVA ---
all_data$buffer_f <- factor(all_data$buffer)

cat("\n================================================================\n")
cat("ANOVA: change_green ~ season * buffer\n")
cat("================================================================\n")
aov_green <- aov(change_green ~ season * buffer_f, data = all_data)
print(summary(aov_green))

cat("\n================================================================\n")
cat("ANOVA: change_imp ~ season * buffer\n")
cat("================================================================\n")
aov_imp <- aov(change_imp ~ season * buffer_f, data = all_data)
print(summary(aov_imp))

# --- 6b. Tukey HSD for buffer within each season ---
cat("\n================================================================\n")
cat("Tukey HSD: Buffer comparisons within each season\n")
cat("================================================================\n")
for (s in unique(all_data$season)) {
  cat(sprintf("\n--- %s ---\n", s))
  sub <- all_data %>% filter(season == s)
  cat("  Greenness:\n")
  print(TukeyHSD(aov(change_green ~ buffer_f, data = sub)))
  cat("  Imperviousness:\n")
  print(TukeyHSD(aov(change_imp ~ buffer_f, data = sub)))
}

# --- 6c. Paired t-tests: season within each buffer ---
cat("\n================================================================\n")
cat("Paired t-tests: Season comparison within each buffer\n")
cat("================================================================\n")
for (b in sort(unique(all_data$buffer))) {
  cat(sprintf("\n--- Buffer: %d ---\n", b))
  sub_wide <- all_data %>%
    filter(buffer == b) %>%
    select(GEOID20, season, change_green, change_imp) %>%
    pivot_wider(names_from = season,
                values_from = c(change_green, change_imp),
                names_sep = "_") %>%
    drop_na()
  cat("  Greenness (Dormant vs Growing):\n")
  print(t.test(sub_wide$change_green_Dormant_Season,
               sub_wide$change_green_Growing_Season, paired = TRUE))
  cat("  Imperviousness (Dormant vs Growing):\n")
  print(t.test(sub_wide$change_imp_Dormant_Season,
               sub_wide$change_imp_Growing_Season, paired = TRUE))
}


# ============================================================
# 7. VISUALIZATIONS (modeled after Lambert et al. 2025)
# ============================================================

theme_lambert <- theme_minimal(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

buffer_labels <- c("500" = "500 m", "1000" = "1000 m", "2000" = "2000 m")

# -------------------------------------------------------
# FIGURE 1: 2010 Regression — NDVI vs log(HU_density)
#   Scatter + OLS line, faceted by season x buffer
#   (analogous to the baseline relationship in Lambert)
# -------------------------------------------------------

fig1 <- ggplot(all_data,
               aes(x = log_HU_density_2010, y = NDVI_2010)) +
  geom_point(alpha = 0.05, size = 0.3, color = "grey40") +
  geom_smooth(method = "lm", se = TRUE, color = "#2c7bb6", linewidth = 0.8) +
  facet_grid(season ~ buffer, labeller = labeller(buffer = buffer_labels)) +
  labs(
    title    = "2010 Baseline: Greenness vs. Housing Density",
    x = "log(Housing Unit Density per ha), 2010",
    y = "NDVI, 2010"
  ) +
  theme_lambert

ggsave(file.path(out_dir, "fig1_baseline_green_regression.png"),
       fig1, width = 10, height = 6, dpi = 300)

# -------------------------------------------------------
# FIGURE 2: 2010 Regression — Imperviousness vs log(HU_density)
# -------------------------------------------------------

fig2 <- ggplot(all_data,
               aes(x = log_HU_density_2010, y = PCT_BUILT_2010)) +
  geom_point(alpha = 0.05, size = 0.3, color = "grey40") +
  geom_smooth(method = "lm", se = TRUE, color = "#d7191c", linewidth = 0.8) +
  facet_grid(season ~ buffer, labeller = labeller(buffer = buffer_labels)) +
  labs(
    title    = "2010 Baseline: Imperviousness vs. Housing Density",
    x = "log(Housing Unit Density per ha), 2010",
    y = "Percent Built (Impervious), 2010"
  ) +
  theme_lambert

ggsave(file.path(out_dir, "fig2_baseline_imp_regression.png"),
       fig2, width = 10, height = 6, dpi = 300)


# -------------------------------------------------------
# FIGURE 3: Change in greenness index vs housing density change
#   (Lambert Figure 3a analog)
#   x = change in log(HU_density), y = change in greenness index
# -------------------------------------------------------

fig3 <- ggplot(all_data,
               aes(x = delta_log_HU, y = change_green)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_point(aes(color = direction_green), alpha = 0.15, size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7) +
  scale_color_manual(
    values = c("Brightened" = "#fdae61", "Dimmed" = "#878787", "No change" = "grey50"),
    name   = "Direction"
  ) +
  facet_grid(season ~ buffer, labeller = labeller(buffer = buffer_labels)) +
  labs(
    title    = "Change in Greenness Index vs. Housing Density Change",
    x = expression(Delta ~ "log(HU density/ha), 2010" %->% "2020"),
    y = expression(Delta ~ "Greenness Index (pseudo-residual change)")
  ) +
  theme_lambert

ggsave(file.path(out_dir, "fig3_change_green_vs_housing.png"),
       fig3, width = 10, height = 6, dpi = 300)

# -------------------------------------------------------
# FIGURE 4: Change in greenness index vs NDVI change
#   (Lambert Figure 3b analog)
# -------------------------------------------------------

fig4 <- ggplot(all_data,
               aes(x = delta_NDVI, y = change_green)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_point(aes(color = direction_green), alpha = 0.15, size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7) +
  scale_color_manual(
    values = c("Brightened" = "#fdae61", "Dimmed" = "#878787", "No change" = "grey50"),
    name   = "Direction"
  ) +
  facet_grid(season ~ buffer, labeller = labeller(buffer = buffer_labels)) +
  labs(
    title    = "Change in Greenness Index vs. Raw NDVI Change",
    x = expression(Delta ~ "NDVI, 2010" %->% "2020"),
    y = expression(Delta ~ "Greenness Index")
  ) +
  theme_lambert

ggsave(file.path(out_dir, "fig4_change_green_vs_ndvi.png"),
       fig4, width = 10, height = 6, dpi = 300)


# -------------------------------------------------------
# FIGURE 5: Change in imperviousness index vs housing density change
# -------------------------------------------------------

fig5 <- ggplot(all_data,
               aes(x = delta_log_HU, y = change_imp)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_point(aes(color = direction_imp), alpha = 0.15, size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7) +
  scale_color_manual(
    values = c("Brightened" = "#fdae61", "Dimmed" = "#878787", "No change" = "grey50"),
    name   = "Direction"
  ) +
  facet_grid(season ~ buffer, labeller = labeller(buffer = buffer_labels)) +
  labs(
    title    = "Change in Imperviousness Index vs. Housing Density Change",
    # subtitle = "Positive = less impervious than expected; Negative = more impervious",
    x = expression(Delta ~ "log(HU density/ha), 2010" %->% "2020"),
    y = expression(Delta ~ "Imperviousness Index (pseudo-residual change)")
  ) +
  theme_lambert

ggsave(file.path(out_dir, "fig5_change_imp_vs_housing.png"),
       fig5, width = 10, height = 6, dpi = 300)

# -------------------------------------------------------
# FIGURE 6: Change in imperviousness index vs raw impervious change
# -------------------------------------------------------

fig6 <- ggplot(all_data,
               aes(x = delta_IMP, y = change_imp)) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_point(aes(color = direction_imp), alpha = 0.15, size = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7) +
  scale_color_manual(
    values = c("Brightened" = "#fdae61", "Dimmed" = "#878787", "No change" = "grey50"),
    name   = "Direction"
  ) +
  facet_grid(season ~ buffer, labeller = labeller(buffer = buffer_labels)) +
  labs(
    title    = "Change in Imperviousness Index vs. Raw Impervious Change",
    # subtitle = "Index change as a function of impervious surface gain/loss",
    x = expression(Delta ~ "PCT_BUILT, 2010" %->% "2020"),
    y = expression(Delta ~ "Imperviousness Index")
  ) +
  theme_lambert

ggsave(file.path(out_dir, "fig6_change_imp_vs_pctbuilt.png"),
       fig6, width = 10, height = 6, dpi = 300)


# -------------------------------------------------------
# FIGURE 9: Distribution of change variables
#   Box + violin plots comparing season x buffer
# -------------------------------------------------------

fig9 <- all_data %>%
  select(season, buffer, buffer_f, change_green, change_imp) %>%
  pivot_longer(cols = c(change_green, change_imp),
               names_to = "variable", values_to = "value") %>%
  mutate(variable = recode(variable,
                           change_green = "Greenness Index Change",
                           change_imp   = "Imperviousness Index Change")) %>%
  ggplot(aes(x = buffer_f, y = value, fill = season)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.2, width = 0.6,
               position = position_dodge(0.7)) +
  facet_wrap(~ variable, scales = "free_y") +
  scale_fill_manual(values = c("Dormant_Season" = "#abd9e9",
                                "Growing_Season" = "#66c164"),
                    name = "Season",
                    labels = c("Dormant", "Growing")) +
  labs(
    title = "Distribution of Index Change by Buffer and Season",
    # title = "Distribution of Pseudo-Residual Change by Buffer and Season",
    subtitle = "Positive = brightening; Negative = dimming (relative to 2010 baseline)",
    x = "Buffer (m)",
    # y = "Index Change (2020 pseudo-residual - 2010 residual)"
    y = "Index Change"
  ) +
  theme_lambert

ggsave(file.path(out_dir, "fig9_boxplot_change_by_buffer_season.png"),
       fig9, width = 10, height = 5, dpi = 300)

# -------------------------------------------------------
# FIGURE 10: Mean change with 95% CI — dot-and-whisker
#   Clean comparison across all season x buffer combos
# -------------------------------------------------------

ci_data <- all_data %>%
  group_by(season, buffer) %>%
  summarise(
    mean_green = mean(change_green),
    se_green   = sd(change_green) / sqrt(n()),
    mean_imp   = mean(change_imp),
    se_imp     = sd(change_imp) / sqrt(n()),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(mean_green, se_green, mean_imp, se_imp),
    names_to = c(".value", "variable"),
    names_pattern = "(mean|se)_(green|imp)"
  ) %>%
  mutate(
    lo = mean - 1.96 * se,
    hi = mean + 1.96 * se,
    variable = recode(variable,
                      green = "Greenness",
                      imp   = "Imperviousness"),
    buffer_label = paste0(buffer, " m")
  )

fig10 <- ggplot(ci_data,
                aes(x = buffer_label, y = mean, color = season)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = lo, ymax = hi),
                  position = position_dodge(width = 0.4), size = 0.5) +
  facet_wrap(~ variable, scales = "free_y") +
  scale_color_manual(values = c("Dormant_Season" = "#2c7bb6",
                                 "Growing_Season" = "#1a9641"),
                     name = "Season",
                     labels = c("Dormant", "Growing")) +
  labs(
    # title = "Mean Change in Pseudo-Residual Index (95% CI)",
    title = "Mean Change in Index (95% CI)",
    # subtitle = "Comparing brightening/dimming across buffer sizes and seasons",
    x = "Buffer",
    y = "Mean Index Change"
  ) +
  theme_lambert

ggsave(file.path(out_dir, "fig10_mean_change_ci.png"),
       fig10, width = 8, height = 5, dpi = 300)


# -------------------------------------------------------
# FIGURE 11: Density ridgeline-style histograms of change
# -------------------------------------------------------

fig11 <- all_data %>%
  mutate(buffer_label = paste0(buffer, " m")) %>%
  ggplot(aes(x = change_green, fill = season)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ buffer_label, ncol = 1) +
  scale_fill_manual(values = c("Dormant_Season" = "#abd9e9",
                                "Growing_Season" = "#66c164"),
                    name = "Season",
                    labels = c("Dormant", "Growing")) +
  labs(
    title = "Distribution of Greenness Index Change",
    x = "Greenness Index Change (brightening +, dimming -)",
    y = "Density"
  ) +
  theme_lambert +
  coord_cartesian(xlim = quantile(all_data$change_green, c(0.005, 0.995)))

ggsave(file.path(out_dir, "fig11_density_green_change.png"),
       fig11, width = 8, height = 7, dpi = 300)

fig12 <- all_data %>%
  mutate(buffer_label = paste0(buffer, " m")) %>%
  ggplot(aes(x = change_imp, fill = season)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ buffer_label, ncol = 1) +
  scale_fill_manual(values = c("Dormant_Season" = "#abd9e9",
                                "Growing_Season" = "#66c164"),
                    name = "Season",
                    labels = c("Dormant", "Growing")) +
  labs(
    title = "Distribution of Imperviousness Index Change",
    x = "Imperviousness Index Change (brightening +, dimming -)",
    y = "Density"
  ) +
  theme_lambert +
  coord_cartesian(xlim = quantile(all_data$change_imp, c(0.005, 0.995)))

ggsave(file.path(out_dir, "fig12_density_imp_change.png"),
       fig12, width = 8, height = 7, dpi = 300)

# ============================================================
# SAVE FULL ANALYTICAL DATASET
# ============================================================

write.csv(all_data,
          file.path(out_dir, "block_pseudo_residual_analysis.csv"),
          row.names = FALSE)

cat("\n================================================================\n")
cat("All outputs saved to:", out_dir, "\n")
cat("  Tables: table1, table2, table3a, table3b CSVs\n")
cat("  Figures: fig1-fig12 PNGs\n")
cat("  Data: block_pseudo_residual_analysis.csv\n")
cat("================================================================\n")
cat("ANALYSIS COMPLETE\n")
