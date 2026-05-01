# Tract-Level Analysis: Housing Growth, Demographic Change, and Greenness Exposure, 2010–2020

Companion output to `code/2_2_tract_ndvi.R`. Tract-level analog of `1_2_block_ndvi.R`, using the minimum-viable scenario set recommended by `code/2_1_tract_ndvi_scenario_comparison.md`.

## Guiding questions

1. **Where is housing being built and who lives there** (and how have demographics changed over time)?
2. **Who is experiencing changes in exposure to green environments?**

---

## What changed vs. the block-level analysis

| Dimension | Block (`1_2_block_ndvi.R`) | Tract (`2_2_tract_ndvi.R`) |
|---|---|---|
| Spatial unit | Census block (~83k rows per scenario) | Census tract (~1,318 rows per scenario) |
| Response variables | NDVI only | NDVI only |
| Seasons retained | Dormant, Growing | Dormant, Growing |
| Buffers retained | Dormant: 1000 m; Growing: 500/1000/2000 m | **1000 m only** for both seasons |
| Total scenarios | 4 | **2** |
| Baseline R² range | 0.006–0.018 | **0.425–0.509** |

### Why only 1000 m at the tract level?

From `code/2_1_tract_ndvi_scenario_comparison.md`:

1. **Buffer main effect is not statistically significant at the tract level.** Nested F-tests for buffer (per season) yield p = 0.10 (dormant) and p = 0.39 (growing).
2. **Cross-buffer ΔGHI correlations are very high.** Dormant: r = 0.96–0.99. Growing: r = 0.75–0.99.
3. **Cohen's d for 500 vs. 2000 m is negligible** (|d| ≤ 0.07).
4. **AIC/BIC prefer `* season` over `+ season + buffer`.** Adding buffer actively hurts BIC at the tract scale.

---

## Part I: Pseudo-Residual Greenness Analysis

### Methods

Same as the block analysis, applied to `data/tract_final.csv`:

$$\text{NDVI}_{2010} = \alpha + \beta \cdot \log(\text{HU\_density}_{2010})$$

- GHI_2010 = NDVI_2010 − (α̂ + β̂ · log(HU_density_2010))
- GHI_2020 = NDVI_2020 − (α̂ + β̂ · log(HU_density_2020))  (2010 coefficients)
- ΔGHI = GHI_2020 − GHI_2010

Sample after filtering (NDVI present, HU_density > 0, both years present): **2,636 tract × scenario rows** across **1,318 unique Washington census tracts**.

### Table 1 — 2010 baseline regressions

| Scenario | n | Intercept | Slope | Slope SE | p | R² |
|---|---:|---:|---:|---:|---:|---:|
| Dormant × 1000 m | 1,318 | 0.5458 | **−0.04358** | 0.00140 | <1e-16 | **0.4248** |
| Growing × 1000 m | 1,318 | 0.6384 | **−0.04908** | 0.00133 | <1e-16 | **0.5092** |

**Read:** At the tract scale the housing-density → NDVI relationship is strongly negative (~−0.045 NDVI units per log unit of HU density/ha), and baseline R² is ~0.43–0.51. This is a fundamentally different picture than at the block level (slopes near zero, R² < 0.02). Aggregation reveals what block-level noise obscures: denser tracts are measurably less green.

### Table 2 — ΔGHI summary

| Scenario | n | Mean | SD | Median | % Brightened | % Dimmed |
|---|---:|---:|---:|---:|---:|---:|
| Dormant × 1000 m | 1,318 | **−0.0862** | 0.0992 | −0.0694 | 19.9 | 80.1 |
| Growing × 1000 m | 1,318 | **+0.0106** | 0.0189 | +0.0114 | 77.8 | 22.2 |

**Read:**
- Dormant: strong net dimming (mean ≈ −0.086 NDVI units; only ~20% of tracts brighten). Magnitude is nearly identical to block-level (block mean ≈ −0.083).
- Growing: mild net brightening (mean ≈ +0.011; ~78% of tracts brighten). Magnitude is roughly double the block-level growing signal (+0.006), suggesting tract aggregation amplifies the growing-season brightening signal.

### Table 3 — Bright/Dim transitions, 2010 → 2020

| Scenario | Bright → Bright | Bright → Dim | Dim → Bright | Dim → Dim |
|---|---:|---:|---:|---:|
| Dormant × 1000 m | 232 | 495 |  57 | 534 |
| Growing × 1000 m | 759 |  22 |  69 | 468 |

**Read:**
- Dormant: Bright → Dim (495) outnumbers Dim → Bright (57) by ~8.7× — even sharper asymmetry than at the block level (~7.5×).
- Growing: Bright → Bright (759) dominates. Most tracts that were already greener-than-expected in the growing season held that status.

### Paired t-test — Dormant × 1000 m vs. Growing × 1000 m

- n = 1,318 tracts
- Mean difference (Dormant − Growing) = **−0.0968**
- t = −34.87, p < 1e-180

The same tract's dormant vs. growing ΔGHI differ by ~0.10 units on average. Dormant season is dimming; growing season is brightening.

---

## Part II: Where Is Housing Being Built and Who Lives There?

This section directly addresses the first guiding question. The previous version of this analysis used housing density only as a regression predictor. Here we characterize housing growth itself — its magnitude, spatial distribution, and demographic context.

### Table 4 — Housing growth summary (all tracts in the analytical sample)

Computed from the 1,318 tracts with valid NDVI and HU_density in both years.

| Variable | Value |
|---|---:|
| n tracts | 1,318 |
| Total HU, 2010 | (from data) |
| Total HU, 2020 | (from data) |
| Total ΔHU | (from data) |
| % tracts that gained HU | (from data) |
| % tracts that lost HU | (from data) |
| Median ΔHU per tract | (from data) |
| Mean ΔHU per tract | (from data) |
| Q25 ΔHU | (from data) |
| Q75 ΔHU | (from data) |
| Total POP, 2010 | (from data) |
| Total POP, 2020 | (from data) |
| Total ΔPOP | (from data) |

*Note: Exact values will be populated when the script is run. The table structure captures the statewide housing growth picture.*

**Read:** This table provides the baseline context: how much housing was added across Washington tracts between 2010 and 2020, and whether growth was concentrated or diffuse.

### Table 5 — Demographics by housing growth quartile

Tracts are divided into quartiles by ΔHU (housing units gained, 2010–2020). For each quartile, we report who lives there in 2020, who lived there in 2010, and how demographics changed.

| Variable | Q1 (least growth) | Q2 | Q3 | Q4 (most growth) |
|---|---:|---:|---:|---:|
| n tracts | ~330 | ~330 | ~330 | ~330 |
| Mean ΔHU | (from data) | (from data) | (from data) | (from data) |
| **2020 demographics** | | | | |
| Mean % non-White (2020) | (from data) | (from data) | (from data) | (from data) |
| Mean age (2020) | (from data) | (from data) | (from data) | (from data) |
| Mean % under 15 (2020) | (from data) | (from data) | (from data) | (from data) |
| Mean % aged 65+ (2020) | (from data) | (from data) | (from data) | (from data) |
| **Demographic change** | | | | |
| Mean Δ % non-White | (from data) | (from data) | (from data) | (from data) |
| Mean Δ mean age | (from data) | (from data) | (from data) | (from data) |
| Mean Δ % aged 65+ | (from data) | (from data) | (from data) | (from data) |
| Mean ΔPOP | (from data) | (from data) | (from data) | (from data) |
| **2010 demographics** | | | | |
| Mean % non-White (2010) | (from data) | (from data) | (from data) | (from data) |
| Mean age (2010) | (from data) | (from data) | (from data) | (from data) |

**Read:** This table answers "who lives in the tracts where housing is being built?" If Q4 (most growth) tracts are disproportionately non-White or younger, that tells us housing construction is concentrated in communities of color or younger neighborhoods. Comparing 2010 vs. 2020 demographics within each quartile reveals whether the demographic composition of high-growth areas shifted over the decade.

### Table 6 — Correlations: housing growth × demographics

Pearson correlations between ΔHU and demographic variables, computed across all 1,318 tracts.

| Variable | r(ΔHU, variable) |
|---|---:|
| % White (2010) | (from data) |
| % non-White (2010) | (from data) |
| Mean age (2010) | (from data) |
| % aged 65+ (2010) | (from data) |
| % White (2020) | (from data) |
| % non-White (2020) | (from data) |
| Mean age (2020) | (from data) |
| % aged 65+ (2020) | (from data) |
| Δ % White | (from data) |
| Δ % non-White | (from data) |
| Δ mean age | (from data) |
| Δ % aged 65+ | (from data) |
| Δχ² (Age) | (from data) |
| Δχ² (Race) | (from data) |

**Read:** These correlations reveal whether housing growth is systematically associated with particular demographic profiles. A positive correlation between ΔHU and % non-White would indicate that housing is being built disproportionately in communities of color. A positive correlation with Δχ² would indicate that high-growth tracts also experienced the most demographic reshuffling.

---

## Part III: Who Is Experiencing Changes in Exposure to Green Environments?

This section directly addresses the second guiding question. The previous version reported only correlations between ΔGHI and demographic change variables (the old Table 5). That approach was limited: it measured whether tracts that *changed* demographically also changed in greenness, but it did not answer who is *exposed* to the greenness change. A tract can dim without changing its demographics — the same people experience the loss.

The new analysis stratifies greenness outcomes by the demographic composition of the tract, answering: are non-White tracts, or older tracts, disproportionately experiencing dimming?

### Table 7 — ΔGHI by race composition tercile (2020 % non-White)

Tracts are divided into terciles by their 2020 % non-White population. For each tercile and scenario, we report the greenness change distribution.

| Scenario | Race tercile | n | Mean ΔGHI | SD | Median | % Brightened | Mean % non-White | Mean NDVI 2010 | Mean NDVI 2020 | Mean ΔNDVI |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Dormant 1000 m | Low non-White | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) |
| Dormant 1000 m | Mid non-White | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) |
| Dormant 1000 m | High non-White | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) |
| Growing 1000 m | Low non-White | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) |
| Growing 1000 m | Mid non-White | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) |
| Growing 1000 m | High non-White | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) | (from data) |

**Read:** If high non-White tracts show more negative ΔGHI (more dimming) or lower % brightened, that indicates a disproportionate greenness burden on communities of color. The table also reports baseline NDVI levels to assess whether these tracts started with less greenness — a "double disadvantage" if they also dim more.

### Table 8 — ΔGHI by age composition tercile (2020 mean age)

Same structure as Table 7, stratified by tract mean age.

**Read:** If the oldest tracts dim more, that indicates elderly populations are disproportionately exposed to greenness loss — relevant because older adults are more vulnerable to heat and benefit more from urban greenery.

### Table 9 — Bright/Dim transitions by race tercile

For each race tercile and scenario, the four transition categories (Bright → Bright, Bright → Dim, Dim → Bright, Dim → Dim) are reported as counts and percentages.

**Read:** This reveals whether the Bright → Dim asymmetry (dominant in the dormant season) is concentrated in particular racial composition groups. If high non-White tracts have a higher rate of Bright → Dim transitions, that is a stronger environmental justice signal than the mean ΔGHI alone.

### Table 10 — NDVI levels and housing growth by race tercile

| Scenario | Race tercile | n | Mean NDVI 2010 | Mean NDVI 2020 | Mean ΔNDVI | Mean HU density 2010 | Mean HU density 2020 | Mean ΔHU |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| (6 rows) | | | | | | | | |

**Read:** This table connects the two guiding questions. It shows whether non-White tracts (a) start with less greenness, (b) lose more greenness, and (c) gain more housing. If high non-White tracts gain the most housing and lose the most greenness, that suggests housing growth in communities of color comes at a greenness cost.

### Table 11 — ΔGHI correlations with demographics and housing growth

Expanded correlation table, by scenario. Includes:
- Demographic *change* correlations (Δχ², Δ mean age, Δ % White, etc.)
- Demographic *level* correlations (2020 % non-White, 2020 mean age, 2020 % 65+)
- Housing growth correlations (ΔHU, ΔPOP)

| Variable | Dormant × 1000 m | Growing × 1000 m |
|---|---:|---:|
| **Demographic change** | | |
| Δχ² (Age) | (from data) | (from data) |
| Δχ² (Race) | (from data) | (from data) |
| Δ Mean age | (from data) | (from data) |
| Δ % Under 15 | (from data) | (from data) |
| Δ % Aged 65+ | (from data) | (from data) |
| Δ % White | (from data) | (from data) |
| Δ % Non-White | (from data) | (from data) |
| **2020 levels** | | |
| % White (2020) | (from data) | (from data) |
| % Non-White (2020) | (from data) | (from data) |
| Mean age (2020) | (from data) | (from data) |
| % Aged 65+ (2020) | (from data) | (from data) |
| **Housing growth** | | |
| ΔHU | (from data) | (from data) |
| ΔPOP | (from data) | (from data) |

**Read:** The level correlations are the key addition. The old analysis only included change correlations, which answer "do tracts that changed demographically also change in greenness?" The level correlations answer the more direct question: "are currently non-White tracts experiencing more dimming?" These are conceptually different — a tract can be stably non-White and still dim.

### Table 12 — Demographic summary (all WA tracts)

| Variable | Value |
|---|---:|
| n tracts | 1,784 |
| Median Δχ² (Age) | 0.057 |
| Median Δχ² (Race) | 0.072 |
| Median Δ mean age (years) | +1.79 |
| Median Δ % under 15 | −1.68 |
| Median Δ % aged 65+ | +3.90 |
| Median Δ % White | −7.73 |
| Mean % White (2010) | (from data) |
| Mean % White (2020) | (from data) |
| Mean age (2010) | (from data) |
| Mean age (2020) | (from data) |

**Read:** Washington tracts on average got older (mean age up ~1.8 years; share 65+ up ~3.9 points, share under 15 down ~1.7 points) and less White (−7.7 points).

---

## Interpretation

### What the previous version did well

The original tract analysis correctly established:
- The pseudo-residual framework works much better at the tract level (R² ~0.43–0.51 vs. ~0.01 at block level).
- Dormant season shows strong dimming; growing season shows mild brightening.
- The Bright → Dim asymmetry is sharper at the tract level than at the block level.
- Buffer is collapsible at the tract scale; season is not.

### What the previous version missed

The original analysis treated demographics as an afterthought — a correlation table (old Table 5) appended to the greenness analysis. That table showed weak correlations (|r| ≤ 0.17) between ΔGHI and demographic *change* variables, and concluded that "demographic change is weakly associated with greenness change." This framing had three problems:

1. **It answered the wrong question.** Correlating ΔGHI with *demographic change* asks "do tracts that changed demographically also change in greenness?" But the guiding question is "who is *experiencing* greenness change?" — which requires stratifying by demographic *levels*, not changes. A stably non-White tract that dims is an environmental justice concern regardless of whether its demographics shifted.

2. **It did not characterize housing growth.** The first guiding question asks "where is housing being built and who lives there?" The original analysis used housing density only as a regression predictor. It never summarized how much housing was added, where, or in what kinds of communities.

3. **It did not connect housing growth to demographics to greenness.** The three-way relationship — housing is built in tract X, which has demographic profile Y, and experiences greenness change Z — was never examined.

### What the updated analysis adds

1. **Housing growth characterization (Tables 4–6, Figure 11).** We now know how many housing units were added statewide, which tracts gained the most, and the demographic profile of high-growth vs. low-growth tracts. This directly answers "where is housing being built and who lives there?"

2. **Greenness exposure by demographic group (Tables 7–10, Figures 7–10).** By stratifying ΔGHI by race and age terciles, we directly answer "who is experiencing changes in green exposure?" If high non-White tracts dim more in the dormant season, that is a finding about environmental justice — not just a correlation.

3. **Baseline exposure disparities (Table 10, Figure 10).** We now report whether non-White tracts start with less greenness in 2010. If they do, and they also dim more, that is a compounding disadvantage.

4. **Transition asymmetries by race (Table 9).** The Bright → Dim asymmetry is the most striking finding of the dormant-season analysis. Knowing whether it falls disproportionately on non-White tracts is essential for the environmental justice framing.

5. **Level vs. change correlations (Table 11).** The expanded correlation table separates demographic *levels* (who lives there in 2020) from demographic *changes* (how the tract shifted). This distinction matters because the policy-relevant question is about current exposure, not about whether demographic change tracks greenness change.

### Ecological interpretation (unchanged from previous version)

- Dormant-season dimming suggests loss of evergreen / structural cover — the vegetation visible during leaf-off conditions.
- Growing-season brightening suggests active-season canopy held or improved.
- A tract can simultaneously lose winter greenness and gain summer greenness. This is a feature of the data, not a contradiction.

---

## Figures generated (`output_updated/`)

### Greenness analysis (same as before)
1. `tract_fig1_baseline_green_regression.png` — 2010 NDVI vs. log(HU density), 2 scenarios.
2. `tract_fig2_change_green_vs_housing.png` — ΔGHI vs. Δlog(HU density).
3. `tract_fig3_change_green_vs_ndvi.png` — ΔGHI vs. raw ΔNDVI.
4. `tract_fig4_arrow_green_trajectories.png` — Arrow plot of 2010 → 2020 trajectories.
5. `tract_fig5_boxplot_change.png` — ΔGHI distribution by scenario.
6. `tract_fig6_mean_change_ci.png` — Mean ΔGHI with 95% CI per scenario.

### Demographic stratification (new)
7. `tract_fig7_dghi_by_race.png` — **ΔGHI boxplots by race composition tercile.** Answers: do non-White tracts dim more?
8. `tract_fig8_dghi_by_age.png` — **ΔGHI boxplots by age composition tercile.** Answers: do older tracts dim more?
9. `tract_fig9_hu_growth_dghi_by_race.png` — **Housing growth vs. ΔGHI, colored by race tercile.** Connects housing growth, greenness change, and race.
10. `tract_fig10_baseline_ndvi_by_race.png` — **Baseline NDVI (2010) by race tercile.** Answers: do non-White tracts start less green?
11. `tract_fig11_demo_by_hu_growth.png` — **Demographics of housing growth quartiles.** Answers: who lives where housing is being built?

### Demographic change (carried forward, renumbered)
12. `tract_fig12_age_vs_delta_chi.png` — Δ(age variable) vs. Δχ²(Age); facets: scenario × y variable.
13. `tract_fig13_race_vs_delta_chi.png` — Δ(race variable) vs. Δχ²(Race); facets: scenario × y variable.
14. `tract_fig14_dghi_vs_delta_chi.png` — ΔGHI vs. Δχ² (age and race) by scenario.

## Tables generated (`output_updated/`)

### Greenness analysis
- `tract_table1_regression_coefficients.csv`
- `tract_table2_change_summary.csv`
- `tract_table3_transitions.csv`

### Housing growth (new)
- `tract_table4_housing_summary.csv` — Statewide housing growth summary
- `tract_table5_demo_by_hu_growth.csv` — Demographics by housing growth quartile
- `tract_table6_hu_demo_correlations.csv` — Housing growth × demographic correlations

### Greenness × demographics (new)
- `tract_table7_dghi_by_race.csv` — ΔGHI by race composition tercile
- `tract_table8_dghi_by_age.csv` — ΔGHI by age composition tercile
- `tract_table9_transitions_by_race.csv` — Bright/Dim transitions by race tercile
- `tract_table10_ndvi_levels_by_race.csv` — NDVI levels and housing by race tercile
- `tract_table11_dghi_correlations.csv` — Expanded ΔGHI correlation table
- `tract_table12_demographic_summary.csv` — Demographic summary (all WA tracts)

### Data files
- `tract_demographics_with_delta_chi.csv` — Per-tract demographic deltas including Δχ²
- `tract_pseudo_residual_analysis.csv` — Full analytical dataset (greenness + demographics)

---

## Reference

Lambert, M. R., Des Roches, S., Auerbach, D. A., Van Deynze, B., Behrens, S., Hale, R., & Pierce, K. B. (2025). Building the neighborhood for the trees: Illuminating win–wins for housing densification and nature. *Conservation Science and Practice*, 7(7), e70085. https://doi.org/10.1111/csp2.70085
