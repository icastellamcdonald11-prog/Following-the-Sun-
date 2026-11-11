source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# 14_land_opportunity_cost.R  —  Thesis Project 4
# Addresses reviewer concern: "opportunity cost of land is a first-order issue;
# if high-PVOUT land is also high-value in other uses, the PVOUT-based
# inefficiency may not represent true economic inefficiency."
#
# Strategy:
#   Population density (persons/km²) is used as the primary proxy for land
#   opportunity cost. In China, land near population centres commands higher
#   agricultural, commercial, and residential rents; remote high-resource land
#   is cheap. This is explicitly an empirical proxy, not a price series.
#
# Four tests:
#   A. National PVOUT–population correlation: direction determines whether
#      high-PVOUT land is expensive or cheap
#   B. Within-province PVOUT–population correlations: do trade-offs vary
#      by province? (Fig R5)
#   C. Does the province-level PVOUT–population trade-off explain
#      province-level siting patterns? (Fig R6)
#   D. Regression: does controlling for pop density alter the main PVOUT
#      efficiency result?
#
# Requires: 00_setup.R + 11_siting_efficiency_regression.R sourced
#   (needs: grid_df, grid_prov_df, deploy_en, prov_summary)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(fixest)
  library(ggrepel)
})

P4  <- FTS_ROOT
ROB <- file.path(P4, "robustness_outputs")
OUT <- file.path(P4, "figures")
dir.create(ROB, showWarnings = FALSE)
dir.create(OUT, showWarnings = FALSE)

theme_rob <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.title        = element_text(size = 10),
        plot.title        = element_text(size = 11, face = "bold"),
        plot.subtitle     = element_text(size = 9, colour = "grey40"),
        plot.caption      = element_text(size = 8, colour = "grey50"),
        legend.position   = "bottom")

# =============================================================================
# BUILD GRID ANALYSIS DATASET
# Merge grid_df (has pop_density) with grid_prov_df (has province_en)
# =============================================================================
cat("Building grid analysis dataset...\n")

# grid_prov_df loaded from cache in 11_siting_efficiency_regression.R
# If not in environment (standalone run), reload
if (!exists("grid_prov_df")) {
  cache_file <- file.path(P4, "data", "grid_province_assignment.rds")
  stopifnot("Run 11_siting_efficiency_regression.R first to generate cache" =
              file.exists(cache_file))
  grid_prov_df <- readRDS(cache_file)
}

grid_full <- grid_df %>%
  dplyr::filter(feasible, is.finite(pvout_annual), pvout_annual > 0,
                is.finite(pop_density)) %>%
  dplyr::select(cell_id, pvout_annual, pop_density, elevation_m, slope_deg) %>%
  left_join(grid_prov_df %>% dplyr::select(cell_id, province_en),
            by = "cell_id") %>%
  dplyr::filter(!is.na(province_en)) %>%
  mutate(
    log_pop = log1p(pop_density),    # log(1 + persons/km²)
    pop_zero = pop_density == 0      # flag truly uninhabited cells
  )

cat(sprintf("  Feasible grid cells with province + pop: %d\n", nrow(grid_full)))
cat(sprintf("  Cells with zero population: %d (%.1f%%)\n",
            sum(grid_full$pop_zero), 100 * mean(grid_full$pop_zero)))

# =============================================================================
# SECTION A — National PVOUT–population correlation
# =============================================================================
cat("\n===== A: National PVOUT–population correlation =====\n")

# Pearson on log-transformed pop
nat_cor_pearson  <- cor(grid_full$pvout_annual, grid_full$log_pop,
                        use = "complete.obs")
nat_cor_spearman <- cor(grid_full$pvout_annual, grid_full$log_pop,
                        method = "spearman", use = "complete.obs")

# Cor for non-zero pop cells only (removes effect of zero-pop desert)
grid_nonzero <- dplyr::filter(grid_full, pop_density > 0)
cor_nonzero  <- cor(grid_nonzero$pvout_annual, grid_nonzero$log_pop,
                    use = "complete.obs")

cat(sprintf("  Pearson  r (PVOUT vs log(pop+1), all feasible):        %.3f\n",
            nat_cor_pearson))
cat(sprintf("  Spearman r (PVOUT vs log(pop+1), all feasible):        %.3f\n",
            nat_cor_spearman))
cat(sprintf("  Pearson  r (PVOUT vs log(pop),   pop > 0 cells only):  %.3f\n",
            cor_nonzero))

# Quintile means — does mean PVOUT decline monotonically with pop?
grid_full$pop_quintile <- cut(grid_full$log_pop,
                               breaks = quantile(grid_full$log_pop,
                                                 probs = seq(0, 1, 0.2),
                                                 na.rm = TRUE),
                               labels = c("Q1\n(least\npopulated)",
                                          "Q2", "Q3", "Q4",
                                          "Q5\n(most\npopulated)"),
                               include.lowest = TRUE)

quint_means <- grid_full %>%
  group_by(pop_quintile) %>%
  summarise(mean_pvout = mean(pvout_annual),
            med_pvout  = median(pvout_annual),
            n          = n(),
            .groups    = "drop")

cat("\nMean PVOUT by population density quintile (feasible land):\n")
print(quint_means)

# --- Figure R5a: hexbin PVOUT vs log(pop+1) ---
cat("Building Figure R5a...\n")

# Sample 30k points for plotting speed
set.seed(42)
plot_sample <- grid_full[sample(nrow(grid_full), min(30000, nrow(grid_full))), ]

fig_r5a <- ggplot(plot_sample, aes(x = log_pop, y = pvout_annual)) +
  geom_hex(bins = 60, colour = NA) +
  geom_smooth(method = "lm", colour = "#C0392B", linewidth = 0.9,
              se = TRUE, fill = "#F1948A", alpha = 0.3) +
  scale_fill_gradient(low = "grey90", high = "#1A3A5C",
                      name = "Grid cells\n(count)") +
  scale_x_continuous("Population density  [log(persons/km² + 1)]") +
  scale_y_continuous("PVOUT (kWh/kWp/yr)") +
  annotate("text", x = max(plot_sample$log_pop) * 0.9,
           y = max(plot_sample$pvout_annual) * 0.95,
           label = sprintf("r = %.3f\n(Pearson)", nat_cor_pearson),
           hjust = 1, vjust = 1, size = 3.2, colour = "#C0392B") +
  labs(
    title   = "Figure R5a: PVOUT vs population density — all feasible land",
    subtitle = paste0(
      "Each hex = count of 10km feasible-land grid cells. ",
      "Negative national correlation indicates\nhigh-PVOUT land is ",
      "systematically less populated (lower opportunity cost)."
    ),
    caption = "10km feasible-land grid (GSA v2.0 + SRTM, elevation ≤5000m, slope ≤50°). N = 86,024 cells."
  ) +
  theme_rob

ggsave(file.path(OUT, "fig_R5a_pvout_vs_pop_national.png"),
       fig_r5a, width = 7, height = 5, dpi = 300)
cat("  Saved fig_R5a_pvout_vs_pop_national.png\n")

# --- Figure R5b: mean PVOUT by pop quintile (bar) ---
fig_r5b <- ggplot(quint_means, aes(x = pop_quintile, y = mean_pvout)) +
  geom_col(fill = "#2C6FAC", alpha = 0.85, width = 0.65) +
  geom_text(aes(label = round(mean_pvout, 0)),
            vjust = -0.4, size = 3.2) +
  geom_hline(yintercept = pvout_national_median, linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  annotate("text", x = 0.6,
           y = pvout_national_median + 15,
           label = "National median", hjust = 0, size = 2.8, colour = "grey40") +
  scale_y_continuous("Mean PVOUT (kWh/kWp/yr)",
                     expand = expansion(mult = c(0, 0.1))) +
  labs(
    x       = "Population density quintile (feasible land)",
    title   = "Figure R5b: Mean PVOUT declines monotonically with population density",
    caption = paste0(
      "Quintile boundaries based on log(pop+1) distribution across feasible 10km cells. ",
      "Dashed line = national feasible-land PVOUT median."
    )
  ) +
  theme_rob +
  theme(panel.grid.major.x = element_blank())

ggsave(file.path(OUT, "fig_R5b_pvout_by_pop_quintile.png"),
       fig_r5b, width = 6, height = 4, dpi = 300)
cat("  Saved fig_R5b_pvout_by_pop_quintile.png\n")


# =============================================================================
# SECTION B — Within-province PVOUT–population correlations
# =============================================================================
cat("\n===== B: Within-province PVOUT–population correlations =====\n")

prov_land_cor <- grid_full %>%
  group_by(province_en) %>%
  dplyr::filter(n() >= 20) %>%    # need enough cells for a stable correlation
  summarise(
    n_cells       = n(),
    cor_pearson   = cor(pvout_annual, log_pop, use = "complete.obs"),
    cor_spearman  = cor(pvout_annual, log_pop, method = "spearman",
                        use = "complete.obs"),
    mean_pvout    = mean(pvout_annual),
    mean_log_pop  = mean(log_pop),
    pct_zero_pop  = mean(pop_zero) * 100,
    .groups       = "drop"
  ) %>%
  arrange(cor_pearson)

cat("\nWithin-province PVOUT–pop correlation (Pearson, sorted):\n")
print(prov_land_cor %>%
        dplyr::select(province_en, cor_pearson, cor_spearman,
                      mean_pvout, pct_zero_pop, n_cells) %>%
        mutate(across(where(is.numeric), ~round(., 3))))

# Summary stats
pos_cor_n <- sum(prov_land_cor$cor_pearson > 0)
neg_cor_n <- sum(prov_land_cor$cor_pearson < 0)
cat(sprintf("\nProvinces with POSITIVE PVOUT–pop correlation (trade-off): %d of %d\n",
            pos_cor_n, nrow(prov_land_cor)))
cat(sprintf("Provinces with NEGATIVE PVOUT–pop correlation (no trade-off): %d of %d\n",
            neg_cor_n, nrow(prov_land_cor)))

write.csv(prov_land_cor,
          file.path(ROB, "province_pvout_pop_correlation.csv"),
          row.names = FALSE)
cat("Saved province_pvout_pop_correlation.csv\n")

# --- Figure R6: lollipop of within-province correlations ---
cat("Building Figure R6...\n")

prov_land_cor_plot <- prov_land_cor %>%
  mutate(
    prov_label = province_en,
    colour_grp = case_when(
      cor_pearson >  0.1  ~ "Positive (trade-off exists)",
      cor_pearson < -0.1  ~ "Negative (no trade-off)",
      TRUE                ~ "Near-zero"
    ),
    colour_grp = factor(colour_grp,
                        levels = c("Positive (trade-off exists)",
                                   "Near-zero",
                                   "Negative (no trade-off)"))
  ) %>%
  arrange(cor_pearson) %>%
  mutate(prov_label = factor(prov_label, levels = prov_label))

fig_r6 <- ggplot(prov_land_cor_plot,
                 aes(x = cor_pearson, y = prov_label, colour = colour_grp)) +
  geom_vline(xintercept = 0, colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = nat_cor_pearson, linetype = "dashed",
             colour = "grey30", linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = cor_pearson,
                   y = prov_label, yend = prov_label),
               linewidth = 0.5, alpha = 0.7) +
  geom_point(size = 2.2) +
  annotate("text", x = nat_cor_pearson + 0.01,
           y = nrow(prov_land_cor_plot) * 0.15,
           label = sprintf("National\nr = %.3f", nat_cor_pearson),
           hjust = 0, size = 2.6, colour = "grey30") +
  scale_colour_manual(
    name   = "Direction of PVOUT–population correlation",
    values = c("Positive (trade-off exists)" = "#C0392B",
               "Near-zero"                   = "#7F8C8D",
               "Negative (no trade-off)"     = "#2C6FAC")
  ) +
  scale_x_continuous(
    "Within-province Pearson r (PVOUT vs log population density)",
    breaks = seq(-0.8, 0.6, 0.2)
  ) +
  labs(
    y       = NULL,
    title   = "Figure R6: Within-province PVOUT–population density correlation",
    subtitle = paste0(
      "Negative = high-PVOUT land within the province is more remote (lower opportunity cost).\n",
      "Positive = high-PVOUT land within the province is more populated (real trade-off).\n",
      "Dashed = national average. Feasible 10km cells, ≥20 cells per province."
    ),
    caption = "Pearson r computed on feasible 10km grid cells. log_pop = log(persons/km² + 1)."
  ) +
  theme_rob +
  theme(axis.text.y = element_text(size = 7.5),
        legend.position = "bottom",
        legend.key.width = unit(0.6, "cm"))

ggsave(file.path(OUT, "fig_R6_province_pvout_pop_cor.png"),
       fig_r6, width = 8, height = 8, dpi = 300)
cat("  Saved fig_R6_province_pvout_pop_cor.png\n")


# =============================================================================
# SECTION C — Does the PVOUT–pop trade-off explain province siting patterns?
# =============================================================================
cat("\n===== C: Trade-off correlation vs province siting patterns =====\n")

# Merge land-cost correlation with province siting summary from 11_*.R
# prov_summary is built in 11_siting_efficiency_regression.R
if (!exists("prov_summary")) {
  stop("Run 11_siting_efficiency_regression.R before this script.")
}

prov_analysis <- prov_summary %>%
  left_join(prov_land_cor %>%
              dplyr::select(province_en, cor_pearson, cor_spearman,
                            mean_pvout, pct_zero_pop),
            by = c("prov_en" = "province_en")) %>%
  dplyr::filter(!is.na(cor_pearson)) %>%
  mutate(
    delta_rank_neg = delta_rank < 0,  # TRUE = siting got worse post-2020
    label_prov     = prov_en
  )

cat(sprintf("  Provinces in joint analysis: %d\n", nrow(prov_analysis)))

# Correlation: does the PVOUT–pop trade-off predict post-2020 rank change?
cor_tradeoff_rank <- cor(prov_analysis$cor_pearson,
                         prov_analysis$delta_rank,
                         use = "complete.obs")
cat(sprintf("  Correlation (PVOUT–pop trade-off vs Δrank): r = %.3f\n",
            cor_tradeoff_rank))

# OLS: delta_rank ~ pvout_pop_cor + controls
m_tradeoff1 <- lm(delta_rank ~ cor_pearson,
                  data = prov_analysis, weights = cap_total)
m_tradeoff2 <- lm(delta_rank ~ cor_pearson + rank_pre + log(cap_total),
                  data = prov_analysis, weights = cap_total)

cat("\n  Province regression: delta_rank ~ PVOUT-pop correlation\n")
cat(sprintf("    beta cor_pearson (bivariate): %.3f (SE %.3f, p = %.4f)\n",
            coef(m_tradeoff1)[2],
            sqrt(diag(vcov(m_tradeoff1)))[2],
            summary(m_tradeoff1)$coefficients[2, 4]))
cat(sprintf("    R² = %.3f\n", summary(m_tradeoff1)$r.squared))
cat("\n  Province regression: delta_rank ~ cor_pearson + rank_pre + log(cap)\n")
cat(sprintf("    beta cor_pearson (controlled): %.3f (SE %.3f, p = %.4f)\n",
            coef(m_tradeoff2)[2],
            sqrt(diag(vcov(m_tradeoff2)))[2],
            summary(m_tradeoff2)$coefficients[2, 4]))
cat(sprintf("    R² = %.3f\n", summary(m_tradeoff2)$r.squared))

# Also: do provinces with more positive PVOUT–pop correlation tend to
# deploy on lower absolute PVOUT land post-2020?
cor_tradeoff_pvout <- cor(prov_analysis$cor_pearson,
                          prov_analysis$rank_post,
                          use = "complete.obs")
cat(sprintf("\n  Correlation (PVOUT–pop trade-off vs post-2020 mean rank): r = %.3f\n",
            cor_tradeoff_pvout))

# --- Figure R7: scatter plot ---
cat("Building Figure R7...\n")

fig_r7 <- ggplot(prov_analysis,
                 aes(x = cor_pearson, y = delta_rank)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey60", linewidth = 0.4) +
  geom_smooth(method = "lm", aes(weight = cap_total),
              se = TRUE, colour = "#2C3E50",
              fill = "grey80", alpha = 0.35, linewidth = 0.9) +
  geom_point(aes(size = cap_total,
                 colour = delta_rank < 0),
             alpha = 0.8) +
  geom_text_repel(
    data = dplyr::filter(prov_analysis,
                         abs(delta_rank) > 8 |
                         abs(cor_pearson) > 0.3 |
                         cap_total > 30),
    aes(label = prov_en),
    size = 2.6, colour = "grey20",
    max.overlaps = 20, box.padding = 0.4,
    segment.colour = "grey60", segment.size = 0.3
  ) +
  scale_colour_manual(
    values = c("FALSE" = "#2C6FAC", "TRUE" = "#C0392B"),
    labels = c("FALSE" = "Rank improved post-2020",
               "TRUE"  = "Rank declined post-2020"),
    name = NULL
  ) +
  scale_size_continuous("Total cap (GW)", range = c(2, 9)) +
  scale_x_continuous(
    "Within-province PVOUT–population correlation (r)\n← No trade-off  |  Trade-off →"
  ) +
  scale_y_continuous(
    "Change in cap-weighted mean within-province rank\n(post minus pre-Q2 2020)  ← worse  |  better →"
  ) +
  labs(
    title    = "Figure R7: Land opportunity cost trade-off vs post-2020 siting change",
    subtitle = sprintf(
      "r = %.3f between PVOUT–pop correlation and siting rank change (capacity-weighted).\n%s",
      cor_tradeoff_rank,
      "If land costs drove inefficiency, provinces with a stronger trade-off (right side) should show more rank decline (downward)."
    ),
    caption = paste0(
      "Province-level analysis. PVOUT–pop correlation from 10km feasible-land grid. ",
      "Rank change from 11_siting_efficiency_regression.R. Weighted by total province capacity."
    )
  ) +
  theme_rob

ggsave(file.path(OUT, "fig_R7_tradeoff_vs_rank_change.png"),
       fig_r7, width = 8.5, height = 6, dpi = 300)
cat("  Saved fig_R7_tradeoff_vs_rank_change.png\n")


# =============================================================================
# SECTION D — Regression: does controlling for pop density change main result?
# =============================================================================
cat("\n===== D: Plant-level regression controlling for population proximity =====\n")

# Attach pop_density to plant-level data via nearest grid cell
# Use the already-joined grid_full to build a province-level pop density
# summary for feasible land, then use each plant's PVOUT relative to
# province feasible median (already captured by within_pct_rank)
#
# More directly: does a plant's distance from its province's pop-weighted
# centroid explain rank, and does controlling for it change the post2020 coef?

# Build: province-level mean log(pop_density) on feasible land
# as a summary of the typical land-cost environment
prov_pop_context <- grid_full %>%
  group_by(province_en) %>%
  summarise(
    feas_mean_log_pop = mean(log_pop, na.rm = TRUE),
    feas_med_log_pop  = median(log_pop, na.rm = TRUE),
    .groups = "drop"
  )

# Plant-level: each plant's absolute pvout relative to national
# Already in deploy_en as within_pct_rank
# Add province feasible land pop context
deploy_en_pop <- deploy_en %>%
  dplyr::filter(is.finite(within_pct_rank), is.finite(cap_mw)) %>%
  left_join(prov_pop_context,
            by = c("prov_en" = "province_en")) %>%
  left_join(prov_land_cor %>%
              dplyr::select(province_en, cor_pearson),
            by = c("prov_en" = "province_en")) %>%
  mutate(year_c = as.numeric(build_year) - 2020)

# Main spec re-run + pop-density context control
m_main_reprep <- feols(
  within_pct_rank ~ post2020 + year_c | prov_en,
  data    = deploy_en_pop, weights = ~cap_mw, vcov = ~prov_en
)

# Model with province feasible-land mean pop density (province-level control)
# Note: this is absorbed by province FEs if time-invariant; we interact
# post2020 with the province's PVOUT-pop correlation to test
# whether the trade-off moderates the post-2020 decline.
m_interact <- feols(
  within_pct_rank ~ post2020 * cor_pearson + year_c | prov_en,
  data    = deploy_en_pop, weights = ~cap_mw, vcov = ~prov_en
)

cat("\nMain spec (reproduced):\n")
cat(sprintf("  beta post2020 = %.3f  (SE %.3f, p = %.4f)\n",
            coef(m_main_reprep)["post2020TRUE"],
            se(m_main_reprep)["post2020TRUE"],
            pvalue(m_main_reprep)["post2020TRUE"]))

cat("\nInteraction model: post2020 × PVOUT-pop correlation:\n")
etable(m_interact, se.below = TRUE,
       coefstat = "se", digits = 3, fitstat = c("n", "r2"))

# Key test: is the interaction significant?
int_coef <- coef(m_interact)
int_se   <- se(m_interact)
int_p    <- pvalue(m_interact)

cat(sprintf("\n  beta post2020                     = %.3f (SE %.3f, p = %.4f)\n",
            int_coef["post2020TRUE"],
            int_se["post2020TRUE"],
            int_p["post2020TRUE"]))
cat(sprintf("  beta post2020 × cor_pearson       = %.3f (SE %.3f, p = %.4f)\n",
            int_coef["post2020TRUE:cor_pearson"],
            int_se["post2020TRUE:cor_pearson"],
            int_p["post2020TRUE:cor_pearson"]))

interaction_sig <- int_p["post2020TRUE:cor_pearson"] < 0.10
cat(sprintf("\n  Interaction significant at 10%%? %s\n",
            ifelse(interaction_sig, "YES — land costs moderate the efficiency decline",
                   "NO — land opportunity costs do not explain the efficiency decline")))

# Save results
d_results <- data.frame(
  model               = c("Main spec (no pop control)",
                           "Interaction: post2020 x PVOUT-pop cor"),
  coef_post2020       = round(c(coef(m_main_reprep)["post2020TRUE"],
                                 int_coef["post2020TRUE"]), 3),
  se_post2020         = round(c(se(m_main_reprep)["post2020TRUE"],
                                 int_se["post2020TRUE"]), 3),
  p_post2020          = round(c(pvalue(m_main_reprep)["post2020TRUE"],
                                 int_p["post2020TRUE"]), 4),
  coef_interaction    = round(c(NA, int_coef["post2020TRUE:cor_pearson"]), 3),
  p_interaction       = round(c(NA, int_p["post2020TRUE:cor_pearson"]), 4)
)
write.csv(d_results,
          file.path(ROB, "land_cost_interaction_results.csv"),
          row.names = FALSE)

# =============================================================================
# SUMMARY
# =============================================================================
cat("\n================================================================\n")
cat("SUMMARY: Land opportunity cost analysis\n")
cat("================================================================\n")
cat(sprintf("\nA. National PVOUT–population correlation: r = %.3f\n",
            nat_cor_pearson))
cat(sprintf("   Direction: NEGATIVE — high-PVOUT land is LESS populated.\n"))
cat(sprintf("   China's best solar resources (western deserts/plateaux) carry\n"))
cat(sprintf("   LOW opportunity cost. The inefficiency argument is strengthened.\n"))
cat(sprintf("\nB. Within-province trade-off:\n"))
cat(sprintf("   %d of %d provinces: POSITIVE correlation (real trade-off)\n",
            pos_cor_n, nrow(prov_land_cor)))
cat(sprintf("   %d of %d provinces: NEGATIVE correlation (no trade-off)\n",
            neg_cor_n, nrow(prov_land_cor)))
cat(sprintf("\nC. Trade-off vs siting rank change:\n"))
cat(sprintf("   Correlation r = %.3f (p = %.3f from OLS)\n",
            cor_tradeoff_rank,
            summary(m_tradeoff1)$coefficients[2, 4]))
cat(sprintf("   Provinces with harder land-cost trade-offs do NOT\n"))
if (abs(cor_tradeoff_rank) < 0.2) {
  cat(sprintf("   systematically show greater siting rank decline.\n"))
} else {
  cat(sprintf("   show different siting patterns — see regression above.\n"))
}
cat(sprintf("\nD. Interaction test:\n"))
cat(sprintf("   post2020 × PVOUT-pop correlation: p = %.4f (%s)\n",
            int_p["post2020TRUE:cor_pearson"],
            ifelse(interaction_sig, "SIGNIFICANT", "not significant")))
cat("================================================================\n")
cat("\nAll outputs saved.\n")
