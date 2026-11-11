source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# robustness_checks.R
# Run AFTER source("00_setup.R")
# Produces all numbers needed for the four reviewer items:
#   1. Causal wording — regression SE + CI for §7.1 OLS
#   2. Sample construction table numbers
#   3. Threshold sensitivity + leave-one-out
#   4. Uncertainty: Wilson CI for 70.3%, bootstrap CI for Xinjiang share
# =============================================================================

library(dplyr)
library(tidyr)

# =============================================================================
# SECTION A — Sample construction table
# =============================================================================
cat("\n===== A: SAMPLE CONSTRUCTION =====\n")

# Total plants in deploy_df (post all filters)
cat("Final deploy_df plants:", nrow(deploy_df), "\n")
cat("Final deploy_df capacity (GW):", round(sum(deploy_df$cap_mw)/1e3, 1), "\n")
cat("Provinces represented:", length(unique(deploy_df$prov)), "\n")
cat("Quarters covered:", paste(range(deploy_df$year_quarter), collapse = " to "), "\n")

# Cumulative N and capacity at key quarters
key_quarters <- c("2017 Q4", "2018 Q4", "2019 Q4", "2020 Q2",
                  "2021 Q4", "2022 Q4", "2023 Q1", "2023 Q4", "2024 Q2")

panel_by_quarter <- deploy_df %>%
  count(year_quarter, wt = cap_mw, name = "cap_added_mw") %>%
  arrange(year_quarter) %>%
  mutate(
    n_added    = deploy_df %>% count(year_quarter) %>% arrange(year_quarter) %>% pull(n),
    cum_n      = cumsum(n_added),
    cum_cap_gw = cumsum(cap_added_mw) / 1e3
  )

cat("\nCumulative plants and capacity at key quarters:\n")
print(panel_by_quarter %>% filter(year_quarter %in% key_quarters))

# Province-level final N
prov_n <- deploy_df %>%
  group_by(prov) %>%
  summarise(n_plants = n(), cap_gw = round(sum(cap_mw)/1e3, 2)) %>%
  arrange(desc(cap_gw))

cat("\nTop 10 provinces by capacity:\n")
print(head(prov_n, 10))

# Unmatched plant check (already handled in setup — report from upv_sf)
if (exists("upv_sf")) {
  unmatched_n <- sum(is.na(upv_sf$prov))
  cat("\nPlants requiring nearest-province assignment:", unmatched_n, "\n")
}


# =============================================================================
# SECTION B — Threshold sensitivity (top 10%, top 25%, top 20%)
# =============================================================================
cat("\n===== B: THRESHOLD SENSITIVITY =====\n")

thresholds <- c(
  "Q75 (top 25%)" = pvout_q75_national,
  "Q80 (top 20%)" = pvout_q80_national,
  "Q90 (top 10%)" = pvout_q90_national
)

# Capacity-weighted % below each threshold
thresh_results <- lapply(names(thresholds), function(nm) {
  thr <- thresholds[[nm]]

  # unweighted share below threshold
  pct_below_uw <- mean(deploy_df$pvout < thr) * 100

  # capacity-weighted share below threshold
  pct_below_cw <- sum(deploy_df$cap_mw[deploy_df$pvout < thr]) /
                  sum(deploy_df$cap_mw) * 100

  # generation shortfall vs benchmark (plants below threshold only)
  # shortfall per plant (MW) = cap_mw * (thr - pvout) [in same pvout units]
  shortfall_mw_pvout <- sum(
    pmax(0, (thr - deploy_df$pvout) * deploy_df$cap_mw), na.rm = TRUE
  )
  # convert: cap_mw * delta_pvout (kWh/kWp/day * 365) -> MWh/yr -> TWh/yr
  # pvout units are kWh/kWp/day; cap_mw is MW = 1000 kWp
  # shortfall_twh = sum(cap_mw * 1000 * delta_pvout_daily * 365) / 1e9
  shortfall_twh <- sum(
    pmax(0, (thr - deploy_df$pvout) * deploy_df$cap_mw * 1000 * 365),
    na.rm = TRUE
  ) / 1e9

  # capex equivalent at 0.20 RMB/kWh LCOE, 3450 RMB/kW unit investment
  # foregone TWh * 1e9 kWh * 0.20 RMB = foregone revenue (RMB)
  # divide by (3450 RMB/kW) to get capex-equivalent kW, /1e9 for bn RMB
  capex_bn_rmb <- shortfall_twh * 1e9 * 0.20 / (3450 * 1e3) / 1e6
  # simplified: shortfall_twh * 1e9 * 0.20 / 3450000 / 1e6
  # = shortfall_twh * 0.20 * 1e9 / 3.45e9 = shortfall_twh * 57.97
  capex_bn_rmb2 <- shortfall_twh * (0.20 * 1e9) / (3450 * 1e6)

  data.frame(
    threshold       = nm,
    pvout_cutoff    = round(thr, 3),
    pct_below_unwtd = round(pct_below_uw, 1),
    pct_below_capwt = round(pct_below_cw, 1),
    shortfall_twh   = round(shortfall_twh, 1),
    capex_bn_rmb    = round(capex_bn_rmb2, 1)
  )
}) %>% bind_rows()

cat("\nThreshold sensitivity table:\n")
print(thresh_results)

# Also report % of capacity on TOP of each threshold (post-2020 only)
cat("\nPost-2020 capacity on top-resource land by threshold:\n")
post2020 <- deploy_df %>%
  filter(build_year >= 2020 | (build_year == 2020 & build_quarter >= 2))

for (nm in names(thresholds)) {
  thr <- thresholds[[nm]]
  pct <- sum(post2020$cap_mw[post2020$pvout >= thr]) / sum(post2020$cap_mw) * 100
  cat(sprintf("  %s: %.1f%% of post-Q2 2020 capacity above threshold\n", nm, pct))
}


# =============================================================================
# SECTION C — Leave-one-province-out
# =============================================================================
cat("\n===== C: LEAVE-ONE-PROVINCE-OUT =====\n")

provinces <- unique(deploy_df$prov)

# For each excluded province, compute quarterly capacity-weighted mean PVOUT
loo_quarterly <- lapply(provinces, function(p) {
  deploy_df %>%
    filter(prov != p) %>%
    group_by(year_quarter) %>%
    summarise(cw_pvout = weighted.mean(pvout, cap_mw), .groups = "drop") %>%
    mutate(excluded = p)
}) %>% bind_rows()

# Full sample quarterly series
full_quarterly <- deploy_df %>%
  group_by(year_quarter) %>%
  summarise(cw_pvout = weighted.mean(pvout, cap_mw), .groups = "drop") %>%
  mutate(excluded = "NONE (full sample)")

# Pre/post Q2 2020 mean for each specification
loo_summary <- bind_rows(full_quarterly, loo_quarterly) %>%
  mutate(period = case_when(
    year_quarter <  "2020 Q2" ~ "pre",
    year_quarter >= "2020 Q2" & year_quarter <= "2023 Q1" ~ "decline",
    year_quarter >  "2023 Q1" ~ "recovery",
    TRUE ~ NA_character_
  )) %>%
  group_by(excluded, period) %>%
  summarise(mean_pvout = mean(cw_pvout), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = mean_pvout) %>%
  mutate(
    decline_change  = round(decline - pre, 4),
    recovery_change = round(recovery - decline, 4)
  ) %>%
  arrange(decline_change)

cat("\nLeave-one-out: effect on Q2 2020 – Q1 2023 decline (sorted by largest decline):\n")
print(loo_summary %>% select(excluded, pre, decline, recovery, decline_change, recovery_change))

# Focus: Xinjiang, Qinghai, both
focus <- c("NONE (full sample)", "Xinjiang", "Xinjiang Uygur", "Qinghai",
           "Inner Mongolia", "Neimenggu", "Hebei")

# Try to match actual province names (may vary by GADM encoding)
xinjiang_name <- provinces[grepl("Xinjiang|Xin|新疆", provinces, ignore.case = TRUE)][1]
qinghai_name  <- provinces[grepl("Qinghai|青海", provinces, ignore.case = TRUE)][1]
cat("\nXinjiang province label in data:", xinjiang_name, "\n")
cat("Qinghai province label in data:", qinghai_name, "\n")

cat("\nKey LOO specifications:\n")
key_specs <- c("NONE (full sample)", xinjiang_name, qinghai_name)
print(loo_summary %>% filter(excluded %in% key_specs))

# Double exclusion: both Xinjiang and Qinghai
both_excl <- deploy_df %>%
  filter(!prov %in% c(xinjiang_name, qinghai_name)) %>%
  group_by(year_quarter) %>%
  summarise(cw_pvout = weighted.mean(pvout, cap_mw), .groups = "drop") %>%
  mutate(excluded = "Both Xinjiang + Qinghai") %>%
  mutate(period = case_when(
    year_quarter <  "2020 Q2" ~ "pre",
    year_quarter >= "2020 Q2" & year_quarter <= "2023 Q1" ~ "decline",
    year_quarter >  "2023 Q1" ~ "recovery",
    TRUE ~ NA_character_
  )) %>%
  group_by(excluded, period) %>%
  summarise(mean_pvout = mean(cw_pvout), .groups = "drop") %>%
  pivot_wider(names_from = period, values_from = mean_pvout) %>%
  mutate(
    decline_change  = round(decline - pre, 4),
    recovery_change = round(recovery - decline, 4)
  )

cat("\nDouble exclusion (Xinjiang + Qinghai):\n")
print(both_excl)


# =============================================================================
# SECTION D — Uncertainty reporting
# =============================================================================
cat("\n===== D: UNCERTAINTY REPORTING =====\n")

# --- D1. Wilson CI for 70.3% (share below national median) ---
n_total   <- nrow(deploy_df)
n_below   <- sum(deploy_df$pvout < pvout_national_median, na.rm = TRUE)
pct_below <- n_below / n_total * 100

# Wilson CI (better than normal approximation for proportions)
wilson_ci <- function(x, n, conf = 0.95) {
  z   <- qnorm(1 - (1 - conf) / 2)
  p   <- x / n
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half   <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(lower = centre - half, upper = centre + half)
}

ci_below <- wilson_ci(n_below, n_total)
cat(sprintf("\nShare of plants below national PVOUT median:\n"))
cat(sprintf("  N below / N total: %d / %d\n", n_below, n_total))
cat(sprintf("  Point estimate: %.1f%%\n", pct_below))
cat(sprintf("  95%% Wilson CI: [%.1f%%, %.1f%%]\n",
            ci_below["lower"] * 100, ci_below["upper"] * 100))

# Capacity-weighted version
cap_below_pct <- sum(deploy_df$cap_mw[deploy_df$pvout < pvout_national_median]) /
                 sum(deploy_df$cap_mw) * 100
cat(sprintf("  Capacity-weighted share below median: %.1f%%\n", cap_below_pct))


# --- D2. Bootstrap CI for Xinjiang's share of national deterioration ---
cat("\nBootstrap CI for Xinjiang share of post-Q2 2020 deterioration:\n")
cat("(1000 resamples — may take ~30 seconds)\n")

post_df <- deploy_df %>% filter(
  build_year > 2020 | (build_year == 2020 & build_quarter >= 2)
)

set.seed(42)
boot_shares <- replicate(1000, {
  # Resample within each province
  boot_sample <- post_df %>%
    group_by(prov) %>%
    slice_sample(prop = 1, replace = TRUE) %>%
    ungroup() %>%
    mutate(gap = pmax(pvout_q75_national - pvout, 0),
           contrib = cap_mw * gap)

  total_contrib <- sum(boot_sample$contrib, na.rm = TRUE)
  xinj_contrib  <- sum(boot_sample$contrib[boot_sample$prov == xinjiang_name],
                       na.rm = TRUE)
  if (total_contrib == 0) return(NA_real_)
  xinj_contrib / total_contrib * 100
})

boot_shares <- boot_shares[is.finite(boot_shares)]
cat(sprintf("  Point estimate (observed): %.1f%%\n",
            sum(post_df$cap_mw[post_df$prov == xinjiang_name] *
                pmax(pvout_q75_national - post_df$pvout[post_df$prov == xinjiang_name], 0)) /
            sum(post_df$cap_mw * pmax(pvout_q75_national - post_df$pvout, 0)) * 100))
cat(sprintf("  Bootstrap mean:            %.1f%%\n", mean(boot_shares)))
cat(sprintf("  95%% bootstrap CI:         [%.1f%%, %.1f%%]\n",
            quantile(boot_shares, 0.025), quantile(boot_shares, 0.975)))


# --- D3. Province-level growth vs PVOUT rank regression (Fig 7.12) ---
cat("\nFig 7.12 regression: growth rate vs within-province PVOUT rank change:\n")

# Pre/post capacity-weighted mean within-province PVOUT rank
rank_by_prov_period <- deploy_df %>%
  group_by(prov) %>%
  mutate(prov_rank = rank(pvout) / n() * 100) %>%    # within-province percentile
  ungroup() %>%
  mutate(period = ifelse(
    build_year < 2020 | (build_year == 2020 & build_quarter < 2),
    "pre", "post"
  )) %>%
  group_by(prov, period) %>%
  summarise(cw_rank = weighted.mean(prov_rank, cap_mw),
            cap_gw  = sum(cap_mw) / 1e3,
            .groups = "drop")

rank_wide <- rank_by_prov_period %>%
  pivot_wider(names_from = period, values_from = c(cw_rank, cap_gw)) %>%
  filter(!is.na(cw_rank_pre), !is.na(cw_rank_post)) %>%
  mutate(
    rank_change    = cw_rank_post - cw_rank_pre,
    log_growth     = log(cap_gw_post / cap_gw_pre)
  ) %>%
  filter(is.finite(log_growth))

# Bivariate regression
m_bivar <- lm(rank_change ~ log_growth, data = rank_wide)
cat("\n  Bivariate (rank_change ~ log_growth):\n")
print(summary(m_bivar)$coefficients)
cat(sprintf("  R²: %.3f\n", summary(m_bivar)$r.squared))
cat("  95%% CI for log_growth coefficient:\n")
print(confint(m_bivar))

# Controlling for pre-2020 rank
m_controlled <- lm(rank_change ~ log_growth + cw_rank_pre, data = rank_wide)
cat("\n  Controlled (+ pre-2020 rank):\n")
print(summary(m_controlled)$coefficients)
cat(sprintf("  R²: %.3f\n", summary(m_controlled)$r.squared))
cat("  95%% CI for log_growth coefficient:\n")
print(confint(m_controlled)["log_growth", , drop = FALSE])


# =============================================================================
# SECTION E — Summary output file
# =============================================================================
cat("\n===== E: SAVING RESULTS =====\n")

results_dir <- file.path(
  FTS_ROOT,
  "robustness_outputs"
)
dir.create(results_dir, showWarnings = FALSE)

write.csv(thresh_results,  file.path(results_dir, "threshold_sensitivity.csv"),  row.names = FALSE)
write.csv(loo_summary,     file.path(results_dir, "leave_one_out.csv"),           row.names = FALSE)
write.csv(panel_by_quarter,file.path(results_dir, "sample_construction_panel.csv"), row.names = FALSE)
write.csv(prov_n,          file.path(results_dir, "province_n_cap.csv"),          row.names = FALSE)
write.csv(rank_wide,       file.path(results_dir, "fig712_regression_data.csv"),  row.names = FALSE)

cat("Results saved to:", results_dir, "\n")
cat("\nAll done. Paste the printed numbers into the thesis draft.\n")
