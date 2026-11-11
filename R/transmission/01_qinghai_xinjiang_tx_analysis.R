source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# 01_qinghai_xinjiang_tx_analysis.R
# Qinghai vs Xinjiang post-shuangtan transmission cost analysis.
#
# Grid nodes: Figshare substation dataset (2020), >=220kV
#   Chen et al. (2026), doi: 10.6084/m9.figshare.31569094
# TX cost formula: Wang et al. (2022), gf = 0.0002 × dist_km
# Electricity price deflator: 0.37 CNY/kWh
# Distance cap: 600 km
#
# Outputs (to this folder):
#   fig_01a_pvout_dist_scatter.pdf/.png
#   fig_01b_ngv_distribution.pdf/.png
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(ggplot2)
  library(scales)
})

P3  <- FTS_PROCESSED
P4  <- FTS_ROOT
OUT <- file.path(P4, "transmission")

GF_PER_KM  <- 0.0002
ELEC_PRICE <- 0.37
DIST_CAP   <- 600
CURT_PROVS <- c("Xinjiang", "Gansu", "Ningxia")
CURT_RATE  <- 0.20

# =============================================================================
# 1. Load plant data
# =============================================================================
cat("Loading plant data...\n")

# Attributes: pvout, cap_mw, prov_name, build_year, build_quarter
deploy_df <- readRDS(file.path(P4, "cache_1km/deploy_df_1km.rds"))
cat(sprintf("  Plants in deploy_df: %d\n", nrow(deploy_df)))

# Geometries (for distance computation)
upv_sf <- readRDS(file.path(P3, "UPV_core_GRW_2024Q2_ge1MW.rds"))
cat(sprintf("  Plants in upv_sf:    %d\n", nrow(upv_sf)))

# Align: deploy_df may have one fewer row — match on build_year + build_quarter + area
# Safer: use row position after ensuring same ordering
# deploy_df was derived from upv_sf in 00_setup.R; drop any extra rows from upv_sf
n <- nrow(deploy_df)
upv_geom <- upv_sf[seq_len(n), ]   # trim to same length

# Get CRS from PVOUT raster
pvout_r  <- rast(file.path(P3, "pvout_annual_10km.tif"))
grid_crs <- crs(pvout_r, proj = TRUE)

# =============================================================================
# 2. Post-shuangtan filter (Q3 2020 onward)
# =============================================================================
post_idx <- which(deploy_df$build_year > 2020 |
                  (deploy_df$build_year == 2020 & deploy_df$build_quarter >= 3))

post_df   <- deploy_df[post_idx, ]
post_geom <- upv_geom[post_idx, ]

cat(sprintf("Post-shuangtan plants: %d\n", nrow(post_df)))
cat(sprintf("  Qinghai:  %d\n", sum(post_df$prov_name == "Qinghai")))
cat(sprintf("  Xinjiang: %d\n", sum(post_df$prov_name == "Xinjiang")))

# =============================================================================
# 3. Distance to nearest 2020 >=220kV substation
# =============================================================================
cat("Computing distances to nearest 2020 grid node...\n")

nodes    <- read.csv(file.path(OUT, "grid_nodes_2020.csv"))
nodes_sf <- st_as_sf(nodes, coords = c("lng", "lat"), crs = "EPSG:4326") %>%
            st_transform(grid_crs)

plants_sf <- post_geom %>%
             st_centroid() %>%
             st_transform(grid_crs)

nearest_idx <- st_nearest_feature(plants_sf, nodes_sf)
dist_km <- pmin(
  as.numeric(st_distance(plants_sf, nodes_sf[nearest_idx, ], by_element = TRUE)) / 1000,
  DIST_CAP
)

post_df <- post_df %>%
  mutate(
    dist_km    = dist_km,
    tx_penalty = (GF_PER_KM * dist_km / ELEC_PRICE) * pvout,
    curt_pen   = ifelse(prov_name %in% CURT_PROVS, CURT_RATE * pvout, 0),
    ngv        = pvout - tx_penalty - curt_pen
  )

cat(sprintf("Distance summary — median: %.1f km  mean: %.1f km  max: %.1f km\n",
            median(dist_km), mean(dist_km), max(dist_km)))

qh <- post_df %>% filter(prov_name == "Qinghai")
xj <- post_df %>% filter(prov_name == "Xinjiang")

# =============================================================================
# 4. Summary statistics
# =============================================================================
cat("\n=== POST-SHUANGTAN SUMMARY (2020 grid nodes, >=220kV) ===\n")

summarise_prov <- function(df, label) {
  cat(sprintf("\n%s (%d plants, %.1f GW):\n",
              label, nrow(df), sum(df$cap_mw) / 1000))
  cat(sprintf("  Cap-weighted PVOUT:      %.0f kWh/kWp/yr\n",
              weighted.mean(df$pvout, df$cap_mw)))
  cat(sprintf("  Cap-weighted dist:       %.1f km\n",
              weighted.mean(df$dist_km, df$cap_mw)))
  cat(sprintf("  Cap-weighted TX penalty: %.1f kWh/kWp/yr\n",
              weighted.mean(df$tx_penalty, df$cap_mw)))
  cat(sprintf("  Cap-weighted NGV:        %.0f kWh/kWp/yr\n",
              weighted.mean(df$ngv, df$cap_mw)))
  cat(sprintf("  Median PVOUT:            %.0f kWh/kWp/yr\n", median(df$pvout)))
  cat(sprintf("  Median dist:             %.1f km\n",          median(df$dist_km)))
}

summarise_prov(qh, "Qinghai")
summarise_prov(xj, "Xinjiang")

pvout_gap_cw <- weighted.mean(qh$pvout, qh$cap_mw) - weighted.mean(xj$pvout, xj$cap_mw)
pvout_gap_uw <- mean(qh$pvout) - mean(xj$pvout)
cat(sprintf("\nPVOUT gap — cap-weighted: %.1f  unweighted: %.1f kWh/kWp/yr\n",
            pvout_gap_cw, pvout_gap_uw))

# =============================================================================
# 5. OLS: PVOUT ~ dist_km
# =============================================================================
cat("\n=== OLS: PVOUT ~ distance ===\n")
lm_qh <- lm(pvout ~ dist_km, data = qh)
lm_xj <- lm(pvout ~ dist_km, data = xj)
lm_qh_cw <- lm(pvout ~ dist_km, data = qh, weights = cap_mw)
lm_xj_cw <- lm(pvout ~ dist_km, data = xj, weights = cap_mw)

cat(sprintf("Qinghai  (unweighted): slope=%+.3f/km  intercept=%.0f  R²=%.3f\n",
            coef(lm_qh)[2], coef(lm_qh)[1], summary(lm_qh)$r.squared))
cat(sprintf("Xinjiang (unweighted): slope=%+.3f/km  intercept=%.0f  R²=%.3f\n",
            coef(lm_xj)[2], coef(lm_xj)[1], summary(lm_xj)$r.squared))
cat(sprintf("Qinghai  (cap-wtd):    slope=%+.3f/km  intercept=%.0f  R²=%.3f\n",
            coef(lm_qh_cw)[2], coef(lm_qh_cw)[1], summary(lm_qh_cw)$r.squared))
cat(sprintf("Xinjiang (cap-wtd):    slope=%+.3f/km  intercept=%.0f  R²=%.3f\n",
            coef(lm_xj_cw)[2], coef(lm_xj_cw)[1], summary(lm_xj_cw)$r.squared))

# =============================================================================
# 6. Oaxaca-Blinder decomposition
# =============================================================================
ob_decomp <- function(df_a, df_b, weighted = FALSE) {
  if (weighted) {
    lm_a <- lm(pvout ~ dist_km, data = df_a, weights = cap_mw)
    lm_b <- lm(pvout ~ dist_km, data = df_b, weights = cap_mw)
    mean_pvout_a <- weighted.mean(df_a$pvout, df_a$cap_mw)
    mean_pvout_b <- weighted.mean(df_b$pvout, df_b$cap_mw)
    mean_dist_a  <- weighted.mean(df_a$dist_km, df_a$cap_mw)
    mean_dist_b  <- weighted.mean(df_b$dist_km, df_b$cap_mw)
  } else {
    lm_a <- lm(pvout ~ dist_km, data = df_a)
    lm_b <- lm(pvout ~ dist_km, data = df_b)
    mean_pvout_a <- mean(df_a$pvout)
    mean_pvout_b <- mean(df_b$pvout)
    mean_dist_a  <- mean(df_a$dist_km)
    mean_dist_b  <- mean(df_b$dist_km)
  }

  gap          <- mean_pvout_a - mean_pvout_b
  constant_eff <- coef(lm_a)[1] - coef(lm_b)[1]
  endowment    <- coef(lm_b)[2] * (mean_dist_a - mean_dist_b)
  coef_eff     <- (coef(lm_a)[2] - coef(lm_b)[2]) * mean_dist_a

  cat(sprintf("\nOaxaca-Blinder (%s):\n", ifelse(weighted, "cap-weighted", "unweighted")))
  cat(sprintf("  Total gap:        %+.1f kWh/kWp/yr\n", gap))
  cat(sprintf("  Constant effect:  %+.1f  (%+.0f%%)\n", constant_eff, constant_eff/gap*100))
  cat(sprintf("  Endowment effect: %+.1f  (%+.0f%%)\n", endowment,    endowment/gap*100))
  cat(sprintf("  Coeff effect:     %+.1f  (%+.0f%%)\n", coef_eff,     coef_eff/gap*100))

  invisible(list(gap=gap, constant=constant_eff, endowment=endowment, coef=coef_eff))
}

cat("\n--- Qinghai vs Xinjiang ---")
ob_cw <- ob_decomp(qh, xj, weighted = TRUE)
ob_uw <- ob_decomp(qh, xj, weighted = FALSE)

endow_diff_pp <- abs(ob_cw$endowment / ob_cw$gap - ob_uw$endowment / ob_uw$gap) * 100
cat(sprintf("\nEndowment share difference across weighting: %.1f pp (robustness check)\n",
            endow_diff_pp))

# =============================================================================
# 7. Counterfactual
# =============================================================================
cat("\n=== COUNTERFACTUAL ===\n")
cf_pvout <- predict(lm_qh_cw, newdata = data.frame(dist_km = xj$dist_km))
cf_gain  <- weighted.mean(cf_pvout - xj$pvout, xj$cap_mw)
cf_twh   <- cf_gain * sum(xj$cap_mw) / 1e6

cat(sprintf("If XJ resource quality matched QH slope at same distances:\n"))
cat(sprintf("  Mean PVOUT gain:    +%.0f kWh/kWp/yr\n", cf_gain))
cat(sprintf("  Annual output gain: ~%.1f TWh/yr (%.1f GW XJ capacity)\n",
            cf_twh, sum(xj$cap_mw)/1000))

# =============================================================================
# 8. Figure A: PVOUT vs distance scatter
# =============================================================================
cat("\nBuilding figures...\n")

plot_df <- post_df %>%
  filter(prov_name %in% c("Qinghai", "Xinjiang")) %>%
  rename(Province = prov_name)

prov_cols <- c("Qinghai" = "#2e9e56", "Xinjiang" = "#e8936a")

fig_scatter <- ggplot(plot_df, aes(x = dist_km, y = pvout, colour = Province)) +
  geom_point(aes(size = cap_mw), alpha = 0.55, shape = 16) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9) +
  scale_colour_manual(values = prov_cols) +
  scale_size_continuous(name = "Capacity (MW)", range = c(0.8, 6),
                        breaks = c(100, 500, 1000, 3000)) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  labs(
    title    = "PVOUT vs. Distance to Nearest Substation (≥220 kV)",
    subtitle = "Post-shuangtan plants (Q3 2020 – Q2 2024), Qinghai and Xinjiang",
    x        = "Distance to nearest ≥220 kV substation (km)",
    y        = "PVOUT (kWh/kWp/yr)",
    caption  = paste0(
      "Grid nodes: Chen et al. (2026), ≥220 kV substations, 2020 vintage (n = ",
      nrow(nodes_sf), ").  ",
      "PVOUT: Global Solar Atlas v2.0 (ESMAP, 2019).\n",
      "TX cost: Wang et al. (2022), gf = 0.0002 × d (km); penalty = (gf / 0.37) × PVOUT.  ",
      "Plant polygons: GRW Q2 2024 (Robinson et al., 2025)."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    plot.subtitle   = element_text(size = 9, colour = "grey30"),
    plot.caption    = element_text(size = 7.5, colour = "grey45", hjust = 0,
                                   margin = margin(t = 8)),
    legend.position = "right",
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin     = margin(8, 8, 8, 8)
  )

ggsave(file.path(OUT, "fig_01a_pvout_dist_scatter.pdf"),
       fig_scatter, width = 10, height = 6.5)
ggsave(file.path(OUT, "fig_01a_pvout_dist_scatter.png"),
       fig_scatter, width = 10, height = 6.5, dpi = 250)

# =============================================================================
# 9. Figure B: NGV distribution
# =============================================================================
ngv_means <- plot_df %>%
  group_by(Province) %>%
  summarise(med = weighted.mean(ngv, cap_mw), .groups = "drop")

fig_ngv <- ggplot(plot_df, aes(x = ngv, fill = Province, weight = cap_mw)) +
  geom_histogram(binwidth = 40, alpha = 0.70, position = "identity", colour = "white") +
  geom_vline(data = ngv_means, aes(xintercept = med, colour = Province),
             linewidth = 1, linetype = "dashed") +
  scale_fill_manual(values  = prov_cols) +
  scale_colour_manual(values = prov_cols) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = unit_format(unit = "GW", scale = 1e-3)) +
  labs(
    title    = "Net Generation Value Distribution",
    subtitle = "Capacity-weighted; NGV = PVOUT − TX penalty − curtailment penalty (Xinjiang: −20% PVOUT)",
    x        = "NGV (kWh/kWp/yr)",
    y        = "Installed capacity",
    caption  = "Dashed lines: capacity-weighted mean NGV per province."
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    plot.subtitle   = element_text(size = 9, colour = "grey30"),
    plot.caption    = element_text(size = 7.5, colour = "grey45", hjust = 0),
    legend.position = "right",
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin     = margin(8, 8, 8, 8)
  )

ggsave(file.path(OUT, "fig_01b_ngv_distribution.pdf"),
       fig_ngv, width = 10, height = 6)
ggsave(file.path(OUT, "fig_01b_ngv_distribution.png"),
       fig_ngv, width = 10, height = 6, dpi = 250)

cat("Saved: fig_01a_pvout_dist_scatter, fig_01b_ngv_distribution\n")
cat("Done.\n")
