source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# 11b_siting_efficiency_1km.R  —  Thesis Project 4
# Replication of 11_siting_efficiency_regression.R at 1 km resolution
#
# Key change vs 11_siting_efficiency_regression.R:
#   - Uses pvout_annual.tif  (~900 m / 0.00833° pixel, 32.1 M cells) instead
#     of pvout_annual_10km.tif (10 km, 97 k cells)
#   - Elevation/slope resampled from 10 km → 1 km for feasibility mask
#   - Province raster built by rasterizing GADM boundaries to 1 km grid
#   - Within-province PVOUT percentile ranks computed on the 1 km distribution
#     (finer bins → more precise ecdf per province)
#   - Plant PVOUT re-extracted from the 1 km source raster
#
# All expensive steps are cached; subsequent runs load from RDS/tif.
#
# Outputs (figures/):
#   fig_11b_1_annual_rank_timeseries.png
#   fig_11b_2_province_pre_post.png
#   fig_11b_3_growth_vs_delta_rank.png
#   fig_11b_4_regression_coefs.png
#   fig_11b_compare_10k_vs_1km.png   — direct comparison of both resolutions
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(fixest)
  library(patchwork)
  library(broom)
})

P4        <- FTS_ROOT
P3proc    <- FTS_PROCESSED
OUT       <- file.path(P4, "figures")
CACHE     <- file.path(P4, "cache_1km")
dir.create(OUT,   showWarnings = FALSE)
dir.create(CACHE, showWarnings = FALSE)

# Paths to source rasters
PVOUT_1KM_PATH  <- file.path(P3proc, "pvout_annual.tif")           # ~1 km WGS84
ELEV_10KM_PATH  <- file.path(P3proc, "elevation_10km.tif")         # 10 km CEA
SLOPE_10KM_PATH <- file.path(P3proc, "slope_deg_10km.tif")         # 10 km CEA

# Cache paths
PROV_RASTER_PATH  <- file.path(CACHE, "province_id_1km.tif")
FEAS_RASTER_PATH  <- file.path(CACHE, "feasibility_mask_1km.tif")
DEPLOY_1KM_PATH   <- file.path(CACHE, "deploy_df_1km.rds")

cat("================================================\n")
cat("11b: Within-province siting efficiency at 1 km\n")
cat("================================================\n\n")

# =============================================================================
# 0. Check 00_setup.R objects exist
# =============================================================================
required <- c("deploy_df", "upv_sf", "china_provinces",
              "pvout_national_median", "pvout_q75_national")
missing_obj <- required[!sapply(required, exists)]
if (length(missing_obj) > 0) {
  stop("Missing objects from 00_setup.R: ", paste(missing_obj, collapse = ", "),
       "\nPlease run source('00_setup.R') first.")
}

# Shuangtan boundary (consistent with rest of thesis)
# Pre  = build_year < 2020  OR  (build_year == 2020 AND build_quarter == 1)
# Post = (build_year == 2020 AND build_quarter >= 2)  OR  build_year > 2020
deploy_df <- deploy_df %>%
  mutate(post2020 = as.integer(
    build_year > 2020 | (build_year == 2020 & build_quarter >= 2)
  ))

# =============================================================================
# 1. Load 1 km PVOUT raster
# =============================================================================
cat("Loading 1 km PVOUT raster...\n")
pvout_1km <- rast(PVOUT_1KM_PATH)
cat(sprintf("  Dimensions : %d x %d = %s cells\n",
            nrow(pvout_1km), ncol(pvout_1km),
            format(ncell(pvout_1km), big.mark = ",")))
cat(sprintf("  Resolution : %.6f deg (~%.0f m at 35N)\n",
            res(pvout_1km)[1], res(pvout_1km)[1] * 111000))
cat(sprintf("  CRS        : EPSG:%s\n", crs(pvout_1km, describe = TRUE)$code))

# =============================================================================
# 2. Build 1 km feasibility mask (cache)
# =============================================================================
if (!file.exists(FEAS_RASTER_PATH)) {
  cat("\nBuilding 1 km feasibility mask (resample from 10 km)...\n")

  elev_10km  <- rast(ELEV_10KM_PATH)
  slope_10km <- rast(SLOPE_10KM_PATH)

  # Reproject pvout template to CEA to match elevation/slope rasters
  pvout_template_cea <- project(pvout_1km, crs(elev_10km), res = res(pvout_1km)[1] * 111000,
                                 method = "near")
  cat("  Resampling elevation to 1 km CEA...\n")
  elev_1km  <- resample(elev_10km,  pvout_template_cea, method = "bilinear")
  cat("  Resampling slope to 1 km CEA...\n")
  slope_1km <- resample(slope_10km, pvout_template_cea, method = "bilinear")

  cat("  Applying feasibility thresholds (elev > 5000m OR slope > 50 deg)...\n")
  feas_cea  <- (elev_1km <= 5000) & (slope_1km <= 50)

  # Reproject feasibility mask back to WGS84 to match PVOUT raster
  cat("  Reprojecting feasibility mask to WGS84...\n")
  feas_wgs <- project(feas_cea, pvout_1km, method = "near")
  feas_wgs <- resample(feas_wgs, pvout_1km, method = "near")

  writeRaster(feas_wgs, FEAS_RASTER_PATH, overwrite = TRUE)
  cat(sprintf("  Cached to %s\n", FEAS_RASTER_PATH))
  rm(elev_10km, slope_10km, pvout_template_cea, elev_1km, slope_1km, feas_cea)
  gc()
} else {
  cat("Loading cached 1 km feasibility mask...\n")
}
feas_1km <- rast(FEAS_RASTER_PATH)

# =============================================================================
# 3. Build province ID raster at 1 km (cache)
# =============================================================================
if (!file.exists(PROV_RASTER_PATH)) {
  cat("\nBuilding 1 km province ID raster (rasterizing GADM boundaries)...\n")
  cat("  This takes ~3-5 minutes for 32 M cells. Caching result.\n")

  # Get province name column (same logic as 00_setup.R)
  prov_col <- names(china_provinces)[
    sapply(names(china_provinces), function(n) {
      x <- china_provinces[[n]]
      is.character(x) | is.factor(x)
    })
  ]
  # Prefer the English 'name' column if it exists
  if ("name" %in% prov_col) {
    prov_col <- "name"
  } else {
    prov_col <- prov_col[prov_col != "geometry"][1]
  }
  cat(sprintf("  Using province column: '%s'\n", prov_col))

  # Assign integer province IDs
  prov_df <- data.frame(
    prov_name = as.character(china_provinces[[prov_col]]),
    prov_id   = seq_len(nrow(china_provinces))
  )
  prov_vect <- vect(china_provinces)
  prov_vect$prov_id <- prov_df$prov_id

  cat("  Rasterizing (this is the slow step)...\n")
  prov_raster <- rasterize(prov_vect, pvout_1km, field = "prov_id")
  writeRaster(prov_raster, PROV_RASTER_PATH, overwrite = TRUE)

  # Also save the ID → name lookup
  saveRDS(prov_df, file.path(CACHE, "province_id_lookup.rds"))
  cat(sprintf("  Cached province raster to %s\n", PROV_RASTER_PATH))
  rm(prov_vect, prov_raster)
  gc()
} else {
  cat("Loading cached province ID raster...\n")
}
prov_raster  <- rast(PROV_RASTER_PATH)
prov_id_lut  <- readRDS(file.path(CACHE, "province_id_lookup.rds"))

# =============================================================================
# 4. Re-extract plant-level PVOUT at 1 km + compute percentile ranks (cache)
# =============================================================================
if (!file.exists(DEPLOY_1KM_PATH)) {
  cat("\nExtracting 1 km PVOUT for each UPV plant centroid...\n")

  # Plant centroids in WGS84
  # upv_sf has 16714 rows; deploy_df has 16713 (one NA-PVOUT plant dropped).
  # Identify the dropped row: whichever upv_sf row has NA 10km PVOUT.
  pvout_10km_r <- rast(file.path(P3proc, "pvout_annual_10km.tif"))
  cents_all    <- st_centroid(upv_sf) %>% st_transform(6933)
  pvout_10km_vals <- terra::extract(pvout_10km_r, vect(cents_all))[[2]]
  keep_rows    <- which(!is.na(pvout_10km_vals))  # length == 16713

  plant_pts <- vect(st_transform(cents_all, 4326))

  pvout_ext <- terra::extract(pvout_1km,   plant_pts)
  prov_ext  <- terra::extract(prov_raster, plant_pts)

  # Filter to the same 16713 rows as deploy_df, then bind
  deploy_1km <- deploy_df %>%
    mutate(pvout_1km   = pvout_ext[[2]][keep_rows],
           prov_id_1km = prov_ext[[2]][keep_rows])

  # Merge province name from lookup
  deploy_1km <- deploy_1km %>%
    left_join(prov_id_lut, by = c("prov_id_1km" = "prov_id"))

  # -----------------------------------------------------------------------
  # Build within-province PVOUT ECDF at 1 km resolution
  # -----------------------------------------------------------------------
  cat("Building within-province 1 km PVOUT distributions...\n")

  # Extract all feasible cell values province-by-province using terra
  # Stack: pvout, feasibility, province ID into a 3-layer raster
  cat("  Stacking rasters for cell-level extraction...\n")
  stack_r <- c(pvout_1km, feas_1km, prov_raster)
  names(stack_r) <- c("pvout", "feasible", "prov_id")

  # Extract values for all non-NA cells (memory-efficient chunk approach)
  cat("  Reading raster values (may take 1-2 min for 32 M cells)...\n")
  vals <- as.data.frame(stack_r, na.rm = FALSE)
  cat(sprintf("  Total cells: %s\n", format(nrow(vals), big.mark = ",")))

  # Keep only feasible, valid cells with province assignments
  vals_feas <- vals %>%
    filter(!is.na(pvout), !is.na(prov_id), feasible == 1) %>%
    select(pvout, prov_id)
  cat(sprintf("  Feasible cells with province: %s\n",
              format(nrow(vals_feas), big.mark = ",")))
  rm(vals); gc()

  # Build ECDF per province (stored as function)
  cat("  Computing per-province ECDF functions...\n")
  prov_ecdfs <- vals_feas %>%
    group_by(prov_id) %>%
    summarise(ecdf_fn = list(ecdf(pvout)),
              n_cells  = n(),
              pvout_med = median(pvout),
              .groups = "drop")
  rm(vals_feas); gc()

  # For each plant, apply its province's ECDF → percentile rank
  cat("  Assigning within-province percentile ranks to plants...\n")
  deploy_1km$pct_rank_1km <- NA_real_

  for (pid in unique(deploy_1km$prov_id_1km[!is.na(deploy_1km$prov_id_1km)])) {
    idx      <- which(deploy_1km$prov_id_1km == pid & !is.na(deploy_1km$pvout_1km))
    ecdf_row <- prov_ecdfs[prov_ecdfs$prov_id == pid, ]
    if (nrow(ecdf_row) == 0 || length(idx) == 0) next
    fn <- ecdf_row$ecdf_fn[[1]]
    deploy_1km$pct_rank_1km[idx] <- fn(deploy_1km$pvout_1km[idx]) * 100
  }

  # Also record province-level feasible cell counts
  deploy_1km <- deploy_1km %>%
    left_join(prov_ecdfs %>% select(prov_id, n_cells_1km = n_cells,
                                    prov_pvout_med_1km = pvout_med),
              by = c("prov_id_1km" = "prov_id"))

  saveRDS(deploy_1km, DEPLOY_1KM_PATH)
  cat(sprintf("  Cached plant-level 1 km data to %s\n", DEPLOY_1KM_PATH))
  rm(stack_r, prov_ecdfs); gc()

} else {
  cat("Loading cached plant-level 1 km data...\n")
  deploy_1km <- readRDS(DEPLOY_1KM_PATH)
}

# Quick coverage check
n_total  <- nrow(deploy_1km)
n_ranked <- sum(!is.na(deploy_1km$pct_rank_1km))
cat(sprintf("  Plants with 1km rank: %d / %d  (%.1f%%)\n",
            n_ranked, n_total, 100 * n_ranked / n_total))

# Working subset: valid rank + capacity
df <- deploy_1km %>%
  filter(!is.na(pct_rank_1km), !is.na(cap_mw), cap_mw > 0,
         !is.na(prov_name), !is.na(build_year))

cat(sprintf("  Analysis subset: %d plants\n\n", nrow(df)))

# =============================================================================
# Helper: capacity-weighted mean & median
# =============================================================================
wt_mean <- function(x, w) sum(x * w, na.rm = TRUE) / sum(w[!is.na(x)], na.rm = TRUE)
wt_med  <- function(x, w) {
  ok <- !is.na(x) & !is.na(w) & w > 0
  x  <- x[ok]; w  <- w[ok] / sum(w[ok])
  x[which.min(abs(cumsum(w[order(x)]) - 0.5))]
}

theme_base <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title   = element_text(size = 11, face = "bold"),
        plot.caption = element_text(size = 8, colour = "grey50"))

# =============================================================================
# T1 — Pre/post summary (descriptive)
# =============================================================================
cat("T1: Pre/post descriptive comparison (1 km ranks)...\n")

pre_df  <- df %>% filter(post2020 == 0)
post_df <- df %>% filter(post2020 == 1)

pre_mean  <- wt_mean(pre_df$pct_rank_1km,  pre_df$cap_mw)
post_mean <- wt_mean(post_df$pct_rank_1km, post_df$cap_mw)
delta     <- post_mean - pre_mean

wx <- wilcox.test(pre_df$pct_rank_1km, post_df$pct_rank_1km,
                  alternative = "greater")

cat(sprintf("  Pre-2020  cap-wt mean rank : %.2f\n", pre_mean))
cat(sprintf("  Post-2020 cap-wt mean rank : %.2f\n", post_mean))
cat(sprintf("  Delta (post - pre)         : %+.2f percentile points\n", delta))
cat(sprintf("  Wilcoxon test p-value      : %.4f\n\n", wx$p.value))

# =============================================================================
# T2 — Annual time series of within-province percentile rank
# =============================================================================
cat("T2: Annual time series...\n")

ann <- df %>%
  group_by(build_year) %>%
  summarise(
    wt_mean_rank = wt_mean(pct_rank_1km, cap_mw),
    wt_med_rank  = wt_med(pct_rank_1km,  cap_mw),
    total_cap_gw = sum(cap_mw) / 1e3,
    n_plants     = n(),
    .groups = "drop"
  )

# Fit separate OLS trends pre/post
pre_ann  <- ann %>% filter(build_year <= 2019)
post_ann <- ann %>% filter(build_year >= 2020)

m_pre  <- lm(wt_mean_rank ~ build_year, data = pre_ann,
              weights = pre_ann$total_cap_gw)
m_post <- lm(wt_mean_rank ~ build_year, data = post_ann,
              weights = post_ann$total_cap_gw)

pre_slope  <- round(coef(m_pre)["build_year"],  2)
post_slope <- round(coef(m_post)["build_year"], 2)

fig_t2 <- ggplot(ann, aes(x = build_year)) +
  # Background shading for post-2020
  annotate("rect", xmin = 2019.5, xmax = Inf,
           ymin = -Inf, ymax = Inf, fill = "#FEF9E7", alpha = 0.7) +
  annotate("text", x = 2020.2, y = max(ann$wt_mean_rank) + 1,
           label = "Post-shuangtan", hjust = 0, size = 3, colour = "grey40") +
  # Capacity-weighted mean rank
  geom_line(aes(y = wt_mean_rank), colour = "#2471A3", linewidth = 1) +
  geom_point(aes(y = wt_mean_rank, size = total_cap_gw),
             colour = "#2471A3", alpha = 0.8) +
  # Trend lines
  geom_smooth(data = pre_ann,  aes(y = wt_mean_rank), method = "lm",
              se = FALSE, linetype = "dashed", colour = "grey30",
              linewidth = 0.6) +
  geom_smooth(data = post_ann, aes(y = wt_mean_rank), method = "lm",
              se = FALSE, linetype = "dashed", colour = "#C0392B",
              linewidth = 0.6) +
  geom_hline(yintercept = 50, linetype = "dotted",
             colour = "grey50", linewidth = 0.5) +
  scale_x_continuous("Commissioning year", breaks = 2017:2024) +
  scale_y_continuous("Cap-wt mean within-province rank (0–100)", limits = c(25, 75)) +
  scale_size_continuous("Annual capacity\nadded (GW)", range = c(2, 8)) +
  theme_base +
  labs(
    title   = "T2 (1 km)  Annual within-province PVOUT percentile rank",
    caption = sprintf("Pre-2020 trend: %+.2f pp/yr.  Post-2020 trend: %+.2f pp/yr.  Horizontal dotted: provincial median (rank 50).",
                      pre_slope, post_slope)
  )

ggsave(file.path(OUT, "fig_11b_1_annual_rank_timeseries.png"),
       fig_t2, width = 8, height = 4.5, dpi = 300)
cat(sprintf("  Pre-2020 trend : %+.2f pp/yr\n", pre_slope))
cat(sprintf("  Post-2020 trend: %+.2f pp/yr\n\n", post_slope))

# =============================================================================
# T3 — Fixed-effects regression
# =============================================================================
cat("T3: Fixed-effects regression...\n")

df$province_f <- factor(df$prov_name)
df$year_c     <- df$build_year - 2020   # centred at cutoff

# Five models (identical structure to 11_siting_efficiency_regression.R)
m1 <- feols(pct_rank_1km ~ post2020,
            data = df, weights = ~cap_mw)
m2 <- feols(pct_rank_1km ~ post2020 | province_f,
            data = df, weights = ~cap_mw, cluster = ~province_f)
m3 <- feols(pct_rank_1km ~ post2020 + year_c | province_f,
            data = df, weights = ~cap_mw, cluster = ~province_f)
m4 <- feols(pct_rank_1km ~ post2020 + year_c + I(year_c^2) | province_f,
            data = df, weights = ~cap_mw, cluster = ~province_f)
m5 <- feols(pct_rank_1km ~ post2020 | province_f,    # unweighted robustness
            data = df, cluster = ~province_f)

coef_tbl <- lapply(list(m1 = m1, m2 = m2, m3 = m3, m4 = m4, m5 = m5), function(m) {
  s <- summary(m)$coeftable
  if (!"post2020" %in% rownames(s)) return(NULL)
  data.frame(
    coef  = s["post2020", "Estimate"],
    se    = s["post2020", "Std. Error"],
    pval  = s["post2020", "Pr(>|t|)"]
  )
}) %>% bind_rows(.id = "model") %>%
  mutate(lo = coef - 1.96 * se, hi = coef + 1.96 * se,
         sig = ifelse(pval < 0.01, "p<0.01",
               ifelse(pval < 0.05, "p<0.05",
               ifelse(pval < 0.10, "p<0.10", "n.s."))),
         label = c("No FE (cap-wt)", "Province FE (cap-wt)",
                   "+ year trend", "+ quadratic trend",
                   "Province FE (unweighted)"),
         model = factor(model, levels = rev(c("m1","m2","m3","m4","m5")),
                        labels = rev(c("No FE (cap-wt)", "Province FE (cap-wt)",
                                       "+ year trend", "+ quadratic trend",
                                       "Province FE (unweighted)"))))

cat("  Model coefficients on post2020:\n")
print(coef_tbl %>% select(model, coef, se, pval, sig), digits = 3)

fig_t3 <- ggplot(coef_tbl, aes(x = coef, y = model, colour = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_colour_manual(
    values = c("p<0.01" = "#1A5276", "p<0.05" = "#2471A3",
               "p<0.10" = "#7FB3D3", "n.s."   = "grey60"),
    name = "Significance") +
  theme_base +
  labs(title   = "T3 (1 km)  Post-2020 coefficient: within-province PVOUT rank",
       x       = "Estimated change in percentile rank (post vs pre 2020)",
       y       = NULL,
       caption = "Error bars: 95% CI. Province FE models: SEs clustered by province.")

ggsave(file.path(OUT, "fig_11b_3_regression_coefs.png"),
       fig_t3, width = 8, height = 4, dpi = 300)
cat("\n")

# =============================================================================
# T4 — Province-level: capacity growth vs change in siting rank
# =============================================================================
cat("T4: Province-level growth vs delta rank...\n")

prov_panel <- df %>%
  group_by(prov_name, post2020) %>%
  summarise(
    wt_rank = wt_mean(pct_rank_1km, cap_mw),
    cap_gw  = sum(cap_mw) / 1e3,
    n       = n(),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(names_from = post2020,
                     values_from = c(wt_rank, cap_gw, n),
                     names_prefix = "") %>%
  rename(rank_pre  = `wt_rank_0`, rank_post = `wt_rank_1`,
         cap_pre   = `cap_gw_0`,  cap_post  = `cap_gw_1`,
         n_pre     = `n_0`,       n_post    = `n_1`) %>%
  filter(!is.na(rank_pre), !is.na(rank_post),
         !is.na(cap_pre),  !is.na(cap_post)) %>%
  mutate(delta_rank   = rank_post - rank_pre,
         log_growth   = log(cap_post / cap_pre),
         total_cap_gw = cap_pre + cap_post)

cat(sprintf("  Provinces with both pre and post plants: %d\n", nrow(prov_panel)))

m_prov1 <- lm(delta_rank ~ log_growth, data = prov_panel,
               weights = total_cap_gw)
m_prov2 <- lm(delta_rank ~ log_growth + rank_pre, data = prov_panel,
               weights = total_cap_gw)

cat(sprintf("  Bivariate β (log_growth):            %+.3f  (p = %.3f)\n",
            coef(m_prov1)["log_growth"],
            summary(m_prov1)$coefficients["log_growth", "Pr(>|t|)"]))
cat(sprintf("  With rank_pre control β (log_growth):%+.3f  (p = %.3f)\n",
            coef(m_prov2)["log_growth"],
            summary(m_prov2)$coefficients["log_growth", "Pr(>|t|)"]))

fig_t4 <- ggplot(prov_panel, aes(x = log_growth, y = delta_rank,
                                  size = total_cap_gw)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
  geom_smooth(method = "lm", aes(weight = total_cap_gw),
              se = TRUE, colour = "#C0392B", fill = "#FADBD8",
              linewidth = 0.8) +
  geom_point(colour = "#2471A3", alpha = 0.75) +
  ggrepel::geom_text_repel(aes(label = prov_name), size = 2.5,
                            colour = "grey30", max.overlaps = 20) +
  scale_size_continuous("Total capacity\n(GW)", range = c(2, 10)) +
  theme_base +
  labs(
    title   = "T4 (1 km)  Capacity growth rate vs change in within-province siting rank",
    x       = "Log(post-2020 capacity / pre-2020 capacity)",
    y       = "Change in cap-wt mean rank (post - pre, pp)",
    caption = sprintf("Regression line (cap-weighted). β = %+.3f (bivariate), %+.3f (controlling for pre-2020 rank).",
                      coef(m_prov1)["log_growth"],
                      coef(m_prov2)["log_growth"])
  )

ggsave(file.path(OUT, "fig_11b_3_growth_vs_delta_rank.png"),
       fig_t4, width = 8, height = 5.5, dpi = 300)
cat("\n")

# =============================================================================
# 5. Province pre/post comparison chart
# =============================================================================
prov_plot <- prov_panel %>%
  filter(!is.na(rank_pre), !is.na(rank_post)) %>%
  arrange(desc(total_cap_gw)) %>%
  head(25) %>%
  mutate(prov_name = factor(prov_name, levels = prov_name[order(delta_rank)]))

fig_prov <- ggplot(prov_plot) +
  geom_segment(aes(x = rank_pre, xend = rank_post,
                   y = prov_name, yend = prov_name,
                   colour = delta_rank > 0),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               linewidth = 0.8) +
  geom_point(aes(x = rank_pre,  y = prov_name), shape = 1, size = 2,
             colour = "grey40") +
  geom_vline(xintercept = 50, linetype = "dotted", colour = "grey50",
             linewidth = 0.5) +
  scale_colour_manual(values = c("TRUE" = "#27AE60", "FALSE" = "#C0392B"),
                      labels = c("TRUE" = "Improved", "FALSE" = "Declined"),
                      name   = "Post-2020\nchange") +
  scale_x_continuous("Cap-wt within-province PVOUT rank (0–100)",
                     limits = c(20, 80), breaks = seq(20, 80, 10)) +
  theme_base +
  labs(title   = "T2b (1 km)  Within-province siting rank: pre vs post 2020 by province",
       y       = NULL,
       caption = "Circle = pre-2020 rank; arrowhead = post-2020 rank. Top 25 provinces by total capacity. Dotted: provincial median (rank 50).")

ggsave(file.path(OUT, "fig_11b_2_province_pre_post.png"),
       fig_prov, width = 8, height = 7, dpi = 300)

# =============================================================================
# 6. Direct comparison: 10 km vs 1 km ranks
# =============================================================================
cat("Building 10 km vs 1 km comparison figure...\n")

# Note: deploy_df has pct_rank from 11_siting_efficiency_regression.R if run;
# if not, we compare pvout values directly as a fallback
if ("pct_rank" %in% names(deploy_df)) {
  comp_df <- deploy_df %>%
    inner_join(deploy_1km %>% select(pct_rank_1km), by = character(0),
               copy = TRUE) %>%
    filter(!is.na(pct_rank), !is.na(pct_rank_1km))
  # Actually do it by row index since they share the same plants
  comp_df <- data.frame(
    pct_rank_10km = deploy_df$pct_rank,
    pct_rank_1km  = deploy_1km$pct_rank_1km,
    post2020      = deploy_df$post2020,
    cap_mw        = deploy_df$cap_mw
  ) %>% filter(!is.na(pct_rank_10km), !is.na(pct_rank_1km))

  fig_comp <- ggplot(comp_df, aes(pct_rank_10km, pct_rank_1km)) +
    geom_hex(bins = 60, alpha = 0.85) +
    geom_abline(slope = 1, intercept = 0, colour = "red",
                linetype = "dashed", linewidth = 0.5) +
    scale_fill_viridis_c(name = "N plants", trans = "log10") +
    theme_base +
    labs(title   = "Resolution comparison: 10 km vs 1 km within-province ranks",
         x       = "10 km percentile rank",
         y       = "1 km percentile rank",
         caption = "Diagonal: perfect agreement. Hex bins coloured by plant count.")

  ggsave(file.path(OUT, "fig_11b_compare_10km_vs_1km.png"),
         fig_comp, width = 6, height = 5.5, dpi = 300)
  cat("  Saved 10 km vs 1 km comparison figure.\n")
}

# =============================================================================
# Summary
# =============================================================================
cat("\n============================================================\n")
cat("SUMMARY: Within-province siting efficiency at 1 km resolution\n")
cat("============================================================\n")
cat(sprintf("Plants analysed           : %d\n", nrow(df)))
cat(sprintf("Pre-2020  mean rank (cw)  : %.2f\n", pre_mean))
cat(sprintf("Post-2020 mean rank (cw)  : %.2f\n", post_mean))
cat(sprintf("Delta (post - pre)        : %+.2f pp\n", delta))
cat(sprintf("Wilcoxon p (pre > post)   : %.4f\n", wx$p.value))
cat(sprintf("Annual trend pre-2020     : %+.2f pp/yr\n", pre_slope))
cat(sprintf("Annual trend post-2020    : %+.2f pp/yr\n", post_slope))
cat(sprintf("Province FE β (post2020)  : %+.3f pp  (Model m2)\n",
            coef_tbl$coef[coef_tbl$model == "Province FE (cap-wt)"]))
cat(sprintf("Province-level β (growth) : %+.3f (bivariate)\n",
            coef(m_prov1)["log_growth"]))
cat("============================================================\n")
cat("\nFigures saved to:", OUT, "\n")
