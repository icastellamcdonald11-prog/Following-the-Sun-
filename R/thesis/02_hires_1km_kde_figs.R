source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# 02_hires_1km_kde_figs.R  —  Thesis Project 4
# Reproduces Figures 6.1–6.3 at ~1 km resolution (0.00833° GSA grid → CEA 1 km)
# instead of the 10 km aggregated grid used in 01_siting_density_figs.R.
#
# PVOUT source : pvout_annual.tif (P3/processed) — original GSA grid, WGS84,
#                0.00833° (~926 m at 30°N), projected here to CEA 1 km.
# Topography   : topography_allowed_mask_10km.tif disaggregated to 1 km (near).
#                Note: this is a coarse approximation — each 10 km block is
#                uniformly marked feasible/infeasible across its 100 sub-cells.
#                Original 1 km SRTM was not found in the project.
# Population   : pop_density_10km_cea.tif disaggregated to 1 km (bilinear),
#                divided by 100 to recover persons/km².  Fine-scale variation
#                within each 10 km block is smoothed but totals are preserved.
# Plants       : PVOUT re-extracted from the 1 km raster at centroid locations.
#
# All intermediate rasters are cached to the Project 4 folder on first run.
# Subsequent runs skip the expensive projection / disaggregation steps.
#
# Requires: 00_setup.R to be sourced first (for upv_sf, deploy_df, etc.)
#           OR run standalone — this script will reload basemap + UPV data.
# =============================================================================

library(terra)
library(sf)
library(dplyr)
library(ggplot2)

P3  <- FTS_PROCESSED
P4  <- FTS_ROOT
OUT <- file.path(P4, "figures")
dir.create(OUT, showWarnings = FALSE)

CEA_CRS <- "+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# =============================================================================
# 1.  PVOUT at 1 km CEA
# =============================================================================
pvout_1km_cache <- file.path(P4, "data/derived/pvout_1km_cea.tif")

if (file.exists(pvout_1km_cache)) {
  cat("Loading cached 1 km PVOUT raster...\n")
  pvout_1km <- rast(pvout_1km_cache)
} else {
  cat("Projecting pvout_annual.tif to CEA at 1 km...\n")
  pvout_wgs <- rast(file.path(P3, "pvout_annual.tif"))

  # Load China nation boundary for crop + mask
  basemap_tmp   <- readRDS(file.path(P3, "china_basemap_bundle.rds"))
  china_cea_tmp <- st_transform(basemap_tmp$china_nation, CEA_CRS)

  pvout_1km <- project(pvout_wgs, CEA_CRS, res = 1000, method = "bilinear")
  pvout_1km <- crop(pvout_1km, vect(china_cea_tmp))
  pvout_1km <- mask(pvout_1km, vect(china_cea_tmp))

  writeRaster(pvout_1km, pvout_1km_cache, overwrite = TRUE)
  rm(pvout_wgs, basemap_tmp, china_cea_tmp); gc()
  cat("  Saved pvout_1km_cea.tif\n")
}

# =============================================================================
# 2.  Topography mask at 1 km CEA
#     Disaggregate 10 km mask (0/1) using "near" so cell labels stay integer.
#     Resample to match pvout_1km extent/alignment exactly.
# =============================================================================
topo_1km_cache <- file.path(P4, "data/derived/topo_mask_1km_cea.tif")

if (file.exists(topo_1km_cache)) {
  cat("Loading cached 1 km topography mask...\n")
  topo_1km <- rast(topo_1km_cache)
} else {
  cat("Disaggregating topography mask to 1 km...\n")
  topo_10km <- rast(file.path(P3, "topography_allowed_mask_10km.tif"))
  topo_1km  <- disagg(topo_10km, fact = 10, method = "near")
  topo_1km  <- resample(topo_1km, pvout_1km, method = "near")
  writeRaster(topo_1km, topo_1km_cache, overwrite = TRUE)
  rm(topo_10km); gc()
  cat("  Saved topo_mask_1km_cea.tif\n")
}

# =============================================================================
# 3.  Population density at 1 km CEA
#     Disaggregate the 10 km sum raster (divide by 100 → persons/km²).
#     Bilinear interpolation preserves smooth spatial gradient while keeping
#     aggregate totals correct to within rounding.
# =============================================================================
pop_1km_cache <- file.path(P4, "data/derived/pop_density_1km_cea.tif")

if (file.exists(pop_1km_cache)) {
  cat("Loading cached 1 km population raster...\n")
  pop_1km <- rast(pop_1km_cache)
} else {
  cat("Disaggregating population to 1 km...\n")
  pop_10km <- rast(file.path(P4, "data/derived/pop_density_10km_cea.tif"))
  # 10 km raster stores sum over 100 cells → divide to get persons/km²
  pop_1km  <- disagg(pop_10km / 100, fact = 10, method = "bilinear")
  pop_1km  <- resample(pop_1km, pvout_1km, method = "bilinear")
  pop_1km  <- clamp(pop_1km, lower = 0, values = FALSE)   # remove negative artefacts
  writeRaster(pop_1km, pop_1km_cache, overwrite = TRUE)
  rm(pop_10km); gc()
  cat("  Saved pop_density_1km_cea.tif\n")
}

# =============================================================================
# 4.  Extract raster value vectors
# =============================================================================
cat("Extracting raster value vectors (~32 M cells — may take ~30 s)...\n")

pvout_v <- values(pvout_1km, na.rm = FALSE)[, 1]
topo_v  <- values(topo_1km,  na.rm = FALSE)[, 1]
pop_v   <- values(pop_1km,   na.rm = FALSE)[, 1]

all_mask  <- is.finite(pvout_v) & pvout_v > 0
feas_mask <- all_mask & !is.na(topo_v) & topo_v == 1

all_land_pvout <- pvout_v[all_mask]
feas_pvout     <- pvout_v[feas_mask]
all_land_med   <- median(all_land_pvout)
feas_med       <- median(feas_pvout)

pop_all_mask  <- all_mask  & !is.na(pop_v) & pop_v >= 0
pop_feas_mask <- feas_mask & !is.na(pop_v) & pop_v >= 0

pop_all_pvout  <- pvout_v[pop_all_mask];  pop_all_wt  <- pop_v[pop_all_mask]
pop_feas_pvout <- pvout_v[pop_feas_mask]; pop_feas_wt <- pop_v[pop_feas_mask]
pop_all_wt[pop_all_wt   == 0] <- 1e-6
pop_feas_wt[pop_feas_wt == 0] <- 1e-6

cat(sprintf("  All-land cells : %s\n",  format(sum(all_mask),  big.mark = ",")))
cat(sprintf("  Feasible cells : %s\n",  format(sum(feas_mask), big.mark = ",")))

# =============================================================================
# 5.  Re-extract plant PVOUT from the 1 km raster
# =============================================================================
cat("Extracting plant PVOUT from 1 km raster...\n")

if (!exists("upv_sf")) {
  sf_use_s2(FALSE)
  upv_sf <- readRDS(file.path(P3, "UPV_core_GRW_2024Q2_ge1MW.rds"))
  upv_sf <- st_make_valid(upv_sf)
}

plant_cents   <- st_centroid(upv_sf)
pvout_ext_1km <- terra::extract(pvout_1km, vect(plant_cents), method = "simple")

dep_pvout <- as.numeric(pvout_ext_1km[[2]])
dep_cap   <- as.numeric(upv_sf$capacity_mw_used)
ok        <- is.finite(dep_pvout) & is.finite(dep_cap) & dep_cap > 0
dep_pvout <- dep_pvout[ok]
dep_cap   <- dep_cap[ok]
dep_w     <- dep_cap / sum(dep_cap)

dep_med_unw   <- median(dep_pvout)
dep_mean_capw <- sum(dep_pvout * dep_w)

cat(sprintf("  Plants with valid 1 km PVOUT : %d / %d\n", sum(ok), nrow(upv_sf)))
cat(sprintf("  Unweighted plant median      : %.1f kWh/kWp/yr\n", dep_med_unw))
cat(sprintf("  Capacity-weighted mean       : %.1f kWh/kWp/yr\n", dep_mean_capw))

# =============================================================================
# 6.  KDE computation
#     For land vectors >500 k cells, draw a stratified random subsample so
#     density() remains fast. The subsample is large enough that the KDE
#     shape is indistinguishable from the full-population curve.
# =============================================================================
MAXN <- 500000L

subsamp_kde <- function(x, w = NULL, n = MAXN) {
  # Subsample uniformly; rescale weights on the subsample.
  if (length(x) > n) {
    idx <- sample.int(length(x), n)
    x   <- x[idx]
    if (!is.null(w)) w <- w[idx]
  }
  list(x = x, w = w)
}

wkde <- function(x, weights = NULL, from, to, bw = "nrd0", n = 2048) {
  if (length(x) < 2) stop("Need \u22652 observations for KDE")
  if (!is.null(weights)) weights <- weights / sum(weights)
  k <- density(x, weights = weights, bw = bw, from = from, to = to, n = n)
  data.frame(pvout = k$x, dens = k$y)
}

# X range — clip tails at 0.05th / 99.95th percentiles to avoid edge artefacts
xfrom_61 <- floor(quantile(all_land_pvout, 0.0005) / 50) * 50
xto_61   <- ceiling(quantile(all_land_pvout, 0.9995) / 50) * 50
xfrom    <- floor(quantile(feas_pvout, 0.0005) / 50) * 50
xto      <- ceiling(quantile(feas_pvout, 0.9995) / 50) * 50

cat("Computing KDEs...\n")

s_al  <- subsamp_kde(all_land_pvout)
s_fe  <- subsamp_kde(feas_pvout)
s_pa  <- subsamp_kde(pop_all_pvout,  pop_all_wt)
s_pf  <- subsamp_kde(pop_feas_pvout, pop_feas_wt)

kde_land_all  <- wkde(s_al$x, from = xfrom_61, to = xto_61)
kde_land_feas <- wkde(s_fe$x, from = xfrom,    to = xto)
kde_dep_unw   <- wkde(dep_pvout,               from = xfrom_61, to = xto_61)
kde_dep_capw  <- wkde(dep_pvout, weights = dep_w, from = xfrom, to = xto)
kde_pop_all   <- wkde(s_pa$x, weights = s_pa$w, from = xfrom_61, to = xto_61)
kde_pop_feas  <- wkde(s_pf$x, weights = s_pf$w, from = xfrom,    to = xto)

cat("  KDEs computed.\n")

# =============================================================================
# 7.  Shared theme and colour palette  (identical to 01_siting_density_figs.R)
# =============================================================================
col_land <- "grey70"
col_pop  <- "#1A3A5C"
col_dep  <- "#C0392B"
lwd_main <- 0.9
lwd_pop  <- 0.75

theme_dens <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor    = element_blank(),
    panel.grid.major.y  = element_blank(),
    axis.line.x   = element_line(colour = "grey30", linewidth = 0.3),
    axis.ticks.x  = element_line(colour = "grey30", linewidth = 0.3),
    axis.title    = element_text(size = 10),
    legend.position   = "bottom",
    legend.key.width  = unit(1.1, "cm"),
    legend.key.height = unit(0.35, "cm"),
    legend.text       = element_text(size = 8.5),
    legend.spacing.x  = unit(0.3, "cm"),
    plot.title    = element_text(size = 11, face = "bold"),
    plot.caption  = element_text(size = 8, colour = "grey50")
  )

xsc <- function(from, to)
  scale_x_continuous(name = "PVOUT (kWh/kWp/year)",
                     limits = c(from, to),
                     breaks = pretty(c(from, to), n = 8),
                     expand = c(0, 0))

ysc <- scale_y_continuous(name   = "Density",
                           expand = expansion(mult = c(0, 0.05)))

res_note <- "1 km resolution (GSA ~926 m grid projected to CEA 1 km)."

# =============================================================================
# Figure 6.1  — All land vs unweighted plants
# =============================================================================
cat("Building Figure 6.1 (1 km)...\n")

ymax_61 <- max(kde_land_all$dens, kde_dep_unw$dens, kde_pop_all$dens) * 1.08

fig_6_1 <- ggplot() +
  geom_area(data = kde_land_all, aes(x = pvout, y = dens),
            fill = col_land, alpha = 0.6, colour = NA) +
  geom_line(data = kde_land_all,
            aes(x = pvout, y = dens, colour = "All land PVOUT"),
            linewidth = lwd_main) +
  geom_line(data = kde_pop_all,
            aes(x = pvout, y = dens, colour = "Population-weighted land"),
            linewidth = lwd_pop, linetype = "dashed") +
  geom_line(data = kde_dep_unw,
            aes(x = pvout, y = dens, colour = "All UPV plants (unweighted)"),
            linewidth = lwd_main) +
  geom_vline(xintercept = all_land_med, linetype = "dashed",
             colour = "grey20", linewidth = 0.5) +
  geom_vline(xintercept = dep_med_unw,  linetype = "dashed",
             colour = col_dep, linewidth = 0.5) +
  annotate("text", x = all_land_med + 15, y = ymax_61 * 0.97,
           label = "Land\nmedian", size = 2.5, colour = "grey20",
           hjust = 0, vjust = 1, lineheight = 0.85) +
  annotate("text", x = dep_med_unw  + 15, y = ymax_61 * 0.97,
           label = "Plant\nmedian", size = 2.5, colour = col_dep,
           hjust = 0, vjust = 1, lineheight = 0.85) +
  scale_colour_manual(
    values = c("All land PVOUT"              = col_land,
               "Population-weighted land"    = col_pop,
               "All UPV plants (unweighted)" = col_dep),
    name  = NULL,
    guide = guide_legend(nrow = 1, override.aes = list(
      linetype  = c("solid", "dashed", "solid"),
      linewidth = c(lwd_main, lwd_pop, lwd_main)))) +
  xsc(xfrom_61, xto_61) + ysc +
  coord_cartesian(ylim = c(0, ymax_61)) +
  theme_dens +
  labs(
    title   = "Figure 6.1  China land PVOUT vs deployed UPV plants (unweighted)",
    caption = paste0(
      "Each plant counted equally regardless of size (n = ", sum(ok), " plants). ",
      "Dashed lines: all-land median and plant median.\n",
      "Navy dashed: population-weighted land PVOUT (WorldPop 2020). ", res_note)
  )

ggsave(file.path(OUT, "fig_6_1_1km_kde_with_pop.png"),
       fig_6_1, width = 8, height = 4.5, dpi = 300)
cat("  Saved fig_6_1_1km_kde_with_pop.png\n")

# =============================================================================
# Figure 6.2  — Feasible land vs capacity-weighted plants
# =============================================================================
cat("Building Figure 6.2 (1 km)...\n")

ymax_62 <- max(kde_land_feas$dens, kde_dep_capw$dens, kde_pop_feas$dens) * 1.08

fig_6_2 <- ggplot() +
  geom_area(data = kde_land_feas, aes(x = pvout, y = dens),
            fill = col_land, alpha = 0.6, colour = NA) +
  geom_line(data = kde_land_feas,
            aes(x = pvout, y = dens, colour = "Feasible land PVOUT"),
            linewidth = lwd_main) +
  geom_line(data = kde_pop_feas,
            aes(x = pvout, y = dens, colour = "Population-weighted (feasible land)"),
            linewidth = lwd_pop, linetype = "dashed") +
  geom_line(data = kde_dep_capw,
            aes(x = pvout, y = dens, colour = "Deployed UPV (capacity-weighted)"),
            linewidth = lwd_main) +
  geom_vline(xintercept = feas_med,      linetype = "dashed",
             colour = "grey20", linewidth = 0.5) +
  geom_vline(xintercept = dep_mean_capw, linetype = "dotted",
             colour = col_dep,  linewidth = 0.5) +
  annotate("text", x = feas_med      + 15, y = ymax_62 * 0.97,
           label = "Feasible\nmedian", size = 2.5, colour = "grey20",
           hjust = 0, vjust = 1, lineheight = 0.85) +
  annotate("text", x = dep_mean_capw + 15, y = ymax_62 * 0.84,
           label = "Cap-wt\nmean", size = 2.5, colour = col_dep,
           hjust = 0, vjust = 1, lineheight = 0.85) +
  scale_colour_manual(
    values = c("Feasible land PVOUT"                 = col_land,
               "Population-weighted (feasible land)" = col_pop,
               "Deployed UPV (capacity-weighted)"    = col_dep),
    name  = NULL,
    guide = guide_legend(nrow = 1, override.aes = list(
      linetype  = c("solid", "dashed", "solid"),
      linewidth = c(lwd_main, lwd_pop, lwd_main)))) +
  xsc(xfrom, xto) + ysc +
  coord_cartesian(ylim = c(0, ymax_62)) +
  theme_dens +
  labs(
    title   = "Figure 6.2  Feasible land PVOUT vs deployed UPV plants (capacity-weighted)",
    caption = paste0(
      "Plants weighted by installed capacity (MW). Total: ",
      round(sum(dep_cap) / 1e3, 1), " GW.\n",
      "Dashed: feasible-land median. Dotted: capacity-weighted mean. ",
      "Navy dashed: population-weighted feasible land PVOUT (WorldPop 2020). ", res_note)
  )

ggsave(file.path(OUT, "fig_6_2_1km_capweighted_kde_with_pop.png"),
       fig_6_2, width = 8, height = 4.5, dpi = 300)
cat("  Saved fig_6_2_1km_capweighted_kde_with_pop.png\n")

# =============================================================================
# Figure 6.3  — Feasible land (uniform grey) vs cap-weighted UPV ≥1 MW
# =============================================================================
cat("Building Figure 6.3 (1 km)...\n")

ymax_63 <- max(kde_land_feas$dens, kde_dep_capw$dens, kde_pop_feas$dens) * 1.08

fig_6_3 <- ggplot() +
  geom_area(data = kde_land_feas, aes(x = pvout, y = dens),
            fill = col_land, alpha = 0.6, colour = NA) +
  geom_line(data = kde_land_feas,
            aes(x = pvout, y = dens, colour = "Feasible land PVOUT"),
            linewidth = lwd_main) +
  geom_line(data = kde_pop_feas,
            aes(x = pvout, y = dens, colour = "Population-weighted (feasible land)"),
            linewidth = lwd_pop, linetype = "dashed") +
  geom_line(data = kde_dep_capw,
            aes(x = pvout, y = dens, colour = "UPV \u22651 MW (capacity-weighted)"),
            linewidth = lwd_main) +
  geom_vline(xintercept = feas_med,      linetype = "dashed",
             colour = "grey20", linewidth = 0.5) +
  geom_vline(xintercept = dep_mean_capw, linetype = "dotted",
             colour = col_dep,  linewidth = 0.5) +
  annotate("text", x = feas_med      + 15, y = ymax_63 * 0.97,
           label = "Feasible\nmedian", size = 2.5, colour = "grey20",
           hjust = 0, vjust = 1, lineheight = 0.85) +
  annotate("text", x = dep_mean_capw + 15, y = ymax_63 * 0.84,
           label = "Cap-wt\nmean", size = 2.5, colour = col_dep,
           hjust = 0, vjust = 1, lineheight = 0.85) +
  scale_colour_manual(
    values = c("Feasible land PVOUT"                 = "grey40",
               "Population-weighted (feasible land)" = col_pop,
               "UPV \u22651 MW (capacity-weighted)"  = col_dep),
    name  = NULL,
    guide = guide_legend(nrow = 1, override.aes = list(
      linetype  = c("solid", "dashed", "solid"),
      linewidth = c(lwd_main, lwd_pop, lwd_main)))) +
  xsc(xfrom, xto) + ysc +
  coord_cartesian(ylim = c(0, ymax_63)) +
  theme_dens +
  labs(
    title   = "Figure 6.3  UPV \u22651 MW distribution vs feasible land PVOUT",
    caption = paste0(
      "Feasible land: elevation \u22645000 m, slope \u226450\u00b0. ",
      "Dashed: feasible-land median. Dotted: capacity-weighted mean.\n",
      "Navy dashed: population-weighted feasible land PVOUT (WorldPop 2020). ", res_note)
  )

ggsave(file.path(OUT, "fig_6_3_1km_upv_feasible_kde_with_pop.png"),
       fig_6_3, width = 8, height = 4.5, dpi = 300)
cat("  Saved fig_6_3_1km_upv_feasible_kde_with_pop.png\n")

# =============================================================================
# Summary statistics
# =============================================================================
cat("\nAll three 1 km figures saved to:", OUT, "\n\n")
cat("Key statistics (1 km resolution):\n")
cat(sprintf("  All-land median PVOUT      : %.1f kWh/kWp/yr\n", all_land_med))
cat(sprintf("  Feasible-land median PVOUT : %.1f kWh/kWp/yr\n", feas_med))
cat(sprintf("  Plant median PVOUT (unwt)  : %.1f kWh/kWp/yr\n", dep_med_unw))
cat(sprintf("  Plant mean PVOUT (cap-wt)  : %.1f kWh/kWp/yr\n", dep_mean_capw))
cat(sprintf("  Share below all-land median: %.1f%%\n",
            100 * mean(dep_pvout < all_land_med)))
