source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# 08_cpvpd_validation.R  —  Thesis Project 4
# Cross-validation of GRW density curves using CPVPD-2024
# (Sun et al. 2025, Scientific Data, DOI: 10.5281/zenodo.15618227)
#
# Produces four figures:
#   Fig C1  — All land vs CPVPD all plants unweighted       (mirrors Fig 6.1)
#   Fig C2  — All land vs CPVPD all plants cap-weighted     (mirrors Fig 6.2)
#   Fig C3  — Feasible land vs CPVPD UPV cap-weighted       (mirrors Fig 6.3)
#   Fig C4  — GRW vs CPVPD overlay: UPV cap-wt on feasible land (comparison)
#
# Speed strategy: read DBF attribute tables only (no polygon geometry),
# using pre-computed lat/lon centroids already in the CPVPD dataset.
# Cached to cpvpd_processed.rds after first run.
#
# Requires: 00_setup.R already sourced (deploy_df, grid_df,
#           pvout_national_median, grid_df$pop_density)
# =============================================================================

library(ggplot2)
library(dplyr)
library(terra)
library(foreign)   # read.dbf()

P4      <- FTS_ROOT
P3_PROC <- FTS_PROCESSED
SHP_DIR <- file.path(FTS_ROOT, "data", "raw", "CPVPD-2024_4326")
OUT     <- file.path(P4, "figures")
CACHE   <- file.path(P4, "data/derived/cpvpd_processed.rds")
dir.create(OUT, showWarnings = FALSE)

# =============================================================================
# 1. Load CPVPD attributes (DBF only — no polygon geometry needed)
#    The Area column (m²) in the DBF IS the polygon area; no geometry required.
# =============================================================================
if (file.exists(CACHE)) {
  cat("Loading cached CPVPD data...\n")
  cpvpd <- readRDS(CACHE)
  cat(sprintf("  %d plants loaded from cache\n", nrow(cpvpd)))
} else {
  cat("Reading CPVPD DBF files...\n")
  exclude   <- c("taiwan", "hong_kong", "macao")
  dbf_files <- list.files(SHP_DIR, pattern = "[.]dbf$", full.names = TRUE)
  dbf_files <- dbf_files[!grepl(paste(exclude, collapse = "|"), basename(dbf_files))]
  cat(sprintf("  Reading %d provincial DBF files...\n", length(dbf_files)))

  dfs <- lapply(dbf_files, function(f) {
    prov <- sub("[.]dbf$", "", basename(f))
    d    <- read.dbf(f, as.is = TRUE)
    d$province_en <- prov
    d
  })
  cpvpd <- do.call(rbind, dfs)
  cat(sprintf("  Total rows: %d\n", nrow(cpvpd)))
  cat(sprintf("  Columns: %s\n", paste(names(cpvpd), collapse = ", ")))

  PD_2024        <- 65.90          # MW/km²
  cpvpd$area_km2 <- cpvpd$Area / 1e6
  cpvpd$cap_mw   <- cpvpd$area_km2 * PD_2024

  cat("Extracting PVOUT from 10 km raster...\n")
  pvout_r <- rast(file.path(P3_PROC, "pvout_annual_10km.tif"))
  pts     <- vect(cbind(cpvpd$longitude, cpvpd$latitude), crs = "EPSG:4326")
  pts_t   <- project(pts, crs(pvout_r))
  ext     <- terra::extract(pvout_r, pts_t)
  cpvpd$pvout <- ext[[2]]

  valid <- is.finite(cpvpd$pvout) & cpvpd$pvout > 0
  cat(sprintf("  Valid PVOUT: %d / %d plants\n", sum(valid), nrow(cpvpd)))
  cpvpd <- cpvpd[valid, ]

  saveRDS(cpvpd, CACHE)
  cat(sprintf("  Cached to %s\n", CACHE))
}

# =============================================================================
# 2. Filters and capacity weights
# =============================================================================
PD_2024     <- 65.90
AREA_1MW_M2 <- 1 / PD_2024 * 1e6   # ~15,190 m²

cpvpd_all <- cpvpd
cpvpd_upv <- cpvpd[cpvpd$Area >= AREA_1MW_M2, ]   # >=1 MW

# Ensure cap_mw exists (in case loaded from cache without it)
if (!"cap_mw" %in% names(cpvpd_all)) {
  cpvpd_all$cap_mw <- (cpvpd_all$Area / 1e6) * PD_2024
  cpvpd_upv$cap_mw <- (cpvpd_upv$Area / 1e6) * PD_2024
}

w_all <- cpvpd_all$cap_mw / sum(cpvpd_all$cap_mw)
w_upv <- cpvpd_upv$cap_mw / sum(cpvpd_upv$cap_mw)

# GRW reference vectors (from 00_setup.R deploy_df)
ok_dep    <- is.finite(deploy_df$pvout) & is.finite(deploy_df$cap_mw) & deploy_df$cap_mw > 0
dep_pvout <- deploy_df$pvout[ok_dep]
dep_cap   <- deploy_df$cap_mw[ok_dep]
dep_w     <- dep_cap / sum(dep_cap)

cat(sprintf("\nCPVPD all plants:  %d  (total cap: %.1f GW)\n",
            nrow(cpvpd_all), sum(cpvpd_all$cap_mw) / 1e3))
cat(sprintf("CPVPD UPV >=1MW:   %d\n", nrow(cpvpd_upv)))
cat(sprintf("GRW  UPV >=1MW:    %d  (total cap: %.1f GW)\n",
            length(dep_pvout), sum(dep_cap) / 1e3))

# =============================================================================
# 3. Summary statistics
# =============================================================================
wmed <- function(x, w) {
  w <- w / sum(w); o <- order(x)
  x[o][which.min(abs(cumsum(w[o]) - 0.5))]
}

# CPVPD scalars
cpvpd_unw_med      <- median(cpvpd_all$pvout)
cpvpd_unw_upv_med  <- median(cpvpd_upv$pvout)
cpvpd_capw_mean_all <- sum(cpvpd_all$pvout * w_all)
cpvpd_capw_mean_upv <- sum(cpvpd_upv$pvout * w_upv)
cpvpd_capw_med_upv  <- wmed(cpvpd_upv$pvout, w_upv)

# GRW scalars
grw_unw_med        <- median(dep_pvout)
grw_capw_mean      <- sum(dep_pvout * dep_w)
grw_capw_med       <- wmed(dep_pvout, dep_w)

cat("\n=== KEY STATISTICS: CPVPD vs GRW ===\n")
cat(sprintf("%-44s  %6s  %8s  %7s\n", "", "Median", "Cap-wt mean", "Pct<1400"))
cat(sprintf("%-44s  %6.0f  %11.0f  %6.1f%%\n", "CPVPD all, unweighted",
            cpvpd_unw_med, mean(cpvpd_all$pvout), 100 * mean(cpvpd_all$pvout < 1400)))
cat(sprintf("%-44s  %6.0f  %11.0f  %6.1f%%\n", "CPVPD all, cap-weighted",
            wmed(cpvpd_all$pvout, w_all), cpvpd_capw_mean_all,
            100 * sum(w_all[cpvpd_all$pvout < 1400])))
cat(sprintf("%-44s  %6.0f  %11.0f  %6.1f%%\n", "CPVPD UPV >=1MW, cap-weighted",
            cpvpd_capw_med_upv, cpvpd_capw_mean_upv,
            100 * sum(w_upv[cpvpd_upv$pvout < 1400])))
cat(sprintf("%-44s  %6.0f  %11.0f  %6.1f%%\n", "GRW  UPV >=1MW, unweighted",
            grw_unw_med, mean(dep_pvout), 100 * mean(dep_pvout < 1400)))
cat(sprintf("%-44s  %6.0f  %11.0f  %6.1f%%\n", "GRW  UPV >=1MW, cap-weighted",
            grw_capw_med, grw_capw_mean, 100 * sum(dep_w[dep_pvout < 1400])))

# =============================================================================
# 4. Land and population vectors (from 00_setup.R grid_df)
# =============================================================================

# All land
all_land_pvout <- grid_df$pvout_annual[is.finite(grid_df$pvout_annual) &
                                         grid_df$pvout_annual > 0]
all_land_med   <- median(all_land_pvout)

# Feasible land
feas_mask      <- grid_df$feasible & is.finite(grid_df$pvout_annual) &
                  grid_df$pvout_annual > 0
feas_pvout     <- grid_df$pvout_annual[feas_mask]
feas_med       <- pvout_national_median

# Population-weighted — all land
pop_all_mask   <- is.finite(grid_df$pvout_annual) & grid_df$pvout_annual > 0 &
                  !is.na(grid_df$pop_density) & grid_df$pop_density >= 0
pop_all_pvout  <- grid_df$pvout_annual[pop_all_mask]
pop_all_wt     <- grid_df$pop_density[pop_all_mask]
pop_all_wt[pop_all_wt == 0] <- 1e-6

# Population-weighted — feasible land
pop_feas_mask  <- feas_mask & !is.na(grid_df$pop_density) & grid_df$pop_density >= 0
pop_feas_pvout <- grid_df$pvout_annual[pop_feas_mask]
pop_feas_wt    <- grid_df$pop_density[pop_feas_mask]
pop_feas_wt[pop_feas_wt == 0] <- 1e-6

# =============================================================================
# 5. KDE helper and x ranges
# =============================================================================
wkde <- function(x, weights = NULL, from, to, bw = "nrd0", n = 2048) {
  if (length(x) < 2) stop("Need >=2 observations for KDE")
  if (!is.null(weights)) weights <- weights / sum(weights)
  k <- density(x, weights = weights, bw = bw, from = from, to = to, n = n)
  data.frame(pvout = k$x, dens = k$y)
}

# Wide range: all land + all cpvpd (for C1, C2)
xlim_w  <- range(c(all_land_pvout, cpvpd_all$pvout))
xfrom_w <- floor(xlim_w[1]   / 50) * 50
xto_w   <- ceiling(xlim_w[2] / 50) * 50

# Narrow range: feasible land + UPV plants (for C3, C4)
xlim_n  <- range(c(feas_pvout, cpvpd_upv$pvout, dep_pvout))
xfrom_n <- floor(xlim_n[1]   / 50) * 50
xto_n   <- ceiling(xlim_n[2] / 50) * 50

# Pre-compute all KDEs
kde_land_all_w  <- wkde(all_land_pvout,  from = xfrom_w, to = xto_w)
kde_land_feas_n <- wkde(feas_pvout,      from = xfrom_n, to = xto_n)
kde_pop_all_w   <- wkde(pop_all_pvout,   weights = pop_all_wt,  from = xfrom_w, to = xto_w)
kde_pop_feas_n  <- wkde(pop_feas_pvout,  weights = pop_feas_wt, from = xfrom_n, to = xto_n)

kde_c_unw_w  <- wkde(cpvpd_all$pvout, from = xfrom_w, to = xto_w)
kde_c_capw_w <- wkde(cpvpd_all$pvout, weights = w_all, from = xfrom_w, to = xto_w)
kde_c_upv_n  <- wkde(cpvpd_upv$pvout, weights = w_upv, from = xfrom_n, to = xto_n)
kde_g_capw_n <- wkde(dep_pvout,        weights = dep_w, from = xfrom_n, to = xto_n)

# Two-tone shading for feasible land
kde_feas_below <- kde_land_feas_n
kde_feas_above <- kde_land_feas_n
kde_feas_below$dens[kde_feas_below$pvout > feas_med] <- 0
kde_feas_above$dens[kde_feas_above$pvout < feas_med] <- 0

# =============================================================================
# 6. Shared theme and scale helpers (matching 01_siting_density_figs.R exactly)
# =============================================================================
col_land <- "grey70"
col_pop  <- "#1A3A5C"    # navy
col_dep  <- "#C0392B"    # red  (CPVPD)
col_grw  <- "#2471A3"    # blue (GRW, used in C4 only)

lwd_main <- 0.9
lwd_pop  <- 0.75

theme_dens <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.line.x  = element_line(colour = "grey30", linewidth = 0.3),
    axis.ticks.x = element_line(colour = "grey30", linewidth = 0.3),
    axis.title   = element_text(size = 10),
    legend.position   = "bottom",
    legend.key.width  = unit(1.1, "cm"),
    legend.key.height = unit(0.35, "cm"),
    legend.text       = element_text(size = 8.5),
    legend.spacing.x  = unit(0.3, "cm"),
    plot.title   = element_text(size = 11, face = "bold"),
    plot.caption = element_text(size = 8, colour = "grey50")
  )

xsc <- function(from, to)
  scale_x_continuous("PVOUT (kWh/kWp/year)", limits = c(from, to),
                     breaks = pretty(c(from, to), n = 8), expand = c(0, 0))

ysc <- scale_y_continuous("Density", expand = expansion(mult = c(0, 0.05)))

# =============================================================================
# Fig C1 — All land vs CPVPD all plants (unweighted)   [mirrors Fig 6.1]
# =============================================================================
cat("\nBuilding Fig C1 (unweighted, all plants — mirrors 6.1)...\n")

ymax_c1 <- max(kde_land_all_w$dens, kde_c_unw_w$dens, kde_pop_all_w$dens) * 1.08

fig_c1 <- ggplot() +
  geom_area(data = kde_land_all_w,
            aes(pvout, dens), fill = col_land, alpha = 0.6, colour = NA) +
  geom_line(data = kde_land_all_w,
            aes(pvout, dens, colour = "All land PVOUT"),
            linewidth = lwd_main) +
  geom_line(data = kde_pop_all_w,
            aes(pvout, dens, colour = "Population-weighted land"),
            linewidth = lwd_pop, linetype = "dashed") +
  geom_line(data = kde_c_unw_w,
            aes(pvout, dens, colour = "All CPVPD plants (unweighted)"),
            linewidth = lwd_main) +
  geom_vline(xintercept = all_land_med,   linetype = "dashed",
             colour = "grey20", linewidth = 0.5) +
  geom_vline(xintercept = cpvpd_unw_med, linetype = "dashed",
             colour = col_dep,  linewidth = 0.5) +
  annotate("text", x = all_land_med   + 15, y = ymax_c1 * 0.97,
           label = "Land\nmedian",  size = 2.5, colour = "grey20",
           hjust = 0, vjust = 1, lineheight = 0.85) +
  annotate("text", x = cpvpd_unw_med + 15, y = ymax_c1 * 0.97,
           label = "Plant\nmedian", size = 2.5, colour = col_dep,
           hjust = 0, vjust = 1, lineheight = 0.85) +
  scale_colour_manual(
    values = c("All land PVOUT"                = col_land,
               "Population-weighted land"      = col_pop,
               "All CPVPD plants (unweighted)" = col_dep),
    name = NULL,
    guide = guide_legend(nrow = 1, override.aes = list(
      linetype  = c("solid", "dashed", "solid"),
      linewidth = c(lwd_main, lwd_pop, lwd_main)))) +
  xsc(xfrom_w, xto_w) + ysc +
  coord_cartesian(ylim = c(0, ymax_c1)) +
  theme_dens +
  labs(
    title   = "Figure C1  China land PVOUT vs all CPVPD plants (unweighted)",
    caption = paste0(
      "All CPVPD-2024 polygons (no size filter), each plant counted equally ",
      "(n = ", nrow(cpvpd_all), " plants). ",
      "Dashed lines: all-land median and plant median.\n",
      "Navy dashed: population-weighted land PVOUT (WorldPop 2020)."
    )
  )

ggsave(file.path(OUT, "fig_c1_cpvpd_unweighted.png"),
       fig_c1, width = 8, height = 4.5, dpi = 300)
cat("  Saved fig_c1_cpvpd_unweighted.png\n")

# =============================================================================
# Fig C2 — All land vs CPVPD all plants (capacity-weighted)  [mirrors Fig 6.2]
# =============================================================================
cat("Building Fig C2 (cap-weighted, all land — mirrors 6.2)...\n")

ymax_c2 <- max(kde_land_all_w$dens, kde_c_capw_w$dens, kde_pop_all_w$dens) * 1.08

fig_c2 <- ggplot() +
  geom_area(data = kde_land_all_w,
            aes(pvout, dens), fill = col_land, alpha = 0.55, colour = NA) +
  geom_line(data = kde_land_all_w,
            aes(pvout, dens, colour = "All land PVOUT (unfiltered)"),
            linewidth = lwd_main) +
  geom_line(data = kde_pop_all_w,
            aes(pvout, dens, colour = "Population-weighted land"),
            linewidth = lwd_pop, linetype = "dashed") +
  geom_line(data = kde_c_capw_w,
            aes(pvout, dens, colour = "All CPVPD plants (capacity-weighted)"),
            linewidth = lwd_main) +
  geom_vline(xintercept = cpvpd_unw_med,      linetype = "dashed",
             colour = col_dep,  linewidth = 0.5) +
  geom_vline(xintercept = all_land_med,        linetype = "dashed",
             colour = "grey20", linewidth = 0.5) +
  geom_vline(xintercept = cpvpd_capw_mean_all, linetype = "dotted",
             colour = col_dep,  linewidth = 0.5) +
  annotate("text", x = cpvpd_unw_med      + 12, y = ymax_c2 * 0.97,
           label = "Plant\nmedian",   size = 2.5, colour = col_dep,
           hjust = 0, vjust = 1, lineheight = 0.85) +
  annotate("text", x = all_land_med       + 12, y = ymax_c2 * 0.85,
           label = "Land\nmedian",    size = 2.5, colour = "grey20",
           hjust = 0, vjust = 1, lineheight = 0.85) +
  annotate("text", x = cpvpd_capw_mean_all + 12, y = ymax_c2 * 0.73,
           label = "Cap-wt\nmean",    size = 2.5, colour = col_dep,
           hjust = 0, vjust = 1, lineheight = 0.85) +
  scale_colour_manual(
    values = c("All land PVOUT (unfiltered)"             = col_land,
               "Population-weighted land"                = col_pop,
               "All CPVPD plants (capacity-weighted)"    = col_dep),
    name = NULL,
    guide = guide_legend(nrow = 1, override.aes = list(
      linetype  = c("solid", "dashed", "solid"),
      linewidth = c(lwd_main, lwd_pop, lwd_main)))) +
  xsc(xfrom_w, xto_w) + ysc +
  coord_cartesian(ylim = c(0, ymax_c2)) +
  theme_dens +
  labs(
    title   = "Figure C2  China solar siting vs solar resource (CPVPD, capacity-weighted)",
    caption = paste0(
      "CPVPD-2024 polygons weighted by installed capacity (MW). Total: ",
      round(sum(cpvpd_all$cap_mw) / 1e3, 1), " GW across ", nrow(cpvpd_all), " plants.\n",
      "Dashed red: unweighted plant median (", round(cpvpd_unw_med), " kWh/kWp/yr). ",
      "Dashed grey: all-land median. Dotted red: cap-weighted mean (",
      round(cpvpd_capw_mean_all), "). Navy dashed: pop-weighted land PVOUT (WorldPop 2020)."
    )
  )

ggsave(file.path(OUT, "fig_c2_cpvpd_capweighted_allland.png"),
       fig_c2, width = 8, height = 4.5, dpi = 300)
cat("  Saved fig_c2_cpvpd_capweighted_allland.png\n")

# =============================================================================
# Fig C3 — Feasible land (split shading) vs CPVPD UPV cap-weighted
#           [mirrors Fig 6.3]
# =============================================================================
cat("Building Fig C3 (UPV cap-weighted vs feasible land — mirrors 6.3)...\n")

ymax_c3 <- max(kde_land_feas_n$dens, kde_c_upv_n$dens, kde_pop_feas_n$dens) * 1.08

fig_c3 <- ggplot() +
  geom_area(data = kde_feas_below,
            aes(pvout, dens), fill = "grey88", alpha = 0.85, colour = NA) +
  geom_area(data = kde_feas_above,
            aes(pvout, dens), fill = "grey45", alpha = 0.85, colour = NA) +
  geom_line(data = kde_land_feas_n,
            aes(pvout, dens, colour = "Feasible land PVOUT"),
            linewidth = lwd_main) +
  geom_line(data = kde_pop_feas_n,
            aes(pvout, dens, colour = "Population-weighted (feasible land)"),
            linewidth = lwd_pop, linetype = "dashed") +
  geom_line(data = kde_c_upv_n,
            aes(pvout, dens, colour = "CPVPD UPV >=1 MW (capacity-weighted)"),
            linewidth = lwd_main) +
  geom_vline(xintercept = feas_med,           linetype = "dashed",
             colour = "grey20", linewidth = 0.5) +
  geom_vline(xintercept = cpvpd_unw_upv_med,  linetype = "dashed",
             colour = col_dep,  linewidth = 0.5) +
  geom_vline(xintercept = cpvpd_capw_mean_upv, linetype = "dotted",
             colour = col_dep,  linewidth = 0.5) +
  annotate("text", x = feas_med            + 12, y = ymax_c3 * 0.97,
           label = "Feasible\nmedian", size = 2.5, colour = "grey20",
           hjust = 0, vjust = 1, lineheight = 0.85) +
  annotate("text", x = cpvpd_unw_upv_med  + 12, y = ymax_c3 * 0.85,
           label = "Plant\nmedian",   size = 2.5, colour = col_dep,
           hjust = 0, vjust = 1, lineheight = 0.85) +
  annotate("text", x = cpvpd_capw_mean_upv + 12, y = ymax_c3 * 0.73,
           label = "Cap-wt\nmean",    size = 2.5, colour = col_dep,
           hjust = 0, vjust = 1, lineheight = 0.85) +
  scale_colour_manual(
    values = c("Feasible land PVOUT"                      = col_land,
               "Population-weighted (feasible land)"      = col_pop,
               "CPVPD UPV >=1 MW (capacity-weighted)"     = col_dep),
    name = NULL,
    guide = guide_legend(nrow = 1, override.aes = list(
      linetype  = c("solid", "dashed", "solid"),
      linewidth = c(lwd_main, lwd_pop, lwd_main)))) +
  xsc(xfrom_n, xto_n) + ysc +
  coord_cartesian(ylim = c(0, ymax_c3)) +
  theme_dens +
  labs(
    title   = "Figure C3  CPVPD UPV >=1 MW distribution vs feasible land PVOUT",
    caption = paste0(
      "Feasible land: elevation <=5000 m, slope <=50 deg. ",
      "Shading: lighter = below feasible median, darker = above.\n",
      "Dashed red: unweighted plant median. Dotted red: cap-weighted mean. ",
      "Navy dashed: population-weighted feasible land PVOUT (WorldPop 2020)."
    )
  )

ggsave(file.path(OUT, "fig_c3_cpvpd_upv_feasible.png"),
       fig_c3, width = 8, height = 4.5, dpi = 300)
cat("  Saved fig_c3_cpvpd_upv_feasible.png\n")

# =============================================================================
# Fig C4 — Direct comparison: GRW vs CPVPD on feasible land background
# =============================================================================
cat("Building Fig C4 (GRW vs CPVPD comparison)...\n")

ymax_c4 <- max(kde_land_feas_n$dens, kde_c_upv_n$dens,
               kde_g_capw_n$dens, kde_pop_feas_n$dens) * 1.08

fig_c4 <- ggplot() +
  geom_area(data = kde_feas_below,
            aes(pvout, dens), fill = "grey88", alpha = 0.85, colour = NA) +
  geom_area(data = kde_feas_above,
            aes(pvout, dens), fill = "grey45", alpha = 0.85, colour = NA) +
  geom_line(data = kde_land_feas_n,
            aes(pvout, dens, colour = "Feasible land PVOUT"),
            linewidth = lwd_main) +
  geom_line(data = kde_pop_feas_n,
            aes(pvout, dens, colour = "Population-weighted (feasible land)"),
            linewidth = lwd_pop, linetype = "dashed") +
  geom_line(data = kde_g_capw_n,
            aes(pvout, dens, colour = "GRW UPV >=1 MW (capacity-weighted)"),
            linewidth = lwd_main, linetype = "dashed") +
  geom_line(data = kde_c_upv_n,
            aes(pvout, dens, colour = "CPVPD UPV >=1 MW (capacity-weighted)"),
            linewidth = lwd_main) +
  # Reference verticals: feasible median (grey), GRW cap-wt mean (blue), CPVPD cap-wt mean (red)
  geom_vline(xintercept = feas_med,            linetype = "dashed",
             colour = "grey20", linewidth = 0.5) +
  geom_vline(xintercept = grw_capw_mean,        linetype = "dotted",
             colour = col_grw,  linewidth = 0.5) +
  geom_vline(xintercept = cpvpd_capw_mean_upv,  linetype = "dotted",
             colour = col_dep,  linewidth = 0.5) +
  annotate("text", x = feas_med             + 12, y = ymax_c4 * 0.97,
           label = "Feasible\nmedian", size = 2.5, colour = "grey20",
           hjust = 0, vjust = 1, lineheight = 0.85) +
  annotate("text", x = grw_capw_mean        + 12, y = ymax_c4 * 0.85,
           label = "GRW\ncap-wt mean",  size = 2.5, colour = col_grw,
           hjust = 0, vjust = 1, lineheight = 0.85) +
  annotate("text", x = cpvpd_capw_mean_upv  + 12, y = ymax_c4 * 0.73,
           label = "CPVPD\ncap-wt mean", size = 2.5, colour = col_dep,
           hjust = 0, vjust = 1, lineheight = 0.85) +
  scale_colour_manual(
    values = c("Feasible land PVOUT"                      = col_land,
               "Population-weighted (feasible land)"      = col_pop,
               "GRW UPV >=1 MW (capacity-weighted)"       = col_grw,
               "CPVPD UPV >=1 MW (capacity-weighted)"     = col_dep),
    name = NULL,
    guide = guide_legend(nrow = 2, override.aes = list(
      linetype  = c("solid", "dashed", "dashed", "solid"),
      linewidth = c(lwd_main, lwd_pop, lwd_main, lwd_main)))) +
  xsc(xfrom_n, xto_n) + ysc +
  coord_cartesian(ylim = c(0, ymax_c4)) +
  theme_dens +
  labs(
    title   = "Figure C4  GRW vs CPVPD-2024: UPV capacity-weighted distributions",
    caption = paste0(
      "Both datasets: UPV >=1 MW, capacity-weighted, against feasible land background.\n",
      "GRW: n=", length(dep_pvout), " plants, ", round(sum(dep_cap)/1e3, 1), " GW. ",
      "CPVPD: n=", nrow(cpvpd_upv), " plants, ",
      round(sum(cpvpd_upv$cap_mw)/1e3, 1), " GW. ",
      "Dotted verticals: cap-weighted means."
    )
  )

ggsave(file.path(OUT, "fig_c4_grw_vs_cpvpd_comparison.png"),
       fig_c4, width = 8, height = 4.5, dpi = 300)
cat("  Saved fig_c4_grw_vs_cpvpd_comparison.png\n")

cat("\nAll CPVPD validation figures saved to:", OUT, "\n")
cat(sprintf("\nKey scalars for comparison:\n"))
cat(sprintf("  All-land median PVOUT          : %.1f\n", all_land_med))
cat(sprintf("  Feasible-land median PVOUT     : %.1f\n", feas_med))
cat(sprintf("  GRW  UPV cap-wt mean           : %.1f\n", grw_capw_mean))
cat(sprintf("  CPVPD UPV cap-wt mean          : %.1f\n", cpvpd_capw_mean_upv))
cat(sprintf("  Difference (GRW - CPVPD)       : %.1f kWh/kWp/yr\n",
            grw_capw_mean - cpvpd_capw_mean_upv))
