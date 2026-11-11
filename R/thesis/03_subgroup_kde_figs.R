source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# 03_subgroup_kde_figs.R  —  Thesis Project 4
# Produces two sets of faceted KDE plots:
#
#   A.  By PROVINCE  — one panel per province, showing:
#         • feasible land PVOUT distribution (gray fill)
#         • population-weighted land PVOUT (navy dashed)
#         • cap-weighted deployed UPV PVOUT (red)
#       Provinces with fewer than MIN_PLANTS plants are dropped.
#
#   B.  By ELEVATION GROUP  — one panel per elevation band:
#         0–200 m | 200–400 m | 400–800 m | 800–1000 m | 1000–1200 m | >1200 m
#       Each panel compares feasible-land + population-weighted + cap-wt plants
#       restricted to cells/plants in that elevation band.
#
# Requires: 00_setup.R sourced (deploy_df, grid_df, upv_sf exist in env).
#           OR run standalone — data will be reloaded from disk.
# =============================================================================

library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(tidyr)

P3  <- FTS_PROCESSED
P4  <- FTS_ROOT
OUT <- file.path(P4, "figures")
dir.create(OUT, showWarnings = FALSE)

# ── Reload core data if not already in environment ───────────────────────────
if (!exists("deploy_df") || !exists("grid_df")) {
  cat("Sourcing 00_setup.R...\n")
  source(file.path(FTS_ROOT, "R", "thesis", "00_setup.R"))
}

CEA_CRS <- "+proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# ── Shared helpers ────────────────────────────────────────────────────────────
BW_ADJUST <- 2.0   # multiplier on nrd0 bandwidth — smooths spikes from narrow distributions

wkde <- function(x, weights = NULL, from, to, bw = "nrd0", n = 1024) {
  x <- x[is.finite(x)]
  if (length(x) < 5) return(NULL)
  if (!is.null(weights)) {
    weights <- weights[is.finite(x)]
    weights <- weights / sum(weights)
  }
  k <- density(x, weights = weights, bw = bw, adjust = BW_ADJUST,
               from = from, to = to, n = n)
  data.frame(pvout = k$x, dens = k$y)
}

col_land <- "grey70"
col_pop  <- "#1A3A5C"
col_dep  <- "#C0392B"
lwd_main <- 0.8
lwd_pop  <- 0.65

theme_facet <- theme_minimal(base_size = 9) +
  theme(
    panel.grid.minor    = element_blank(),
    panel.grid.major.y  = element_blank(),
    panel.border        = element_rect(colour = "grey85", fill = NA, linewidth = 0.3),
    axis.line.x   = element_line(colour = "grey40", linewidth = 0.25),
    axis.ticks.x  = element_line(colour = "grey40", linewidth = 0.25),
    axis.text     = element_text(size = 7),
    axis.title    = element_text(size = 8),
    strip.text    = element_text(size = 7.5, face = "bold"),
    legend.position   = "bottom",
    legend.key.width  = unit(0.9, "cm"),
    legend.key.height = unit(0.3, "cm"),
    legend.text       = element_text(size = 7.5),
    plot.title   = element_text(size = 10, face = "bold"),
    plot.caption = element_text(size = 7, colour = "grey50")
  )

# Shared x limits across all facets (use feasible-land range, clipped)
xfrom_global <- floor(quantile(grid_df$pvout_annual[grid_df$feasible &
                                 is.finite(grid_df$pvout_annual)], 0.001) / 50) * 50
xto_global   <- ceiling(quantile(grid_df$pvout_annual[grid_df$feasible &
                                 is.finite(grid_df$pvout_annual)], 0.999) / 50) * 50

# =============================================================================
# A.  PROVINCE-LEVEL KDE
# =============================================================================
cat("\n── Part A: Province KDE plots ──\n")

# Provinces that installed ≥1 GW post-shuangtan — the focus set for this figure
FOCUS_PROVS <- c(
  "Xinjiang", "Guizhou", "Hubei", "Anhui", "Shandong",
  "Guangdong", "Shaanxi", "Guangxi", "Jiangxi", "Jiangsu",
  "Hebei", "Hunan", "Zhejiang", "Henan", "Hainan"
)

MIN_PLANTS <- 30    # safety floor — drop any focus province that somehow has fewer

# ── A1. Chinese → English province name lookup ────────────────────────────────
# Covers all 34 provincial-level units. Uses substring matching so it works
# whether 地名 stores full formal names (e.g. "新疆维吾尔自治区") or short
# forms (e.g. "新疆").
zh_en_map <- c(
  "北京"   = "Beijing",       "天津"   = "Tianjin",
  "上海"   = "Shanghai",      "重庆"   = "Chongqing",
  "河北"   = "Hebei",         "山西"   = "Shanxi",
  "辽宁"   = "Liaoning",      "吉林"   = "Jilin",
  "黑龙江" = "Heilongjiang",  "江苏"   = "Jiangsu",
  "浙江"   = "Zhejiang",      "安徽"   = "Anhui",
  "福建"   = "Fujian",        "江西"   = "Jiangxi",
  "山东"   = "Shandong",      "河南"   = "Henan",
  "湖北"   = "Hubei",         "湖南"   = "Hunan",
  "广东"   = "Guangdong",     "海南"   = "Hainan",
  "四川"   = "Sichuan",       "贵州"   = "Guizhou",
  "云南"   = "Yunnan",        "陕西"   = "Shaanxi",
  "甘肃"   = "Gansu",         "青海"   = "Qinghai",
  "内蒙古" = "Inner Mongolia","广西"   = "Guangxi",
  "西藏"   = "Tibet",         "宁夏"   = "Ningxia",
  "新疆"   = "Xinjiang",      "香港"   = "Hong Kong",
  "澳门"   = "Macao",         "台湾"   = "Taiwan"
)

translate_prov <- function(x) {
  # Try exact key match first, then substring
  out <- zh_en_map[x]
  unmatched <- is.na(out)
  if (any(unmatched)) {
    for (i in which(unmatched)) {
      hit <- names(zh_en_map)[sapply(names(zh_en_map), grepl, x = x[i], fixed = TRUE)]
      out[i] <- if (length(hit) > 0) zh_en_map[hit[1]] else x[i]
    }
  }
  unname(out)
}

# Apply translation to deploy_df
deploy_df$prov_en <- translate_prov(deploy_df$prov)

# ── A2. Assign province to grid_df via spatial join if needed ─────────────────
if (!"province_en" %in% names(grid_df) ||
    all(is.na(grid_df$province_en)) ||
    (!"province" %in% names(grid_df) && !"province_zh" %in% names(grid_df))) {
  cat("  Province column missing from grid_df — assigning via spatial join...\n")
  sf_use_s2(FALSE)
  grid_sf     <- st_as_sf(grid_df, coords = c("x", "y"), crs = CEA_CRS, remove = FALSE)
  prov_sf_tmp <- st_transform(china_provinces, CEA_CRS)
  cp_cols     <- setdiff(names(prov_sf_tmp), "geometry")
  prov_name_col <- if ("\u5730\u540d" %in% cp_cols) "\u5730\u540d" else
                   if ("name"  %in% cp_cols) "name"  else cp_cols[1]
  names(prov_sf_tmp)[names(prov_sf_tmp) == prov_name_col] <- "province_zh"
  joined_grid <- st_join(grid_sf, prov_sf_tmp[, c("province_zh", "geometry")],
                         join = st_within, left = TRUE)
  unmatched <- which(is.na(joined_grid$province_zh))
  if (length(unmatched) > 0) {
    nearest_idx <- st_nearest_feature(grid_sf[unmatched, ], prov_sf_tmp)
    joined_grid$province_zh[unmatched] <- prov_sf_tmp$province_zh[nearest_idx]
  }
  grid_df$province_zh <- joined_grid$province_zh
  grid_df$province_en <- translate_prov(grid_df$province_zh)
  cat(sprintf("  Assigned province to %d grid cells.\n", nrow(grid_df)))
}

provs_with_plants <- deploy_df %>%
  filter(!is.na(prov_en), prov_en %in% FOCUS_PROVS) %>%
  count(prov_en) %>%
  filter(n >= MIN_PLANTS) %>%
  pull(prov_en)

cat(sprintf("  Focus provinces present with \u2265%d plants: %d\n", MIN_PLANTS, length(provs_with_plants)))

# Match grid province names to deploy_df province names.
# Both come from the same basemap (地名 column), so names should match.
grid_prov_col <- if ("province" %in% names(grid_df)) "province" else
                  names(grid_df)[grepl("prov", names(grid_df), ignore.case=TRUE)][1]

# ── A3. Build province facet data (keyed on English name) ────────────────────
build_province_kde <- function(prov_en_name) {
  # Land — match on English province name
  land_mask <- !is.na(grid_df$province_en) &
               grid_df$province_en == prov_en_name &
               grid_df$feasible &
               is.finite(grid_df$pvout_annual)
  land_pvout <- grid_df$pvout_annual[land_mask]
  if (length(land_pvout) < 5) return(NULL)

  # Population weights on feasible land
  pop_wt <- if ("pop_density" %in% names(grid_df)) {
    w <- grid_df$pop_density[land_mask]
    w[!is.finite(w) | w < 0] <- 1e-6
    w
  } else NULL

  # Plants — match on English province name
  plant_mask <- !is.na(deploy_df$prov_en) &
                deploy_df$prov_en == prov_en_name &
                is.finite(deploy_df$pvout) &
                is.finite(deploy_df$cap_mw) &
                deploy_df$cap_mw > 0
  plant_pvout <- deploy_df$pvout[plant_mask]
  plant_cap   <- deploy_df$cap_mw[plant_mask]
  if (length(plant_pvout) < 5) return(NULL)
  plant_w <- plant_cap / sum(plant_cap)

  kde_land  <- wkde(land_pvout, from = xfrom_global, to = xto_global)
  kde_plant <- wkde(plant_pvout, weights = plant_w,
                    from = xfrom_global, to = xto_global)
  kde_pop   <- if (!is.null(pop_wt))
    wkde(land_pvout, weights = pop_wt, from = xfrom_global, to = xto_global)
  else NULL

  if (is.null(kde_land) || is.null(kde_plant)) return(NULL)

  bind_rows(
    mutate(kde_land,  group = "Feasible land",         province = prov_en_name),
    mutate(kde_plant, group = "Deployed UPV (cap-wt)", province = prov_en_name),
    if (!is.null(kde_pop))
      mutate(kde_pop, group = "Population-weighted",   province = prov_en_name)
  )
}

cat("  Computing per-province KDEs...\n")
prov_kde_list <- lapply(provs_with_plants, build_province_kde)
prov_kde_df   <- bind_rows(Filter(Negate(is.null), prov_kde_list))

# Facet reference lines keyed on English names
prov_medians <- grid_df %>%
  filter(!is.na(province_en), province_en %in% provs_with_plants,
         feasible, is.finite(pvout_annual)) %>%
  group_by(province = province_en) %>%
  summarise(land_med = median(pvout_annual), .groups = "drop")

prov_plant_means <- deploy_df %>%
  filter(!is.na(prov_en), prov_en %in% provs_with_plants,
         is.finite(pvout), is.finite(cap_mw), cap_mw > 0) %>%
  group_by(province = prov_en) %>%
  summarise(cap_wt_mean = sum(pvout * cap_mw) / sum(cap_mw), .groups = "drop")

prov_refs <- left_join(prov_medians, prov_plant_means, by = "province")

# ── A4. Multi-page province plots ─────────────────────────────────────────────
# Layout: 3 columns × 2 rows = 6 provinces per page — larger, more readable panels.
# Sorted descending by cap-weighted mean PVOUT so the best-resourced provinces
# appear first; this creates a natural narrative across pages.

N_COLS <- 3L
N_ROWS <- 3L   # 9 per page → 2 pages for 15 provinces
PER_PAGE <- N_COLS * N_ROWS

# Province-level summary stats for sorting and strip labels
prov_summary <- deploy_df %>%
  filter(prov_en %in% provs_with_plants,
         is.finite(pvout), is.finite(cap_mw), cap_mw > 0) %>%
  group_by(prov_en) %>%
  summarise(
    capwt_mean = sum(pvout * cap_mw) / sum(cap_mw),
    total_gw   = sum(cap_mw) / 1e3,
    .groups = "drop"
  ) %>%
  arrange(desc(capwt_mean))   # highest cap-wt mean PVOUT first

# Reorder provs_with_plants by descending cap-wt mean
provs_with_plants <- prov_summary$prov_en

# Build a named vector of informative strip labels:
#   "Xinjiang\n1 523 kWh/kWp/yr | 87.4 GW"
prov_labels <- setNames(
  sprintf("%s\n%s kWh/kWp/yr | %.1f GW",
          prov_summary$prov_en,
          formatC(round(prov_summary$capwt_mean), format = "d", big.mark = ","),
          prov_summary$total_gw),
  prov_summary$prov_en
)

# Re-level factors with the new sorted order
prov_kde_df$province <- factor(prov_kde_df$province, levels = provs_with_plants)
prov_refs$province   <- factor(prov_refs$province,   levels = provs_with_plants)

n_pages <- ceiling(length(provs_with_plants) / PER_PAGE)

make_prov_plot <- function(page_provs, page_num, n_pages) {
  pd  <- filter(prov_kde_df, province %in% page_provs)
  ref <- filter(prov_refs,   province %in% page_provs)
  pd$province  <- factor(pd$province,  levels = page_provs)
  ref$province <- factor(ref$province, levels = page_provs)

  ggplot() +
    geom_area(data = filter(pd, group == "Feasible land"),
              aes(x = pvout, y = dens),
              fill = col_land, alpha = 0.55, colour = NA) +
    geom_line(data = filter(pd, group == "Feasible land"),
              aes(x = pvout, y = dens, colour = group),
              linewidth = lwd_main * 0.8) +
    geom_line(data = filter(pd, group == "Population-weighted"),
              aes(x = pvout, y = dens, colour = group),
              linewidth = lwd_pop, linetype = "dashed") +
    geom_line(data = filter(pd, group == "Deployed UPV (cap-wt)"),
              aes(x = pvout, y = dens, colour = group),
              linewidth = lwd_main) +
    geom_vline(data = ref, aes(xintercept = land_med),
               linetype = "dashed", colour = "grey35", linewidth = 0.4) +
    geom_vline(data = ref, aes(xintercept = cap_wt_mean),
               linetype = "dotted", colour = col_dep, linewidth = 0.4) +
    facet_wrap(~ province, ncol = N_COLS, scales = "free_y",
               labeller = as_labeller(prov_labels)) +
    scale_colour_manual(
      values = c("Feasible land"         = col_land,
                 "Population-weighted"   = col_pop,
                 "Deployed UPV (cap-wt)" = col_dep),
      name  = NULL,
      guide = guide_legend(nrow = 1, override.aes = list(
        linetype  = c("solid", "dashed", "solid"),
        linewidth = c(lwd_main * 0.8, lwd_pop, lwd_main)))) +
    scale_x_continuous(name   = "PVOUT (kWh/kWp/yr)",
                       limits = c(xfrom_global, xto_global),
                       breaks = pretty(c(xfrom_global, xto_global), n = 4),
                       expand = c(0, 0)) +
    scale_y_continuous(name   = "Density",
                       expand = expansion(mult = c(0, 0.08))) +
    theme_facet +
    theme(strip.text = element_text(size = 8, face = "bold", lineheight = 1.15)) +
    labs(
      title   = sprintf(
        "PVOUT distribution by province: feasible land, population-weighted, and deployed UPV  [page %d of %d]",
        page_num, n_pages),
      caption = paste0(
        "Provinces sorted by cap-weighted mean PVOUT (highest first). ",
        "Strip subtitle: cap-wt mean PVOUT | total deployed capacity.\n",
        "Gray fill: feasible land. Red: UPV capacity-weighted. Navy dashed: population-weighted. ",
        "Dashed vline: feasible-land median. Dotted vline: cap-wt plant mean. ",
        "X-axis fixed across all panels; y-axis free per province. Provinces with \u22651 GW installed post-shuangtan.")
    )
}

# ── A4b. High-performing provinces — dedicated 2×2 panel ─────────────────────
HIGH_PERF <- c("Qinghai", "Inner Mongolia", "Ningxia", "Tibet")

cat("  Building high-performer panel...\n")
hp_kde_list <- lapply(HIGH_PERF, build_province_kde)
hp_kde_df   <- bind_rows(Filter(Negate(is.null), hp_kde_list))

hp_summary <- deploy_df %>%
  filter(prov_en %in% HIGH_PERF,
         is.finite(pvout), is.finite(cap_mw), cap_mw > 0) %>%
  group_by(prov_en) %>%
  summarise(
    capwt_mean = sum(pvout * cap_mw) / sum(cap_mw),
    total_gw   = sum(cap_mw) / 1e3,
    .groups = "drop"
  ) %>%
  arrange(desc(capwt_mean))

hp_labels <- setNames(
  sprintf("%s\n%s kWh/kWp/yr | %.1f GW",
          hp_summary$prov_en,
          formatC(round(hp_summary$capwt_mean), format = "d", big.mark = ","),
          hp_summary$total_gw),
  hp_summary$prov_en
)

hp_order <- hp_summary$prov_en
hp_kde_df$province <- factor(hp_kde_df$province, levels = hp_order)

hp_refs <- prov_refs %>% filter(province %in% HIGH_PERF) %>%
  mutate(province = factor(province, levels = hp_order))

fig_hp <- ggplot() +
  geom_area(data = filter(hp_kde_df, group == "Feasible land"),
            aes(x = pvout, y = dens),
            fill = col_land, alpha = 0.55, colour = NA) +
  geom_line(data = filter(hp_kde_df, group == "Feasible land"),
            aes(x = pvout, y = dens, colour = group),
            linewidth = lwd_main * 0.8) +
  geom_line(data = filter(hp_kde_df, group == "Population-weighted"),
            aes(x = pvout, y = dens, colour = group),
            linewidth = lwd_pop, linetype = "dashed") +
  geom_line(data = filter(hp_kde_df, group == "Deployed UPV (cap-wt)"),
            aes(x = pvout, y = dens, colour = group),
            linewidth = lwd_main) +
  geom_vline(data = hp_refs, aes(xintercept = land_med),
             linetype = "dashed", colour = "grey35", linewidth = 0.4) +
  geom_vline(data = hp_refs, aes(xintercept = cap_wt_mean),
             linetype = "dotted", colour = col_dep, linewidth = 0.4) +
  facet_wrap(~ province, ncol = 2, scales = "free_y",
             labeller = as_labeller(hp_labels)) +
  scale_colour_manual(
    values = c("Feasible land"         = col_land,
               "Population-weighted"   = col_pop,
               "Deployed UPV (cap-wt)" = col_dep),
    name  = NULL,
    guide = guide_legend(nrow = 1, override.aes = list(
      linetype  = c("solid", "dashed", "solid"),
      linewidth = c(lwd_main * 0.8, lwd_pop, lwd_main)))) +
  scale_x_continuous(name   = "PVOUT (kWh/kWp/yr)",
                     limits = c(xfrom_global, xto_global),
                     breaks = pretty(c(xfrom_global, xto_global), n = 4),
                     expand = c(0, 0)) +
  scale_y_continuous(name   = "Density",
                     expand = expansion(mult = c(0, 0.08))) +
  theme_facet +
  theme(strip.text = element_text(size = 9, face = "bold", lineheight = 1.15)) +
  labs(
    title   = "PVOUT distribution — high-resource provinces: feasible land, population-weighted, and deployed UPV",
    caption = paste0(
      "Provinces sorted by cap-weighted mean PVOUT (highest first). ",
      "Strip subtitle: cap-wt mean PVOUT | total deployed capacity.\n",
      "Gray fill: feasible land. Red: UPV capacity-weighted. Navy dashed: population-weighted. ",
      "Dashed vline: feasible-land median. Dotted vline: cap-wt plant mean. ",
      "X-axis fixed across all panels; y-axis free per province.")
  )

ggsave(file.path(OUT, "fig_A_province_kde_highperf.png"),
       fig_hp, width = 10, height = 8, dpi = 300)
cat("  Saved fig_A_province_kde_highperf.png\n")

# ── A4c. Multi-page: focus provinces ─────────────────────────────────────────
cat(sprintf("  Saving %d pages of province KDE plots...\n", n_pages))
for (pg in seq_len(n_pages)) {
  idx        <- ((pg - 1) * PER_PAGE + 1) : min(pg * PER_PAGE, length(provs_with_plants))
  page_provs <- provs_with_plants[idx]
  n_rows_pg  <- ceiling(length(page_provs) / N_COLS)
  fig_h      <- 2.0 + n_rows_pg * 3.5   # taller panels for legibility

  p <- make_prov_plot(page_provs, pg, n_pages)
  fname <- sprintf("fig_A_province_kde_p%d.png", pg)
  ggsave(file.path(OUT, fname), p,
         width = N_COLS * 4.0, height = fig_h, dpi = 300)
  cat(sprintf("    Saved %s  (%d provinces)\n", fname, length(page_provs)))
}

# =============================================================================
# B.  ELEVATION-GROUP KDE
# =============================================================================
cat("\n── Part B: Elevation group KDE plots ──\n")

# ── B1. Add elevation to grid_df and deploy_df ───────────────────────────────

# Grid already has elevation_m
if (!"elevation_m" %in% names(grid_df)) {
  cat("  Extracting elevation for grid cells...\n")
  elev_r <- rast(file.path(P3, "elevation_10km.tif"))
  grid_pts <- st_as_sf(grid_df, coords = c("x","y"), crs = st_crs(elev_r))
  grid_df$elevation_m <- terra::extract(elev_r, vect(grid_pts))[[2]]
}

# Plants: extract elevation from 10km raster at centroid locations
if (!"elevation_m" %in% names(deploy_df)) {
  cat("  Extracting elevation for plant centroids...\n")
  sf_use_s2(FALSE)
  if (!exists("upv_sf")) upv_sf <- readRDS(file.path(P3, "UPV_core_GRW_2024Q2_ge1MW.rds"))
  upv_sf <- st_make_valid(upv_sf)
  elev_r   <- rast(file.path(P3, "elevation_10km.tif"))
  cents    <- st_centroid(upv_sf)
  elev_ext <- terra::extract(elev_r, vect(cents))
  elev_all <- as.numeric(elev_ext[[2]])   # length = nrow(upv_sf) = 16714

  # deploy_df was built in 00_setup.R with:
  #   filter(is.finite(pvout), is.finite(cap_mw), cap_mw > 0, !is.na(prov))
  # Reconstruct the same logical mask so elev_all aligns with deploy_df rows.
  cap_all <- as.numeric(upv_sf$capacity_mw_used)
  if ("pvout" %in% names(upv_sf) && "prov" %in% names(upv_sf)) {
    keep <- is.finite(as.numeric(upv_sf$pvout)) &
            is.finite(cap_all) & cap_all > 0 &
            !is.na(upv_sf$prov)
  } else {
    # Fallback if upv_sf was reloaded without pvout/prov columns
    keep <- is.finite(cap_all) & cap_all > 0
  }

  elev_filtered <- elev_all[keep]
  if (length(elev_filtered) != nrow(deploy_df)) {
    stop(sprintf(
      "Row mismatch after filter: elevation has %d rows, deploy_df has %d. ",
      length(elev_filtered), nrow(deploy_df),
      "Check that upv_sf has pvout and prov columns (source 00_setup.R first)."))
  }
  deploy_df$elevation_m <- elev_filtered
  cat(sprintf("  Elevation assigned to %d plants.\n", nrow(deploy_df)))
}

# ── B2. Define elevation bands ───────────────────────────────────────────────
elev_breaks <- c(0, 200, 400, 800, 1000, 1200, Inf)
elev_labels <- c("0–200 m", "200–400 m", "400–800 m",
                 "800–1000 m", "1000–1200 m", ">1200 m")

grid_df <- grid_df %>%
  mutate(elev_group = cut(elevation_m, breaks = elev_breaks,
                          labels = elev_labels, include.lowest = TRUE,
                          right = FALSE))

deploy_df <- deploy_df %>%
  mutate(elev_group = cut(elevation_m, breaks = elev_breaks,
                          labels = elev_labels, include.lowest = TRUE,
                          right = FALSE))

# ── B3. Build elevation-group facet data ─────────────────────────────────────
build_elev_kde <- function(eg) {
  # Land
  land_mask <- !is.na(grid_df$elev_group) &
               grid_df$elev_group == eg &
               grid_df$feasible &
               is.finite(grid_df$pvout_annual)
  land_pvout <- grid_df$pvout_annual[land_mask]
  if (length(land_pvout) < 5) return(NULL)

  # Population weights
  pop_col <- if ("pop_density" %in% names(grid_df)) "pop_density" else NULL
  pop_wt  <- if (!is.null(pop_col)) {
    w <- grid_df[[pop_col]][land_mask]
    w[!is.finite(w) | w < 0] <- 1e-6
    w
  } else NULL

  # Plants
  plant_mask <- !is.na(deploy_df$elev_group) &
                deploy_df$elev_group == eg &
                is.finite(deploy_df$pvout) &
                is.finite(deploy_df$cap_mw) &
                deploy_df$cap_mw > 0
  plant_pvout <- deploy_df$pvout[plant_mask]
  plant_cap   <- deploy_df$cap_mw[plant_mask]
  if (length(plant_pvout) < 5) return(NULL)
  plant_w <- plant_cap / sum(plant_cap)

  kde_land  <- wkde(land_pvout, from = xfrom_global, to = xto_global)
  kde_plant <- wkde(plant_pvout, weights = plant_w,
                    from = xfrom_global, to = xto_global)
  kde_pop   <- if (!is.null(pop_wt))
    wkde(land_pvout, weights = pop_wt, from = xfrom_global, to = xto_global)
  else NULL

  if (is.null(kde_land) || is.null(kde_plant)) return(NULL)

  bind_rows(
    mutate(kde_land,  group = "Feasible land",          elev_group = eg),
    mutate(kde_plant, group = "Deployed UPV (cap-wt)",  elev_group = eg),
    if (!is.null(kde_pop))
      mutate(kde_pop, group = "Population-weighted", elev_group = eg)
  )
}

cat("  Computing per-elevation-band KDEs...\n")
elev_kde_list <- lapply(elev_labels, build_elev_kde)
elev_kde_df   <- bind_rows(Filter(Negate(is.null), elev_kde_list))

# Keep factor order for faceting
elev_kde_df$elev_group <- factor(elev_kde_df$elev_group, levels = elev_labels)

# Per-band reference lines
elev_refs <- grid_df %>%
  filter(!is.na(elev_group), feasible, is.finite(pvout_annual)) %>%
  group_by(elev_group) %>%
  summarise(land_med = median(pvout_annual), .groups = "drop") %>%
  left_join(
    deploy_df %>%
      filter(!is.na(elev_group), is.finite(pvout), is.finite(cap_mw), cap_mw > 0) %>%
      group_by(elev_group) %>%
      summarise(cap_wt_mean = sum(pvout * cap_mw) / sum(cap_mw), .groups = "drop"),
    by = "elev_group"
  )

elev_refs$elev_group <- factor(elev_refs$elev_group, levels = elev_labels)

# Plant counts per band (for caption / strip annotation)
plant_counts <- deploy_df %>%
  filter(!is.na(elev_group), is.finite(pvout), is.finite(cap_mw), cap_mw > 0) %>%
  count(elev_group) %>%
  mutate(label = paste0("n = ", scales::comma(n)))

# ── B4. Plot ──────────────────────────────────────────────────────────────────
fig_elev <- ggplot() +
  geom_area(data  = filter(elev_kde_df, group == "Feasible land"),
            aes(x = pvout, y = dens),
            fill = col_land, alpha = 0.55, colour = NA) +
  geom_line(data  = filter(elev_kde_df, group == "Feasible land"),
            aes(x = pvout, y = dens, colour = group),
            linewidth = lwd_main * 0.8) +
  geom_line(data  = filter(elev_kde_df, group == "Population-weighted"),
            aes(x = pvout, y = dens, colour = group),
            linewidth = lwd_pop, linetype = "dashed") +
  geom_line(data  = filter(elev_kde_df, group == "Deployed UPV (cap-wt)"),
            aes(x = pvout, y = dens, colour = group),
            linewidth = lwd_main) +
  geom_vline(data = elev_refs,
             aes(xintercept = land_med),
             linetype = "dashed", colour = "grey35", linewidth = 0.4) +
  geom_vline(data = elev_refs,
             aes(xintercept = cap_wt_mean),
             linetype = "dotted", colour = col_dep, linewidth = 0.4) +
  facet_wrap(~ elev_group, ncol = 3) +   # fixed axes (scales = "fixed" is default)
  scale_colour_manual(
    values = c("Feasible land"           = col_land,
               "Population-weighted"     = col_pop,
               "Deployed UPV (cap-wt)"   = col_dep),
    name  = NULL,
    guide = guide_legend(nrow = 1, override.aes = list(
      linetype  = c("solid", "dashed", "solid"),
      linewidth = c(lwd_main * 0.8, lwd_pop, lwd_main)))) +
  scale_x_continuous(name   = "PVOUT (kWh/kWp/yr)",
                     limits = c(xfrom_global, xto_global),
                     breaks = pretty(c(xfrom_global, xto_global), n = 5),
                     expand = c(0, 0)) +
  scale_y_continuous(name   = "Density",
                     expand = expansion(mult = c(0, 0.08))) +
  theme_facet +
  labs(
    title   = "PVOUT distribution by elevation band: feasible land, population-weighted, and deployed UPV",
    caption = paste0(
      "Gray fill: feasible land within elevation band. Red: UPV plants (capacity-weighted). ",
      "Navy dashed: population-weighted land PVOUT.\n",
      "Dashed vline: band feasible-land median. Dotted vline: band cap-wt plant mean. ",
      "Elevation from SRTM (10 km aggregated). Population: WorldPop 2020.")
  )

ggsave(file.path(OUT, "fig_B_elevation_group_kde.png"),
       fig_elev, width = 11, height = 7, dpi = 300)
cat("  Saved fig_B_elevation_group_kde.png\n")

# =============================================================================
# Summary
# =============================================================================
cat("\nDone. Output files:\n")
cat("  figures/fig_A_province_kde.png\n")
cat("  figures/fig_B_elevation_group_kde.png\n")

cat("\nPlant counts by elevation band:\n")
print(
  deploy_df %>%
    filter(!is.na(elev_group), is.finite(pvout), is.finite(cap_mw), cap_mw > 0) %>%
    group_by(elev_group) %>%
    summarise(n_plants = n(),
              total_gw = round(sum(cap_mw) / 1e3, 1),
              med_pvout = round(median(pvout), 0),
              .groups = "drop")
)
