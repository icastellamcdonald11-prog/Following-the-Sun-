source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# make_figures.R: Figs 1-5 for the NCC-format manuscript (Phase 7)
#
# Run from the project root:
#   Rscript paper/figures/make_figures.R          # all figures
#   Rscript paper/figures/make_figures.R 1 3      # selected figures
#
# Inputs (no analysis is recomputed except where stated):
#   fig_cache/plants_quad.rds   plants: centroid PVOUT, capacity, province,
#                               2020 substation distance, penalty, quadrant
#                               (built by thesis_figures.R)
#   fig_cache/grid_access.rds   10 km land grid: PVOUT, feasibility, province,
#                               2020 grid access (built by thesis_figures.R)
#   Project 3 processed data    basemap bundle, GRW polygons, 10 km elevation
#   pvout_1km_cea.tif, topo_mask_1km_cea.tif  (maps in Fig. 4c,d)
#   Transmission Analysis refined/grid_nodes_2025.csv  (Fig. 5)
#
# Style (NCC): 183 mm width, Helvetica 7-8 pt, bold lowercase panel letters,
# no in-panel titles or captions (legends carry them), one colour-blind-safe
# palette, quadrant colours identical in Figs 2 and 5. Output: vector PDF.
# =============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(scales); library(patchwork)
  library(ggrepel); library(sf); library(terra)
})
sf_use_s2(FALSE)
invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))  # Latin-1 glyphs (Ü) in the pdf device

P4  <- FTS_ROOT
P3  <- FTS_PROCESSED
OUT <- file.path(FTS_ROOT, "outputs", "figures")
setwd(P4)

# ----------------------------------------------------------------- style ----
PT_TEXT  <- 7      # axis text, legend text, annotations
PT_TITLE <- 7.5    # axis titles
PT_TAG   <- 8      # panel letters
LW_THIN  <- 0.25; LW_MED <- 0.4; LW_THICK <- 0.7
W_DOUBLE <- 183

# Quadrant colours: Okabe-Ito blue/orange pair, CVD-tested (thesis_figures.R).
QUAD_LEVELS <- c("Premium", "Stranded", "Grid-convenient", "Dominated")
QUAD_COLS <- c(Premium = "#0060A0", Stranded = "#7FC7EF",
               `Grid-convenient` = "#F2C166", Dominated = "#9A4200")
# Okabe-Ito accents for lines (no red/green, per supervisor)
OI_BLUE <- "#0072B2"; OI_ORANGE <- "#E69F00"; OI_SKY <- "#56B4E9"
OI_BLACK <- "#000000"; OI_VERMIL <- "#D55E00"
GREY_DARK <- "#3B3B3B"; GREY_MID <- "#7A7A7A"; GREY_LITE <- "#C8C8C8"

PROV_EN <- c(Neimenggu = "Inner Mongolia", Xizang = "Tibet", Aomen = "Macao",
             HongKong = "Hong Kong")
prov_en <- function(x) ifelse(x %in% names(PROV_EN), PROV_EN[x], x)

theme_ncc <- function() {
  theme_classic(base_size = PT_TEXT, base_family = "Helvetica") +
    theme(
      text = element_text(colour = OI_BLACK, size = PT_TEXT),
      axis.text = element_text(colour = OI_BLACK, size = PT_TEXT),
      axis.title = element_text(size = PT_TITLE),
      axis.line = element_line(linewidth = LW_THIN),
      axis.ticks = element_line(linewidth = LW_THIN),
      axis.ticks.length = unit(1.5, "pt"),
      legend.text = element_text(size = PT_TEXT),
      legend.title = element_text(size = PT_TEXT),
      legend.key.size = unit(8, "pt"),
      legend.background = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      plot.tag = element_text(size = PT_TAG, face = "bold"),
      plot.tag.position = c(0, 1),
      plot.margin = margin(6, 6, 4, 6)
    )
}
theme_map <- function() {
  theme_ncc() + theme(axis.line = element_blank(), axis.text = element_blank(),
                      axis.ticks = element_blank(), axis.title = element_blank())
}
lab_pvout <- expression("PVOUT (kWh kWp"^-1*" yr"^-1*")")

save_pdf <- function(p, n, h_mm) {
  f <- file.path(OUT, sprintf("fig%d.pdf", n))
  ggsave(f, p, device = grDevices::pdf, width = W_DOUBLE, height = h_mm,
         units = "mm", useDingbats = FALSE)
  ggsave(sub("\\.pdf$", "_preview.png", f), p, width = W_DOUBLE, height = h_mm,
         units = "mm", dpi = 200)
  cat(sprintf("  wrote %s (%d x %d mm, %.0f KB)\n", basename(f), W_DOUBLE,
              as.integer(h_mm), file.size(f) / 1024))
}

# ------------------------------------------------------------------ data ----
plants  <- readRDS("fig_cache/plants_quad.rds")
grid_df <- readRDS("fig_cache/grid_access.rds")
plants$quadrant <- factor(plants$quadrant, levels = QUAD_LEVELS)
PVOUT_MED <- median(grid_df$pvout_annual[grid_df$feasible & is.finite(grid_df$pvout_annual)])
TX_MED    <- median(plants$tx_penalty)
qn <- function(x) as.integer(substr(x, 1, 4)) + (as.integer(substr(x, 6, 6)) - 1) / 4
plants$t <- qn(plants$year_quarter)
cat(sprintf("Data: %d plants (%.1f GW); PVOUT median %.1f; penalty median %.2f\n",
            nrow(plants), sum(plants$cap_mw) / 1e3, PVOUT_MED, TX_MED))
fmt <- function(x, d = 0) formatC(x, format = "f", digits = d, big.mark = ",")

# Kernel densities on a common grid (bandwidth as in thesis_figures.R)
BW <- 45
dens <- function(x, w = NULL, rng) {
  if (!is.null(w)) w <- w / sum(w)
  d <- density(x, weights = w, bw = BW, from = rng[1], to = rng[2], n = 512)
  data.frame(x = d$x, y = d$y)
}
DENS_COLS <- c("Feasible land" = GREY_LITE, "Grid-accessible land" = OI_ORANGE,
               "Deployed plants" = OI_BLUE, "Deployed capacity" = OI_BLUE)

density_panel <- function(land, land_acc, dep_x, dep_w, dep_label, rng, ylim, tag,
                          median_lab = TRUE) {
  df_l <- dens(land, rng = rng); df_a <- dens(land_acc, rng = rng)
  df_d <- dens(dep_x, dep_w, rng)
  p <- ggplot() +
    geom_area(data = df_l, aes(x, y, fill = "Feasible land"), alpha = 0.8) +
    geom_line(data = df_a, aes(x, y, colour = "Grid-accessible land"),
              linewidth = LW_THICK, linetype = "22") +
    geom_line(data = df_d, aes(x, y, colour = dep_label), linewidth = LW_THICK) +
    geom_vline(xintercept = PVOUT_MED, linetype = "dashed", linewidth = LW_MED,
               colour = GREY_DARK) +
    scale_fill_manual(values = DENS_COLS, name = NULL,
                      guide = guide_legend(order = 1)) +
    scale_colour_manual(values = DENS_COLS, name = NULL,
                        breaks = c("Grid-accessible land", dep_label),
                        guide = guide_legend(order = 2)) +
    scale_x_continuous(breaks = seq(800, 2200, 200), labels = label_comma()) +
    scale_y_continuous(expand = expansion(c(0, 0.04)), breaks = NULL) +
    coord_cartesian(xlim = rng, ylim = ylim) +
    labs(x = lab_pvout, y = "Density", tag = tag) +
    theme_ncc() +
    theme(legend.position = "inside", legend.position.inside = c(0.01, 0.99),
          legend.justification = c(0, 1), legend.spacing.y = unit(0, "pt"))
  if (median_lab)
    p <- p + annotate("text", x = PVOUT_MED, y = ylim[2], hjust = -0.05, vjust = 1,
                      size = PT_TEXT / .pt, colour = GREY_DARK,
                      label = sprintf("National median\n%s", fmt(PVOUT_MED)))
  p
}

# ============================================================== Figure 1 ====
fig1 <- function() {
  feas <- grid_df |> filter(feasible, is.finite(pvout_annual))
  land <- feas$pvout_annual; acc <- feas$pvout_annual[feas$grid_access]
  rng <- c(700, 2300)
  ymax <- max(dens(land, rng = rng)$y, dens(acc, rng = rng)$y,
              dens(plants$pvout, rng = rng)$y,
              dens(plants$pvout, plants$cap_mw, rng)$y) * 1.12
  a <- density_panel(land, acc, plants$pvout, NULL, "Deployed plants", rng,
                     c(0, ymax), "a")
  b <- density_panel(land, acc, plants$pvout, plants$cap_mw, "Deployed capacity",
                     rng, c(0, ymax), "b")
  a | b
}

# ============================================================== Figure 2 ====
fig2 <- function(dat = plants) {
  tot <- sum(dat$cap_mw)
  shr <- tapply(dat$cap_mw, dat$quadrant, sum) / tot * 100
  xr <- c(700, 2250); yr <- c(0, quantile(dat$tx_penalty, 0.995))
  lab <- data.frame(
    quadrant = factor(QUAD_LEVELS, levels = QUAD_LEVELS),
    x = c(xr[2], xr[2], xr[1], xr[1]), y = c(0.6, yr[2], 0.6, yr[2]),
    hj = c(1, 1, 0, 0), vj = c(0, 1, 0, 1),
    txt = sprintf("%s\n%.1f%% of capacity", QUAD_LEVELS, shr[QUAD_LEVELS]))
  ggplot(dat, aes(pvout, tx_penalty)) +
    geom_density_2d(colour = GREY_MID, linewidth = LW_THIN, alpha = 0.6, bins = 6) +
    geom_point(aes(colour = quadrant, size = cap_mw), alpha = 0.55, stroke = 0) +
    geom_hline(yintercept = TX_MED, linetype = "dashed", colour = GREY_DARK,
               linewidth = LW_MED) +
    geom_vline(xintercept = PVOUT_MED, linetype = "dashed", colour = GREY_DARK,
               linewidth = LW_MED) +
    annotate("label", x = PVOUT_MED, y = yr[2] * 0.62, label = fmt(PVOUT_MED),
             size = PT_TEXT / .pt, linewidth = 0, fill = "white",
             label.padding = unit(1, "pt")) +
    annotate("label", x = 2200, y = TX_MED, label = fmt(TX_MED, 1),
             size = PT_TEXT / .pt, linewidth = 0, fill = "white",
             label.padding = unit(1, "pt")) +
    geom_label(data = lab, aes(x, y, label = txt, hjust = hj, vjust = vj,
                               colour = quadrant),
               size = PT_TEXT / .pt, lineheight = 0.9, fontface = "bold",
               linewidth = 0, fill = alpha("white", 0.85),
               label.padding = unit(1.5, "pt"), show.legend = FALSE) +
    scale_colour_manual(values = QUAD_COLS, guide = "none") +
    scale_size_area(max_size = 5, breaks = c(100, 1000, 5000),
                    labels = c("100", "1,000", "5,000"), name = "Plant capacity (MW)") +
    scale_y_sqrt(breaks = c(0, 5, 25, 100, 250, 500), expand = expansion(c(0, 0.02))) +
    scale_x_continuous(breaks = seq(800, 2200, 200), labels = label_comma()) +
    coord_cartesian(xlim = xr, ylim = yr) +
    guides(size = guide_legend(override.aes = list(colour = GREY_MID, alpha = 0.9,
                                                   size = c(1.4, 3, 5)),
                               nrow = 1, title.position = "left")) +
    labs(x = lab_pvout,
         y = expression("Transmission cost penalty (kWh kWp"^-1*" yr"^-1*")")) +
    theme_ncc() + theme(legend.position = "bottom")
}

# ============================================================== Figure 3 ====
# Cumulative capacity-weighted PVOUT of the installed stock, by quarter.
stock_series <- function(dat) {
  qs <- sort(unique(dat$t))
  data.frame(t = qs, pv = sapply(qs, function(q) {
    s <- dat[dat$t <= q, ]; weighted.mean(s$pvout, s$cap_mw) }))
}
# Equation 12 benchmark: capacity-weighted mean PVOUT of capacity installed
# before Q2 2020 (1,510.6 kWh kWp-1 yr-1; DECISIONS D9).
DET_BENCH <- with(plants[plants$t < 2020.25, ], weighted.mean(pvout, cap_mw))

fig3 <- function() {
  s_all <- stock_series(plants)
  xb <- scale_x_continuous(breaks = 2018:2024, limits = c(2017.7, 2024.5))
  a <- ggplot(s_all, aes(t, pv)) +
    geom_vline(xintercept = 2020.25, linetype = "dotted", colour = GREY_MID,
               linewidth = LW_MED) +
    annotate("text", x = 2020.25, y = Inf, vjust = 1.2, hjust = -0.05,
             label = "Q2 2020", size = PT_TEXT / .pt, colour = GREY_DARK) +
    geom_line(colour = OI_BLACK, linewidth = LW_MED) +
    geom_point(colour = OI_BLACK, size = 0.9) +
    xb + scale_y_continuous(labels = label_comma()) +
    labs(x = "Quarter", y = expression("Capacity-weighted PVOUT (kWh kWp"^-1*" yr"^-1*")"),
         tag = "a") + theme_ncc()

  # post-2020 = commissioned in Q2 2020 or later (Methods; analysis/deterioration_shares.R)
  post <- plants |> filter(t >= 2020.25) |> group_by(prov) |>
    summarise(cap_gw = sum(cap_mw) / 1e3, wpv = weighted.mean(pvout, cap_mw),
              .groups = "drop") |>
    mutate(gap = wpv - DET_BENCH, contrib = cap_gw * gap,
           dir = ifelse(gap >= 0, "Above benchmark", "Below benchmark"),
           name = prov_en(prov)) |>
    filter(cap_gw > 0.15)
  lab_b <- post |> arrange(desc(abs(contrib))) |> slice_head(n = 12)
  b <- ggplot(post, aes(cap_gw, gap)) +
    geom_hline(yintercept = 0, colour = GREY_DARK, linewidth = LW_MED) +
    geom_point(aes(size = abs(contrib), fill = dir), shape = 21, colour = "white",
               stroke = 0.25, alpha = 0.9) +
    geom_text_repel(data = lab_b, aes(label = name), size = PT_TEXT / .pt, seed = 2,
                    min.segment.length = 0.1, segment.size = LW_THIN,
                    segment.colour = GREY_MID, box.padding = 0.3, max.overlaps = 30) +
    scale_fill_manual(values = c("Above benchmark" = QUAD_COLS[["Premium"]],
                                 "Below benchmark" = QUAD_COLS[["Dominated"]]), name = NULL) +
    scale_size_area(max_size = 8, guide = "none") +
    labs(x = "Capacity added from Q2 2020 (GW)",
         y = expression(atop("PVOUT of additions minus benchmark",
                             "(kWh kWp"^-1*" yr"^-1*")")), tag = "b") +
    theme_ncc() +
    theme(legend.position = "inside", legend.position.inside = c(0.99, 0.99),
          legend.justification = c(1, 1))

  s <- bind_rows(
    mutate(s_all, series = "All provinces"),
    mutate(stock_series(filter(plants, prov != "Qinghai")), series = "Excluding Qinghai"),
    mutate(stock_series(filter(plants, prov != "Xinjiang")), series = "Excluding Xinjiang"))
  s$series <- factor(s$series, levels = c("All provinces", "Excluding Qinghai",
                                          "Excluding Xinjiang"))
  c_ <- ggplot(s, aes(t, pv, colour = series, linetype = series, shape = series)) +
    geom_hline(yintercept = PVOUT_MED, colour = GREY_MID, linewidth = LW_MED) +
    annotate("text", x = 2017.75, y = PVOUT_MED, vjust = -0.4, hjust = 0,
             label = sprintf("National feasible-land median (%s)", fmt(PVOUT_MED)),
             size = PT_TEXT / .pt, colour = GREY_DARK) +
    geom_line(linewidth = LW_MED) + geom_point(size = 1.1, fill = "white") +
    scale_colour_manual(values = c(OI_BLACK, OI_BLUE, OI_ORANGE), name = NULL) +
    scale_linetype_manual(values = c("solid", "longdash", "dotted"), name = NULL) +
    scale_shape_manual(values = c(16, 24, 22), name = NULL) +
    xb + scale_y_continuous(labels = label_comma()) +
    labs(x = "Quarter", y = expression("Capacity-weighted PVOUT (kWh kWp"^-1*" yr"^-1*")"),
         tag = "c") +
    theme_ncc() +
    theme(legend.position = "inside", legend.position.inside = c(0.99, 0.99),
          legend.justification = c(1, 1), legend.key.width = unit(22, "pt"))

  attr_out <- list(all = s_all, s = s)
  p <- (a | b) / c_ + plot_layout(heights = c(1, 1))
  attr(p, "series") <- attr_out
  p
}

# ============================================================== shared geo ====
basemap <- readRDS(file.path(P3, "china_basemap_bundle.rds"))
provs   <- basemap$china_provinces[, c("name", "geometry")]
prov_poly <- function(nm) st_union(provs[provs$name == nm, ])

# ============================================================== Figure 4 ====
# a,b: provincial densities on shared axes; c,d: maps on a shared PVOUT scale.
province_density <- function(nm, rng, ylim, tag) {
  g <- grid_df |> filter(feasible, is.finite(pvout_annual), prov == nm)
  p <- plants |> filter(prov == nm)
  density_panel(g$pvout_annual, g$pvout_annual[g$grid_access], p$pvout, p$cap_mw,
                "Deployed capacity", rng, ylim, tag, median_lab = TRUE) +
    annotate("text", x = rng[2], y = ylim[2], hjust = 1, vjust = 1,
             label = prov_en(nm), size = PT_TITLE / .pt, fontface = "bold")
}

CITIES <- data.frame(name = c("Xining", "Ürümqi"), prov = c("Qinghai", "Xinjiang"),
                     lon = c(101.778, 87.617), lat = c(36.617, 43.825))
MAP_CRS <- c(Qinghai  = "+proj=laea +lat_0=35.8 +lon_0=96 +datum=WGS84 +units=m",
             Xinjiang = "+proj=laea +lat_0=41.5 +lon_0=85 +datum=WGS84 +units=m")
PV_LIMS <- c(1000, 2200)

map_panel <- function(nm, tag, upv, pv1k, mask1k, elev10k, show_legend) {
  crs_m <- MAP_CRS[[nm]]
  poly  <- st_transform(prov_poly(nm), crs_m)
  poly_v <- vect(st_transform(prov_poly(nm), crs(pv1k)))
  crop_proj <- function(r, method) {
    r <- mask(crop(r, poly_v), poly_v)
    project(r, crs_m, method = method, res = 1000)
  }
  pv <- crop_proj(pv1k, "bilinear")
  mk <- crop_proj(mask1k, "near")
  pv_df <- as.data.frame(pv, xy = TRUE, na.rm = TRUE); names(pv_df)[3] <- "pv"
  mk_df <- as.data.frame(mk, xy = TRUE, na.rm = TRUE); names(mk_df)[3] <- "ok"
  ex_df <- mk_df[mk_df$ok == 0, ]
  el <- project(mask(crop(elev10k, vect(st_transform(prov_poly(nm), crs(elev10k)))),
                     vect(st_transform(prov_poly(nm), crs(elev10k)))), crs_m)
  lv <- seq(1000, 6000, 1000)
  ct <- st_as_sf(as.contour(el, levels = lv[lv > minmax(el)[1] & lv < minmax(el)[2]]))
  fp <- upv[st_intersects(upv, st_transform(prov_poly(nm), st_crs(upv)), sparse = FALSE)[, 1], ]
  fp <- st_buffer(st_transform(st_centroid(fp), crs_m), 3000)
  city <- st_transform(st_as_sf(CITIES[CITIES$prov == nm, ], coords = c("lon", "lat"),
                                crs = 4326), crs_m)
  bb <- st_bbox(poly)
  sb_len <- 200000; sb_x <- bb[["xmin"]] + 0.04 * (bb[["xmax"]] - bb[["xmin"]])
  sb_y <- bb[["ymin"]] + 0.04 * (bb[["ymax"]] - bb[["ymin"]])
  ggplot() +
    geom_raster(data = pv_df, aes(x, y, fill = pv)) +
    geom_raster(data = ex_df, aes(x, y), fill = GREY_MID, alpha = 0.85) +
    geom_sf(data = ct, colour = "white", linewidth = 0.15, alpha = 0.7) +
    geom_sf(data = fp, aes(colour = "UPV plant (3 km buffer)"), fill = OI_BLACK,
            linewidth = 0.1) +
    geom_sf(data = poly, fill = NA, colour = OI_BLACK, linewidth = LW_MED) +
    geom_sf(data = city, shape = 21, fill = "white", colour = OI_BLACK, size = 1.6,
            stroke = 0.5) +
    geom_sf_text(data = city, aes(label = CITIES$name[CITIES$prov == nm]),
                 nudge_y = 45000, size = PT_TEXT / .pt, fontface = "bold") +
    annotate("rect", xmin = sb_x, xmax = sb_x + sb_len, ymin = sb_y,
             ymax = sb_y + 0.012 * (bb[["ymax"]] - bb[["ymin"]]), fill = OI_BLACK) +
    annotate("text", x = sb_x + sb_len / 2, y = sb_y, vjust = -0.6,
             label = "200 km", size = PT_TEXT / .pt) +
    annotate("text", x = bb[["xmax"]], y = bb[["ymax"]], hjust = 1, vjust = 1,
             label = prov_en(nm), size = PT_TITLE / .pt, fontface = "bold") +
    scale_fill_viridis_c(option = "cividis", limits = PV_LIMS, oob = squish,
                         breaks = c(1000, 1400, 1800, 2200),
                         labels = label_comma(), name = lab_pvout,
                         guide = guide_colourbar(barwidth = unit(110, "pt"),
                                                 barheight = unit(5, "pt"),
                                                 title.position = "top")) +
    scale_colour_manual(values = c("UPV plant (3 km buffer)" = OI_BLACK), name = NULL) +
    coord_sf(crs = crs_m, datum = NA, expand = FALSE) +
    labs(tag = tag) + theme_map() +
    theme(legend.position = if (show_legend) "bottom" else "none")
}

fig4 <- function() {
  provs2 <- c("Qinghai", "Xinjiang")
  rng <- c(800, 2300)
  ymax <- 0
  for (nm in provs2) {
    g <- grid_df |> filter(feasible, is.finite(pvout_annual), prov == nm)
    p <- plants |> filter(prov == nm)
    ymax <- max(ymax, dens(g$pvout_annual, rng = rng)$y,
                dens(p$pvout, p$cap_mw, rng)$y)
  }
  yl <- c(0, ymax * 1.12)
  a <- province_density("Qinghai", rng, yl, "a")
  b <- province_density("Xinjiang", rng, yl, "b") + theme(legend.position = "none")
  upv  <- readRDS(file.path(P3, "UPV_core_GRW_2024Q2_ge1MW.rds"))
  pv1k <- rast(file.path(P4, "data/derived/pvout_1km_cea.tif"))
  mk1k <- rast(file.path(P4, "data/derived/topo_mask_1km_cea.tif")); crs(mk1k) <- crs(pv1k)
  el10 <- rast(file.path(P3, "elevation_10km.tif"))
  c_ <- map_panel("Qinghai", "c", upv, pv1k, mk1k, el10, TRUE)
  d  <- map_panel("Xinjiang", "d", upv, pv1k, mk1k, el10, TRUE)
  (a | b) / ((c_ | d) + plot_layout(guides = "collect") &
               theme(legend.position = "bottom")) +
    plot_layout(heights = c(0.75, 1.25))
}

# ============================================================== Figure 5 ====
# Feasible 10 km land classified on the two thresholds of equation 11:
# PVOUT median of feasible land and the plant-median transmission penalty
# (11.4), with distance to the nearest >=220 kV substation in the 2025 snapshot.
nearest_m <- function(xy, nxy) {
  out <- numeric(nrow(xy))
  for (i in seq_len(nrow(xy)))
    out[i] <- sqrt(min((nxy[, 1] - xy[i, 1])^2 + (nxy[, 2] - xy[i, 2])^2))
  out
}
fig5_data <- function() {
  f <- file.path(OUT, "fig5_land.rds")
  if (file.exists(f)) return(readRDS(f))
  pvr <- rast(file.path(P3, "pvout_annual_10km.tif"))
  # Taiwan has no substations in the Fu et al. dataset, so it is left unclassified
  g <- grid_df |> filter(feasible, is.finite(pvout_annual), prov != "Taiwan")
  nd <- read.csv(file.path(P4, "transmission", "grid_nodes_2025.csv"))
  nd <- nd[is.finite(nd$lng) & is.finite(nd$lat), ]
  # great-circle distance: cell centres to lon/lat, haversine to nodes
  ll <- st_coordinates(st_transform(st_as_sf(g, coords = c("x", "y"), crs = crs(pvr)), 4326))
  R <- 6371; la <- ll[, 2] * pi / 180; lo <- ll[, 1] * pi / 180
  nla <- nd$lat * pi / 180; nlo <- nd$lng * pi / 180; cn <- cos(nla)
  d <- vapply(seq_along(la), function(i) {
    a <- sin((nla - la[i]) / 2)^2 + cos(la[i]) * cn * sin((nlo - lo[i]) / 2)^2
    min(2 * R * asin(pmin(1, sqrt(a)))) }, numeric(1))
  g$dist_km <- pmin(d, 600)
  g$tx <- (0.0002 * g$dist_km / 0.37) * g$pvout_annual
  g$quadrant <- factor(with(g, ifelse(pvout_annual >= PVOUT_MED & tx < TX_MED, "Premium",
                         ifelse(pvout_annual >= PVOUT_MED, "Stranded",
                         ifelse(tx < TX_MED, "Grid-convenient", "Dominated")))),
                       levels = QUAD_LEVELS)
  saveRDS(g, f); g
}

fig5 <- function() {
  g <- fig5_data()
  pvr <- rast(file.path(P3, "pvout_annual_10km.tif"))
  shares <- round(100 * table(g$quadrant) / nrow(g), 1)
  cat("  Fig 5 land shares (% of feasible cells):",
      paste(names(shares), shares, collapse = "; "), "\n")
  cat(sprintf("  Premium land: %.1f%% of feasible area, %s km2\n", shares[["Premium"]],
              fmt(sum(g$quadrant == "Premium") * 100)))
  bp <- function(q) sort(round(100 * prop.table(table(prov_en(g$prov[g$quadrant == q]))), 1),
                         decreasing = TRUE)
  summ <- list(shares = shares, premium_km2 = sum(g$quadrant == "Premium") * 100,
               premium_by_prov = bp("Premium"), stranded_by_prov = bp("Stranded"))
  print(head(summ$premium_by_prov, 8)); print(head(summ$stranded_by_prov, 8))
  saveRDS(summ, file.path(OUT, "fig5_summary.rds"))
  crs_cn <- "+proj=aea +lat_1=25 +lat_2=47 +lat_0=0 +lon_0=105 +datum=WGS84 +units=m"
  r <- rast(pvr); r[] <- NA
  r[cellFromXY(r, cbind(g$x, g$y))] <- as.integer(g$quadrant)
  r <- project(r, crs_cn, method = "near", res = 10000)
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE); names(df)[3] <- "q"
  df$quadrant <- factor(QUAD_LEVELS[df$q], levels = QUAD_LEVELS)
  pr  <- st_transform(provs, crs_cn)
  nat <- st_transform(basemap$china_nation, crs_cn)
  ndl <- st_transform(basemap$nine_dotted_line, crs_cn)
  base_layers <- list(
    geom_sf(data = pr, fill = "grey92", colour = NA),
    geom_raster(data = df, aes(x, y, fill = quadrant)),
    geom_sf(data = pr, fill = NA, colour = "white", linewidth = 0.15),
    geom_sf(data = nat, fill = NA, colour = GREY_DARK, linewidth = LW_MED),
    geom_sf(data = ndl, colour = GREY_DARK, linewidth = LW_MED))
  # densify the lon/lat box before projecting, so curved parallels are respected
  box <- function(b) st_bbox(st_transform(st_segmentize(st_as_sfc(st_bbox(b, crs = 4326)), 0.5), crs_cn))
  lim_main <- st_bbox(st_transform(nat, crs_cn))
  lim_main[["ymin"]] <- box(c(xmin = 105, ymin = 17.8, xmax = 112, ymax = 25))[["ymin"]]
  lim_main[["xmax"]] <- lim_main[["xmax"]] + 0.10 * (lim_main[["xmax"]] - lim_main[["xmin"]])
  lim_scs  <- box(c(xmin = 106, ymin = 3, xmax = 123, ymax = 24))
  main <- ggplot() + base_layers +
    scale_fill_manual(values = QUAD_COLS, name = NULL, drop = FALSE,
                      guide = guide_legend(nrow = 1, override.aes = list(colour = NA))) +
    coord_sf(crs = crs_cn, datum = NA, expand = FALSE,
             xlim = lim_main[c("xmin", "xmax")], ylim = lim_main[c("ymin", "ymax")]) +
    theme_map() + theme(legend.position = "bottom",
                        panel.border = element_rect(fill = NA, colour = NA))
  inset <- ggplot() + base_layers +
    scale_fill_manual(values = QUAD_COLS, guide = "none") +
    coord_sf(crs = crs_cn, datum = NA, expand = FALSE,
             xlim = lim_scs[c("xmin", "xmax")], ylim = lim_scs[c("ymin", "ymax")]) +
    theme_map() + theme(panel.border = element_rect(fill = "white", colour = GREY_DARK,
                                                    linewidth = LW_THIN),
                        plot.margin = margin(0, 0, 0, 0))
  main + inset_element(inset, left = 0.845, bottom = 0.0, right = 1.0, top = 0.30,
                       align_to = "panel")
}

# ------------------------------------------------------------------ main ----
args <- commandArgs(trailingOnly = TRUE)
which_figs <- if (length(args)) as.integer(args) else 1:5
if (1 %in% which_figs) save_pdf(fig1(), 1, 75)
if (2 %in% which_figs) save_pdf(fig2(), 2, 110)
if (3 %in% which_figs) {
  p3 <- fig3(); save_pdf(p3, 3, 150)
  sr <- attr(p3, "series")
  saveRDS(sr, file.path(OUT, "fig3_series.rds"))
}
if (4 %in% which_figs) save_pdf(fig4(), 4, 160)
if (5 %in% which_figs) save_pdf(fig5(), 5, 135)
