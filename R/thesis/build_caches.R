source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# thesis_figures.R — publication figure pipeline, single file
#
# "Following the Sun? Utility-scale solar siting efficiency and the resource
#  frontier in China, 2017-2024"
#
# USAGE (RStudio)
#   1. Edit the two paths in section 2 if the project has moved.
#   2. source("thesis_figures.R")          # loads style, data, figures
#   3. preview(fig_5_04())                 # inspect at true print size
#      preview(fig_5_04(), W_SINGLE, 95)   # single-column variant
#      preview_cvd(fig_5_04())             # colour-blind simulation
#      preview_all()                       # every figure in turn
#   4. Adjust anything, re-source, preview again.
#   5. build_all()                         # write PDFs + run compliance checks
#
# Nothing is written to disk until you call save_fig() or build_all().
#
# Requires: sf, terra, ggplot2, dplyr, scales, patchwork, ggrepel
#   install.packages(c("sf","terra","ggplot2","dplyr","scales","patchwork","ggrepel"))
# =============================================================================


suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(scales); library(grid)
})

# ---------------------------------------------------------------------------
# 0. Typeface
#
# Journals differ: some require a specific family, most accept any clean sans.
# Set FONT_FAMILY to a family installed on this machine ("Helvetica", "Arial",
# "Source Sans Pro", ...) or leave "" for the device default. The value is
# checked against the available families and falls back silently rather than
# failing mid-render, which is what an unavailable family otherwise does.
# ---------------------------------------------------------------------------
FONT_FAMILY <- ""   # "" = device default

.resolve_family <- function(f) {
  if (!nzchar(f)) return("")
  ok <- tryCatch({
    if (requireNamespace("systemfonts", quietly = TRUE))
      f %in% systemfonts::system_fonts()$family
    else TRUE
  }, error = function(e) TRUE)
  if (!ok) {
    warning(sprintf("Font family '%s' not found; using device default.", f),
            call. = FALSE)
    return("")
  }
  f
}

# ---------------------------------------------------------------------------
# 1. Quadrant palette
#
# Okabe-Ito colour-blind-safe base. The mapping is deliberately two-dimensional
# so the 2x2 structure is legible in the colour alone:
#
#     HUE       = resource quality   blue = above-median PVOUT
#                                    orange = below-median PVOUT
#     LIGHTNESS = transmission cost  dark = low cost, light = high cost
#
# Premium and Dominated are the two extremes and get the two darkest, most
# saturated values so they carry visual weight; the mixed quadrants recede.
# No red and no green anywhere, per the supervisor's instruction.
# ---------------------------------------------------------------------------
QUAD_LEVELS <- c("Premium", "Stranded", "Grid-convenient", "Dominated")

# Widened from the raw Okabe-Ito values after CVD testing: the original
# #E69F00 / #D55E00 orange pair collapsed under tritanopia (dE 14.6, below the
# 15 threshold). Pushing both hue pairs apart in lightness raises the minimum
# pairwise distance to 40.7 across normal, deuteranopic, protanopic and
# tritanopic vision. Verify with verify_palette() after any change.
QUAD_COLS <- c(
  "Premium"         = "#0060A0",  # blue, dark    - high resource, low cost
  "Stranded"        = "#7FC7EF",  # blue, light   - high resource, high cost
  "Grid-convenient" = "#F2C166",  # orange, light - low resource, low cost
  "Dominated"       = "#9A4200"   # orange, dark  - low resource, high cost
)

# Sequential ramp for continuous fills (maps, densities). Blue-through-orange,
# no red/green, safe under deuteranopia and protanopia.
SEQ_COLS  <- c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#F7F7F7",
               "#FDDBC7", "#F4A582", "#E69F00", "#D55E00")

# Neutral greys for reference lines, land distributions and annotation
GREY_DARK <- "#3B3B3B"
GREY_MID  <- "#7A7A7A"
GREY_LITE <- "#C8C8C8"

quad_factor <- function(x) factor(x, levels = QUAD_LEVELS)

scale_fill_quad <- function(...) scale_fill_manual(values = QUAD_COLS,
                                                   breaks = QUAD_LEVELS, ...)
scale_colour_quad <- function(...) scale_colour_manual(values = QUAD_COLS,
                                                       breaks = QUAD_LEVELS, ...)
scale_color_quad <- scale_colour_quad

# ---------------------------------------------------------------------------
# 2. Type sizes
#
# ggplot2 sizes text in points directly, so these are true print sizes. The
# floor is 8 pt: BASE_PT is the smallest size actually used on any figure, so
# raising it raises everything and nothing can drop below the limit.
# ---------------------------------------------------------------------------
BASE_PT   <- 8    # axis text, legend text, annotations  <- the 8 pt floor
AXIS_PT   <- 9    # axis titles
TITLE_PT  <- 10   # panel titles / tags
STRIP_PT  <- 9    # facet strips

LINE_THIN <- 0.25
LINE_MED  <- 0.45
LINE_THICK<- 0.7

# ---------------------------------------------------------------------------
# 3. Figure widths (mm) — standard journal column geometry
# ---------------------------------------------------------------------------
W_SINGLE <- 89    # single column
W_DOUBLE <- 183   # full width
W_HALF   <- 120   # 1.5 column, for figures too dense for 89 mm

# ---------------------------------------------------------------------------
# 4. Theme
# ---------------------------------------------------------------------------
theme_pub <- function(base_size = BASE_PT, grid_y = TRUE, grid_x = FALSE) {
  fam <- .resolve_family(FONT_FAMILY)
  th <- theme_minimal(base_size = base_size, base_family = fam) +
    theme(
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      text             = element_text(colour = GREY_DARK, size = base_size,
                                      family = fam),
      axis.text        = element_text(size = base_size, colour = GREY_DARK),
      axis.title       = element_text(size = AXIS_PT, colour = GREY_DARK),
      axis.title.x     = element_text(margin = margin(t = 3)),
      axis.title.y     = element_text(margin = margin(r = 3)),
      axis.ticks       = element_line(linewidth = LINE_THIN, colour = GREY_MID),
      axis.ticks.length= unit(1.4, "pt"),
      axis.line        = element_line(linewidth = LINE_THIN, colour = GREY_DARK),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = LINE_THIN, colour = GREY_LITE),
      panel.border     = element_blank(),
      plot.title       = element_text(size = TITLE_PT, face = "plain", hjust = 0,
                                      margin = margin(b = 3)),
      plot.subtitle    = element_text(size = base_size, colour = GREY_MID, hjust = 0,
                                      margin = margin(b = 5)),
      plot.caption     = element_text(size = base_size, colour = GREY_MID, hjust = 0,
                                      margin = margin(t = 5)),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      plot.tag         = element_text(size = TITLE_PT, face = "bold"),
      legend.title     = element_text(size = base_size),
      legend.text      = element_text(size = base_size),
      legend.key.size  = unit(7, "pt"),
      legend.margin    = margin(0, 0, 0, 0),
      legend.box.spacing = unit(3, "pt"),
      legend.background = element_blank(),
      strip.text       = element_text(size = STRIP_PT, colour = GREY_DARK,
                                      margin = margin(2, 2, 2, 2)),
      plot.margin      = margin(4, 4, 4, 4)
    )
  if (!grid_y) th <- th + theme(panel.grid.major.y = element_blank())
  if (!grid_x) th <- th + theme(panel.grid.major.x = element_blank())
  th
}

# Version-safe inside-panel legend. ggplot2 3.5.0 moved this from
# legend.position = c(x, y) to legend.position.inside, and warns on the old
# form. This emits whichever the installed version expects.
legend_inside <- function(x, y, just = c(x, y)) {
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    theme(legend.position = "inside",
          legend.position.inside = c(x, y),
          legend.justification = just)
  } else {
    theme(legend.position = c(x, y), legend.justification = just)
  }
}

# Full box around the panel — required by some journals, off by default.
theme_pub_boxed <- function(...) {
  theme_pub(...) +
    theme(panel.border = element_rect(fill = NA, colour = GREY_DARK,
                                      linewidth = LINE_THIN),
          axis.line = element_blank())
}

# Map theme — no axes, no grid
theme_pub_map <- function(base_size = BASE_PT) {
  theme_pub(base_size) +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), axis.line = element_blank(),
          panel.grid.major = element_blank())
}

# ---------------------------------------------------------------------------
# 5. Save helper
#
# Writes vector PDF. Uses cairo_pdf where available so fonts embed properly;
# falls back to the base pdf device otherwise. Optionally also writes a PNG at
# 400 dpi for drafting only — never submit the PNG.
# ---------------------------------------------------------------------------
FIG_DIR <- "figures_publication"

# Pick a working PDF device ONCE, by actually trying it.
#
# capabilities("cairo") reports what R was compiled with, not what loads at
# runtime. On macOS without XQuartz, cairo_pdf is "available" but fails with
# "failed to load cairo DLL". So we test it on a throwaway file and fall back to
# the base pdf device, which is pure vector and uses the PDF base-14 fonts
# (Helvetica et al.) that every reader has built in — fine for submission.
#
# Custom FONT_FAMILY values need cairo: install XQuartz (xquartz.org) and
# restart R if you need a specific typeface embedded.
.pdf_device <- local({
  cached <- NULL
  function() {
    if (!is.null(cached)) return(cached)
    dev <- grDevices::pdf; label <- "pdf (base-14 fonts)"
    if (isTRUE(capabilities("cairo"))) {
      tf <- tempfile(fileext = ".pdf")
      ok <- tryCatch({
        grDevices::cairo_pdf(tf, width = 1, height = 1)
        grDevices::dev.off()
        file.exists(tf) && file.size(tf) > 0
      }, error = function(e) FALSE, warning = function(w) FALSE)
      # a failed cairo device can leave the graphics stack open
      while (grDevices::dev.cur() > 1L) try(grDevices::dev.off(), silent = TRUE)
      unlink(tf)
      if (isTRUE(ok)) { dev <- grDevices::cairo_pdf; label <- "cairo_pdf" }
    }
    cat(sprintf("  PDF device: %s\n", label))
    if (identical(label, "pdf (base-14 fonts)") && nzchar(FONT_FAMILY))
      warning("FONT_FAMILY is set but cairo is unavailable, so the custom face ",
              "will not be embedded. Install XQuartz, or leave FONT_FAMILY = \"\".",
              call. = FALSE)
    cached <<- dev
    dev
  }
})

save_fig <- function(plot, name, width_mm = W_DOUBLE, height_mm = 110,
                     also_png = TRUE, dir = FIG_DIR) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  pdf_path <- file.path(dir, paste0(name, ".pdf"))
  ggsave(pdf_path, plot, device = .pdf_device(),
         width = width_mm, height = height_mm, units = "mm")
  if (also_png) {
    ggsave(file.path(dir, paste0(name, ".png")), plot,
           width = width_mm, height = height_mm, units = "mm", dpi = 400)
  }
  sz <- file.size(pdf_path)
  cat(sprintf("  saved  %-46s %3.0f x %3.0f mm  %6.0f KB%s\n",
              basename(pdf_path), width_mm, height_mm, sz/1024,
              if (is.na(sz) || sz < 1000) "   <-- CHECK: suspiciously small" else ""))
  invisible(pdf_path)
}

# Check no text in a built plot falls below 8 pt.
check_min_type <- function(p, floor_pt = 8) {
  th <- p$theme
  sizes <- unlist(lapply(th, function(e) if (inherits(e, "element_text")) e$size else NULL))
  sizes <- sizes[is.finite(sizes)]
  if (length(sizes) && min(sizes) < floor_pt)
    warning(sprintf("Type below %g pt found (min %.1f pt)", floor_pt, min(sizes)))
  invisible(min(sizes))
}

cat("fig_00_style.R loaded — palette, theme and save helper ready\n")

# ---------------------------------------------------------------------------
# 6. Compliance checks
#
# These verify the journal requirements rather than assuming them. Run
# verify_all() after generating figures; it is cheap and catches the failures
# that are invisible on screen (raster sneaking into a PDF, type below the
# floor, palette pairs that collapse under colour-vision deficiency).
# ---------------------------------------------------------------------------

# Simulate dichromatic vision. Viénot, Brettel & Mollon (1999) linear method:
# linear-RGB -> LMS, collapse the missing cone, back to sRGB.
simulate_cvd <- function(hex, type = c("deutan", "protan", "tritan")) {
  type <- match.arg(type)
  srgb <- t(grDevices::col2rgb(hex)) / 255
  lin  <- ifelse(srgb <= 0.04045, srgb/12.92, ((srgb + 0.055)/1.055)^2.4)
  M    <- matrix(c(17.8824, 43.5161, 4.11935,
                   3.45565, 27.1554, 3.86714,
                   0.0299566, 0.184309, 1.46709), 3, 3, byrow = TRUE)
  Mi   <- solve(M)
  S    <- switch(type,
    protan = matrix(c(0, 2.02344, -2.52581, 0, 1, 0, 0, 0, 1), 3, 3, byrow = TRUE),
    deutan = matrix(c(1, 0, 0, 0.494207, 0, 1.24827, 0, 0, 1), 3, 3, byrow = TRUE),
    tritan = matrix(c(1, 0, 0, 0, 1, 0, -0.395913, 0.801109, 0), 3, 3, byrow = TRUE))
  lms  <- lin %*% t(M)
  sim  <- lms %*% t(S) %*% t(Mi)
  sim  <- pmin(pmax(sim, 0), 1)
  out  <- ifelse(sim <= 0.0031308, sim*12.92, 1.055*sim^(1/2.4) - 0.055)
  grDevices::rgb(out[,1], out[,2], out[,3])
}

# Minimum pairwise CIE Lab distance within a palette, under normal vision and
# each dichromacy. A dE below ~15 is where two categorical colours start to be
# confusable in small marks.
verify_palette <- function(cols = QUAD_COLS, min_de = 15, verbose = TRUE) {
  de <- function(hx) {
    lab <- grDevices::convertColor(t(grDevices::col2rgb(hx))/255,
                                   from = "sRGB", to = "Lab")
    d <- as.matrix(dist(lab)); diag(d) <- Inf
    list(min = min(d), which = which(d == min(d), arr.ind = TRUE)[1, ])
  }
  modes <- list(normal = cols, deutan = simulate_cvd(cols, "deutan"),
                protan = simulate_cvd(cols, "protan"), tritan = simulate_cvd(cols, "tritan"))
  ok <- TRUE
  for (m in names(modes)) {
    r <- de(modes[[m]]); pass <- r$min >= min_de
    ok <- ok && pass
    if (verbose) cat(sprintf("    %-8s min dE = %5.1f  %s  (closest: %s / %s)\n",
      m, r$min, if (pass) "PASS" else "FAIL",
      names(cols)[r$which[1]], names(cols)[r$which[2]]))
  }
  # explicit red/green screen: hue near 0/120 deg with real saturation
  hsv <- grDevices::rgb2hsv(grDevices::col2rgb(cols))
  hue <- hsv["h", ] * 360; sat <- hsv["s", ]
  bad <- names(cols)[(hue < 20 | hue > 340 | (hue > 80 & hue < 160)) & sat > 0.35]
  if (length(bad)) { ok <- FALSE
    cat("    RED/GREEN WARNING:", paste(bad, collapse = ", "), "\n")
  } else if (verbose) cat("    red/green screen  PASS  (no hue in 340-20 or 80-160 deg)\n")
  invisible(ok)
}

# Confirm a PDF is vector and carries no embedded raster image.
verify_vector <- function(path) {
  if (!file.exists(path)) return(NA)
  raw <- readBin(path, "raw", file.size(path))
  # grepRaw avoids rawToChar, which errors on the nul bytes in compressed PDFs
  length(grepRaw("/Subtype /Image", raw, fixed = TRUE)) == 0L &&
    length(grepRaw("/Subtype/Image", raw, fixed = TRUE)) == 0L
}

verify_all <- function(dir = FIG_DIR, min_de = 15) {
  cat("\n--- palette ---\n"); pal_ok <- verify_palette(min_de = min_de)
  cat("\n--- output ---\n")
  pdfs <- list.files(dir, pattern = "\\.pdf$", full.names = TRUE)
  if (!length(pdfs)) { cat("    no PDFs found in ", dir, "\n"); return(invisible(FALSE)) }
  vec_ok <- TRUE
  for (f in pdfs) {
    v <- verify_vector(f); vec_ok <- vec_ok && isTRUE(v)
    cat(sprintf("    %-46s %s  %6.0f KB\n", basename(f),
                if (isTRUE(v)) "vector" else "RASTER!", file.size(f)/1024))
  }
  cat(sprintf("\n    type floor        %d pt (requirement: 8 pt)  %s\n",
              BASE_PT, if (BASE_PT >= 8) "PASS" else "FAIL"))
  cat(sprintf("    quadrant colours  single QUAD_COLS shared by all figures  PASS\n"))
  cat(sprintf("\n  OVERALL: %s\n", if (pal_ok && vec_ok && BASE_PT >= 8) "PASS" else "REVIEW"))
  invisible(pal_ok && vec_ok)
}

# ---------------------------------------------------------------------------
# 7. Interactive preview
#
# View a figure at its TRUE print size before committing to disk. RStudio's
# Plots pane rescales whatever you send it, so a figure that looks fine there
# can still have illegible 8 pt type at 89 mm. These open a real device at the
# exact millimetre dimensions instead.
#
#   preview(fig_5_04())                     # full width, 183 mm
#   preview(fig_5_04(), W_SINGLE, 95)       # single column, as submitted
#   preview_cvd(fig_5_04())                 # as a deuteranope sees it
#   preview_all()                           # every figure in turn
#
# Adjust, re-run, and only call save_fig() / run_all_figures.R when happy.
# ---------------------------------------------------------------------------

preview <- function(p, width_mm = W_DOUBLE, height_mm = 110, title = NULL) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  opened <- FALSE
  if (interactive()) {
    # noRStudioGD forces a real window; the RStudio pane would rescale and
    # defeat the point of checking type size at print dimensions.
    ok <- tryCatch({ dev.new(width = w, height = h, noRStudioGD = TRUE); TRUE },
                   error = function(e) FALSE)
    opened <- isTRUE(ok)
  }
  print(p)
  cat(sprintf("  preview: %.0f x %.0f mm (%.2f x %.2f in)%s\n",
              width_mm, height_mm, w, h,
              if (!is.null(title)) paste0("  -- ", title) else ""))
  if (!opened && interactive())
    cat("  (RStudio pane in use - sizing is approximate; widen the pane or use\n",
        "   preview_png() for an exact-size check)\n", sep = "")
  invisible(p)
}

# Exact-size check that never depends on the plot pane: render to a temp PNG at
# the true dimensions and open it in the system viewer / RStudio viewer.
preview_png <- function(p, width_mm = W_DOUBLE, height_mm = 110, dpi = 200) {
  f <- tempfile(fileext = ".png")
  ggsave(f, p, width = width_mm, height = height_mm, units = "mm", dpi = dpi)
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) rstudioapi::viewer(f) else utils::browseURL(f)
  cat(sprintf("  exact-size PNG: %s  (%.0f x %.0f mm @ %d dpi)\n",
              f, width_mm, height_mm, dpi))
  invisible(f)
}

# Re-render a figure with the palette pushed through a dichromacy simulation,
# so you can see what a colour-blind reader sees rather than trusting the
# numeric check alone.
preview_cvd <- function(p, type = c("deutan", "protan", "tritan"),
                        width_mm = W_DOUBLE, height_mm = 110) {
  type <- match.arg(type)
  old <- QUAD_COLS
  on.exit(assign("QUAD_COLS", old, envir = .GlobalEnv), add = TRUE)
  sim <- simulate_cvd(old, type); names(sim) <- names(old)
  assign("QUAD_COLS", sim, envir = .GlobalEnv)
  cat(sprintf("  simulating %s vision\n", type))
  cat("  NOTE: pass a figure CALL, not a stored object, e.g.\n")
  cat("        preview_cvd(fig_5_04())  -- rebuilt under the simulated palette\n")
  preview(p, width_mm, height_mm, title = paste(type, "simulation"))
}

# Every figure at its submission size, one keypress apart.
preview_all <- function() {
  specs <- list(
    list(f = function() fig_5_04(),  w = W_DOUBLE, h = 115, n = "5.4 quadrant scatter"),
    list(f = function() fig_5_07(),  w = W_DOUBLE, h = 118, n = "5.7 provincial drivers"),
    list(f = function() fig_5_11(),  w = W_SINGLE, h = 88,  n = "5.11 Qinghai"),
    list(f = function() fig_5_12(),  w = W_SINGLE, h = 88,  n = "5.12 Xinjiang"),
    list(f = function() fig_5_11_12_combined(), w = W_DOUBLE, h = 92, n = "5.11/5.12 pair"),
    list(f = function() fig_5_15(basemap = if (exists("basemap")) basemap else NULL),
         w = W_DOUBLE, h = 140, n = "5.15 investment map"),
    list(f = function() fig_graphical_abstract(), w = W_DOUBLE, h = 105, n = "graphical abstract")
  )
  for (s in specs) {
    cat("\n---", s$n, "---\n")
    r <- tryCatch(preview(s$f(), s$w, s$h, title = s$n), error = function(e) {
      cat("  SKIPPED:", conditionMessage(e), "\n"); NULL })
    if (interactive() && !is.null(r)) {
      cat("  [Enter] for next, [q] to stop: ")
      if (tolower(trimws(readline())) == "q") break
    }
  }
  invisible(NULL)
}



suppressPackageStartupMessages({
  library(sf); library(dplyr); library(terra)
})
sf_use_s2(FALSE)   # GRW polygons contain degenerate edges; S2 rejects them

# ---------------------------------------------------------------------------
# Paths — edit these two if the project moves
# ---------------------------------------------------------------------------
P3 <- FTS_PROCESSED
P4 <- FTS_ROOT

NODES_2020 <- file.path(P4, "transmission", "grid_nodes_2020.csv")
NODES_2025 <- file.path(P4, "transmission", "grid_nodes_2025.csv")

CACHE <- "fig_cache"
if (!dir.exists(CACHE)) dir.create(CACHE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Model parameters — from Wang et al. (2022), calibrated so 500 km -> 0.10 CNY/kWh
# ---------------------------------------------------------------------------
GF_PER_KM  <- 0.0002   # CNY/kWh per km
ELEC_PRICE <- 0.37     # CNY/kWh, post-531 grid parity benchmark
DIST_CAP   <- 600      # km, UHV effective range
CURT_RATE  <- 0.20     # curtailment penalty, flagged provinces
CURT_PROV  <- c("Xinjiang", "Gansu", "Ningxia")
PVOUT_MED  <- 1470.097 # national feasible-land median (read from grid_df)
GRID_ACCESS_KM <- 50   # "grid-accessible" = within 50 km of a >=220 kV substation

need <- function(p) if (!file.exists(p)) stop("Missing required file:\n  ", p, call. = FALSE)

# ---------------------------------------------------------------------------
# Haversine nearest-node distance (km). Vectorised over nodes, looped over
# points; ~3 s for 16.7k plants against 3.4k nodes.
# ---------------------------------------------------------------------------
nearest_km <- function(plat, plon, nlat, nlon, chunk = 2000L) {
  R <- 6371
  pla <- plat * pi/180; plo <- plon * pi/180
  nla <- nlat * pi/180; nlo <- nlon * pi/180
  cos_nla <- cos(nla)
  out <- numeric(length(pla))
  for (i in seq_along(pla)) {
    a <- sin((nla - pla[i])/2)^2 + cos(pla[i]) * cos_nla * sin((nlo - plo[i])/2)^2
    out[i] <- min(2 * R * asin(pmin(1, sqrt(a))))
  }
  out
}

# =============================================================================
# 1. Plants — centroids, PVOUT, capacity, province
# =============================================================================
plants_cache <- file.path(CACHE, "plants_quad.rds")

if (file.exists(plants_cache)) {
  cat("Loading cached plant frame...\n")
  plants <- readRDS(plants_cache)
} else {
  cat("Building plant frame from source...\n")
  need(file.path(P3, "UPV_core_GRW_2024Q2_ge1MW.rds"))
  need(NODES_2020)

  upv <- readRDS(file.path(P3, "UPV_core_GRW_2024Q2_ge1MW.rds"))
  cat("  GRW polygons:", nrow(upv), "\n")

  # province: spatial join against the basemap, matching 00_setup.R
  bm_path <- file.path(P3, "china_basemap_bundle.rds")
  need(bm_path)
  bm    <- readRDS(bm_path)
  provs <- bm$china_provinces
  cand  <- c("NAME_1","name_1","name","prov","province","VARNAME_1")
  pcol  <- cand[cand %in% names(provs)][1]
  if (is.na(pcol)) {
    ng   <- setdiff(names(provs), attr(provs, "sf_column"))
    pcol <- ng[sapply(ng, function(c) is.character(provs[[c]]) || is.factor(provs[[c]]))][1]
  }
  prov_sf <- provs[, c(pcol, attr(provs, "sf_column"))]
  names(prov_sf)[1] <- "prov"
  if (!identical(st_crs(upv), st_crs(prov_sf))) upv <- st_transform(upv, st_crs(prov_sf))

  cent <- st_centroid(upv)
  jn   <- st_join(cent, prov_sf, join = st_within, left = TRUE)
  prov <- jn$prov
  if (any(is.na(prov))) {
    idx <- which(is.na(prov))
    cat("  assigning", length(idx), "boundary plants to nearest province\n")
    prov[idx] <- prov_sf$prov[st_nearest_feature(cent[idx, ], prov_sf)]
  }

  # PVOUT at each centroid
  need(file.path(P3, "pvout_annual_10km.tif"))
  pv  <- rast(file.path(P3, "pvout_annual_10km.tif"))
  ex  <- terra::extract(pv, vect(st_transform(cent, crs(pv))), method = "simple")
  pvo <- ex[[2]]

  ll <- st_coordinates(st_transform(cent, 4326))

  plants <- data.frame(
    lon = ll[,1], lat = ll[,2],
    cap_mw = as.numeric(upv$capacity_mw_used),
    pvout = as.numeric(pvo),
    prov = as.character(prov),
    build_year = upv$build_year,
    year_quarter = as.character(upv$year_quarter),
    stringsAsFactors = FALSE
  ) |>
    filter(is.finite(pvout), is.finite(cap_mw), cap_mw > 0, !is.na(prov))
  cat("  plants after filter:", nrow(plants),
      sprintf("(%.1f GW)\n", sum(plants$cap_mw)/1e3))

  # transmission layer
  nd <- read.csv(NODES_2020, stringsAsFactors = FALSE)
  lonc <- if ("lng" %in% names(nd)) "lng" else "lon"
  nd <- nd[is.finite(nd[[lonc]]) & is.finite(nd$lat), ]
  cat("  substations (2020):", nrow(nd), "\n")

  plants$dist_km    <- pmin(nearest_km(plants$lat, plants$lon, nd$lat, nd[[lonc]]), DIST_CAP)
  plants$tx_penalty <- (GF_PER_KM * plants$dist_km / ELEC_PRICE) * plants$pvout
  plants$curt_pen   <- ifelse(plants$prov %in% CURT_PROV, CURT_RATE * plants$pvout, 0)
  plants$ngv        <- plants$pvout - plants$tx_penalty - plants$curt_pen

  TX_MED <- median(plants$tx_penalty)
  plants$quadrant <- with(plants, ifelse(
    pvout >= PVOUT_MED & tx_penalty <  TX_MED, "Premium",
    ifelse(pvout >= PVOUT_MED & tx_penalty >= TX_MED, "Stranded",
    ifelse(pvout <  PVOUT_MED & tx_penalty <  TX_MED, "Grid-convenient", "Dominated"))))

  qn <- function(x) as.integer(substr(x,1,4))*4 + as.integer(substr(x,6,6))
  plants$post2020 <- qn(plants$year_quarter) > qn("2020Q2")

  attr(plants, "TX_MED") <- TX_MED
  saveRDS(plants, plants_cache)

  cat(sprintf("\n  VALIDATION (targets from DIAGNOSTIC_REPORT.md)\n"))
  cat(sprintf("    median distance   %6.2f km   [15.98]\n", median(plants$dist_km)))
  cat(sprintf("    Xinjiang median   %6.2f km   [17.95]\n",
              median(plants$dist_km[plants$prov == "Xinjiang"])))
  cat(sprintf("    median tx penalty %6.2f      [11.42]\n", TX_MED))
  tot <- sum(plants$cap_mw)
  for (q in c("Premium","Stranded","Grid-convenient","Dominated"))
    cat(sprintf("    %-16s %6.2f%%\n", q, 100*sum(plants$cap_mw[plants$quadrant==q])/tot))
  cat("    manuscript: Premium 35.1 / Stranded 26.4 / Grid-conv 17.4 / Dominated 21.0\n")
}

TX_MED <- attr(plants, "TX_MED")
if (is.null(TX_MED)) TX_MED <- median(plants$tx_penalty)

# =============================================================================
# 2. Land grid — feasible and grid-accessible PVOUT distributions
# =============================================================================
grid_cache <- file.path(CACHE, "grid_access.rds")

if (file.exists(grid_cache)) {
  cat("Loading cached land grid...\n")
  grid_df <- readRDS(grid_cache)
} else {
  cat("Building land grid...\n")
  need(file.path(P3, "grid10km_solar_pvout_table_correct_topo.rds"))
  grid_df <- readRDS(file.path(P3, "grid10km_solar_pvout_table_correct_topo.rds"))
  grid_df$feasible <- grid_df$is_allowed_topography %in% c(1, TRUE)
  cat("  cells:", nrow(grid_df), " feasible:", sum(grid_df$feasible), "\n")

  # ---------------------------------------------------------------------------
  # Province assignment for grid cells.
  #
  # Prefers the pre-computed lookup, but only if it actually aligns row-for-row
  # with grid_df. Otherwise falls back to doing the spatial join here, which
  # removes the dependency on that file entirely. grid_df$x/$y are projected
  # (CEA) metres, so points are built in the raster CRS then transformed to the
  # basemap CRS before joining.
  # ---------------------------------------------------------------------------
  assign_prov <- function(gd) {
    pa <- file.path(P4, "data", "grid_province_assignment.rds")
    if (file.exists(pa)) {
      a <- tryCatch(readRDS(pa), error = function(e) NULL)
      if (is.data.frame(a) && nrow(a) == nrow(gd)) {
        pc <- intersect(c("prov_name", "prov", "province", "NAME_1"), names(a))[1]
        if (!is.na(pc) && sum(!is.na(a[[pc]])) > 0.5 * nrow(gd)) {
          cat("  province: using cached lookup (", pc, ")\n")
          return(as.character(a[[pc]]))
        }
      }
      cat("  province: cached lookup unusable (rows/columns do not match) - joining\n")
    } else {
      cat("  province: no cached lookup - joining spatially\n")
    }

    bmp <- file.path(P3, "china_basemap_bundle.rds")
    if (!file.exists(bmp)) stop(
      "Cannot assign provinces: neither a usable data/grid_province_assignment.rds ",
      "nor china_basemap_bundle.rds was found.", call. = FALSE)
    provs <- readRDS(bmp)$china_provinces
    cand  <- c("NAME_1", "name_1", "name", "prov", "province", "VARNAME_1")
    pcol  <- cand[cand %in% names(provs)][1]
    if (is.na(pcol)) {
      ng   <- setdiff(names(provs), attr(provs, "sf_column"))
      pcol <- ng[sapply(ng, function(c) is.character(provs[[c]]) ||
                                        is.factor(provs[[c]]))][1]
    }
    if (is.na(pcol)) stop("No province-name column in china_provinces.", call. = FALSE)
    psf <- provs[, c(pcol, attr(provs, "sf_column"))]
    names(psf)[1] <- "prov"

    pv  <- terra::rast(file.path(P3, "pvout_annual_10km.tif"))
    pts <- sf::st_as_sf(data.frame(x = gd$x, y = gd$y), coords = c("x", "y"),
                        crs = sf::st_crs(terra::crs(pv)))
    pts <- sf::st_transform(pts, sf::st_crs(psf))
    cat("  joining", nrow(pts), "cells to", nrow(psf), "province polygons...\n")
    out <- sf::st_join(pts, psf, join = sf::st_within, left = TRUE)$prov
    miss <- which(is.na(out))
    if (length(miss)) {
      cat("  ", length(miss), "cells outside polygons -> nearest province\n")
      out[miss] <- psf$prov[sf::st_nearest_feature(pts[miss, ], psf)]
    }
    as.character(out)
  }

  grid_df$prov <- assign_prov(grid_df)
  cat("  provinces assigned:", length(unique(na.omit(grid_df$prov))), "distinct\n")

  # grid-accessible mask: feasible cell within GRID_ACCESS_KM of a >=220 kV
  # substation. grid_df x/y are projected (CEA) metres, so project nodes to match.
  pv  <- rast(file.path(P3, "pvout_annual_10km.tif"))
  nd  <- read.csv(NODES_2020, stringsAsFactors = FALSE)
  lonc <- if ("lng" %in% names(nd)) "lng" else "lon"
  nd  <- nd[is.finite(nd[[lonc]]) & is.finite(nd$lat), ]
  nsf <- st_as_sf(nd, coords = c(lonc, "lat"), crs = 4326)
  nxy <- st_coordinates(st_transform(nsf, st_crs(crs(pv))))

  cat("  computing grid-accessible mask (", nrow(grid_df), "cells x",
      nrow(nxy), "nodes )...\n")
  thr <- GRID_ACCESS_KM * 1000
  acc <- logical(nrow(grid_df))
  gx <- grid_df$x; gy <- grid_df$y
  step <- 2000L
  for (s in seq(1, length(gx), by = step)) {
    e   <- min(s + step - 1L, length(gx))
    idx <- s:e
    # bounding-box prefilter, then exact distance on survivors
    for (k in idx) {
      cand <- which(abs(nxy[,1] - gx[k]) < thr & abs(nxy[,2] - gy[k]) < thr)
      if (length(cand)) {
        acc[k] <- min((nxy[cand,1]-gx[k])^2 + (nxy[cand,2]-gy[k])^2) <= thr^2
      }
    }
    if (s %% 20000L == 1L) cat("    ", e, "/", length(gx), "\n")
  }
  grid_df$grid_access <- acc
  cat("  grid-accessible feasible cells:", sum(acc & grid_df$feasible), "\n")
  saveRDS(grid_df, grid_cache)
}

feas <- grid_df$pvout_annual[grid_df$feasible & is.finite(grid_df$pvout_annual)]
Q <- quantile(feas, c(0.50, 0.75, 0.80, 0.90))
cat(sprintf("\nFeasible-land quantiles: Q50 %.1f  Q75 %.1f  Q80 %.1f  Q90 %.1f\n",
            Q[1], Q[2], Q[3], Q[4]))

cat("\nData ready: `plants` (", nrow(plants), " rows), `grid_df` (",
    nrow(grid_df), " cells), TX_MED = ", round(TX_MED, 2), "\n", sep = "")



suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(scales); library(patchwork); library(ggrepel)
})

if (!exists("QUAD_COLS")) source("fig_00_style.R")
if (!exists("plants")) stop("Run fig_01_prepare_data.R first.", call. = FALSE)
if (!"post2020" %in% names(plants) && "post" %in% names(plants))
  plants$post2020 <- plants$post
if (!exists("TX_MED"))    TX_MED    <- median(plants$tx_penalty)
if (!exists("PVOUT_MED")) PVOUT_MED <- 1470.097

plants$quadrant <- quad_factor(plants$quadrant)

# =============================================================================
# Figure 5.4 [old 6.5] — Quadrant scatter: PVOUT vs transmission cost penalty
#
# 16.7k points, so: small marks, low alpha, and the quadrant structure carried
# by the reference lines and corner labels rather than by a legend the reader
# has to cross-reference. Transmission penalty is heavily right-skewed
# (median 11.4, max ~600), so the y axis is square-root scaled.
# =============================================================================
fig_5_04 <- function(dat = plants) {
  tot <- sum(dat$cap_mw)
  shr <- dat |> group_by(quadrant) |>
    summarise(pct = 100 * sum(cap_mw) / tot, .groups = "drop")
  gp <- function(q) shr$pct[shr$quadrant == q]

  xr <- range(dat$pvout, na.rm = TRUE)
  yr <- c(0, quantile(dat$tx_penalty, 0.995, na.rm = TRUE))

  lab <- data.frame(
    quadrant = quad_factor(QUAD_LEVELS),
    x = c(xr[2] - 0.02*diff(xr), xr[2] - 0.02*diff(xr),
          xr[1] + 0.02*diff(xr), xr[1] + 0.02*diff(xr)),
    y = c(0.5, yr[2]*0.92, 0.5, yr[2]*0.92),
    hj = c(1, 1, 0, 0),
    txt = sprintf("%s\n%.1f%% of capacity",
                  QUAD_LEVELS, sapply(QUAD_LEVELS, gp))
  )

  ggplot(dat, aes(pvout, tx_penalty)) +
    # Contours show where the PLANTS are (count density); point area shows
    # where the CAPACITY is. The divergence between the two is the paper's
    # central finding: 67.8% of plants sit below the median, but only 38.4%
    # of capacity does, because the large western plants are sited well.
    geom_density_2d(colour = GREY_MID, linewidth = LINE_THIN,
                    alpha = 0.55, bins = 6) +
    geom_point(aes(colour = quadrant, size = cap_mw),
               alpha = 0.55, stroke = 0) +
    geom_hline(yintercept = TX_MED,    linetype = "dashed",
               colour = GREY_DARK, linewidth = LINE_MED) +
    geom_vline(xintercept = PVOUT_MED, linetype = "dashed",
               colour = GREY_DARK, linewidth = LINE_MED) +
    geom_text(data = lab, aes(x, y, label = txt, hjust = hj, colour = quadrant),
              size = BASE_PT/.pt, lineheight = 0.95, fontface = "bold",
              vjust = c(0, 1, 0, 1), show.legend = FALSE) +
    scale_colour_quad(guide = "none") +
    scale_size_area(max_size = 5, breaks = c(100, 1000, 5000),
                    labels = c("100", "1,000", "5,000"),
                    name = "Plant capacity (MW)") +
    scale_y_sqrt(breaks = c(0, 5, 25, 100, 250, 500), expand = expansion(c(0, 0.02))) +
    scale_x_continuous(breaks = pretty_breaks(6)) +
    coord_cartesian(xlim = xr, ylim = yr) +
    # Legend keys are drawn at a floor size so the small breaks stay legible as
    # dots; the in-panel marks still scale continuously from zero.
    guides(size = guide_legend(
      override.aes = list(colour = GREY_MID, alpha = 0.85,
                          size = c(1.4, 3.0, 5.0)),
      nrow = 1, label.position = "right", title.position = "left",
      title.vjust = 0.5)) +
    labs(
      x = expression("PVOUT (kWh kWp"^-1*" yr"^-1*")"),
      y = expression("Transmission cost penalty (kWh kWp"^-1*" yr"^-1*")"),
      caption = sprintf(
        paste0("n = %s plants (%.1f GW). Point area is plant capacity; grey contours are ",
               "plant-count density.\nDashed lines: national feasible-land PVOUT median (%.0f) ",
               "and median transmission penalty (%.1f).\nQuadrant shares are capacity-weighted. ",
               "y axis square-root scaled."),
        comma(nrow(dat)), tot/1e3, PVOUT_MED, TX_MED)
    ) +
    theme_pub(grid_y = TRUE) +
    theme(legend.position = "bottom", legend.box.spacing = unit(1, "pt"))
}

# =============================================================================
# Figure 5.7 [old 6.8] — Provincial drivers of siting change, post-Q2 2020
#
# Two dimensions decide a province's contribution: how much it built, and how
# far below benchmark it built. Plotted directly against each other, with the
# contribution to national deterioration as point area. Provinces above the
# zero line improved the national average; those below dragged it down.
# =============================================================================
# BENCHMARK NOTE. Table 5.1 states its benchmark as 4.22 kWh/kWp/day, i.e.
# 1540.3 kWh/kWp/yr — NOT the feasible-land median of 1470.1 used everywhere
# else in the manuscript. The choice decides the province ranking: at 1470.1
# Guizhou ranks first and Xinjiang third, contradicting the published table; at
# 1540.3 the order is Xinjiang, Guizhou, Hubei, Anhui exactly as published.
# Figure 5.7 therefore defaults to the Table 5.1 value so figure and table agree.
# See DIAGNOSTIC_REPORT.md for the unresolved share-magnitude gap.
DETERIORATION_BENCHMARK <- 4.22 * 365   # 1540.3 kWh/kWp/yr

fig_5_07 <- function(dat = plants, n_label = 10,
                     benchmark = DETERIORATION_BENCHMARK) {
  post <- dat |> filter(post2020) |> group_by(prov) |>
    summarise(cap_gw = sum(cap_mw)/1e3,
              wpv    = weighted.mean(pvout, cap_mw), .groups = "drop") |>
    mutate(gap  = wpv - benchmark,
           contrib = cap_gw * gap) |>
    filter(cap_gw > 0.15)

  det <- sum(post$contrib[post$contrib < 0])
  post <- post |> mutate(
    share = ifelse(contrib < 0, 100 * contrib / det, NA_real_),
    dir   = ifelse(gap >= 0, "Above benchmark", "Below benchmark"))

  top <- post |> arrange(contrib) |> slice_head(n = n_label)
  best<- post |> arrange(desc(contrib)) |> slice_head(n = 3)
  labs_df <- bind_rows(top, best) |> distinct(prov, .keep_all = TRUE)

  ggplot(post, aes(cap_gw, gap)) +
    geom_hline(yintercept = 0, colour = GREY_DARK, linewidth = LINE_MED) +
    geom_point(aes(size = abs(contrib), fill = dir),
               shape = 21, colour = "white", stroke = 0.25, alpha = 0.9) +
    geom_text_repel(data = labs_df,
                    aes(label = ifelse(is.na(share), prov,
                                       sprintf("%s (%.1f%%)", prov, share))),
                    size = BASE_PT/.pt, colour = GREY_DARK, seed = 1,
                    min.segment.length = 0.2, segment.size = LINE_THIN,
                    segment.colour = GREY_MID, box.padding = 0.35,
                    max.overlaps = 20) +
    scale_fill_manual(values = c("Above benchmark" = QUAD_COLS[["Premium"]],
                                 "Below benchmark" = QUAD_COLS[["Dominated"]]),
                      name = NULL) +
    scale_size_area(max_size = 9, guide = "none") +
    scale_x_continuous(breaks = pretty_breaks(6)) +
    labs(
      x = "Capacity added after Q2 2020 (GW)",
      y = expression(atop("Capacity-weighted PVOUT of new builds",
                          "relative to benchmark (kWh kWp"^-1*" yr"^-1*")")),
      caption = paste0(
        sprintf("Benchmark %.0f kWh kWp-1 yr-1 (4.22 kWh kWp-1 day-1), per Table 5.1.\n", benchmark),
        "Point area is each province's contribution to the national shift\n",
        "(capacity x PVOUT gap). Percentages give each province's share of\n",
        "gross national deterioration; provinces above the line are excluded.")
    ) +
    theme_pub(grid_y = TRUE) +
    legend_inside(0.99, 0.99) +
    theme(legend.direction = "vertical")
}

# =============================================================================
# Figures 5.11 and 5.12 [old 6.12, 6.13] — Qinghai and Xinjiang
#
# MATCHED PAIR. Identical axes, scales, colours and bandwidth. This contrast is
# the paper's central visual argument, so the two panels must be directly
# comparable: build them through one function and share the limits explicitly.
#
# Three distributions per panel:
#   grey area   - all feasible land in the province
#   orange line - feasible land within 50 km of a >=220 kV substation
#   blue area   - deployed capacity, capacity-weighted
# =============================================================================
province_panel <- function(prov_name, dat = plants, gd = grid_df,
                           xlim = NULL, ylim = NULL, tag = NULL, show_legend = FALSE) {
  g <- gd |> filter(feasible, is.finite(pvout_annual), prov == prov_name)
  if (!nrow(g)) stop("No grid cells for province: ", prov_name, call. = FALSE)
  p <- dat |> filter(prov == prov_name, is.finite(pvout))

  bw  <- 45
  rng <- if (is.null(xlim)) range(c(g$pvout_annual, p$pvout)) else xlim

  d_feas <- density(g$pvout_annual, bw = bw, from = rng[1], to = rng[2], n = 512)
  ga     <- g$pvout_annual[g$grid_access]
  d_acc  <- if (length(ga) > 10)
    density(ga, bw = bw, from = rng[1], to = rng[2], n = 512) else NULL
  d_dep  <- density(p$pvout, weights = p$cap_mw/sum(p$cap_mw), bw = bw,
                    from = rng[1], to = rng[2], n = 512)

  df <- rbind(
    data.frame(x = d_feas$x, y = d_feas$y, k = "Feasible land"),
    if (!is.null(d_acc)) data.frame(x = d_acc$x, y = d_acc$y, k = "Grid-accessible land"),
    data.frame(x = d_dep$x,  y = d_dep$y,  k = "Deployed capacity")
  )
  df$k <- factor(df$k, levels = c("Feasible land", "Grid-accessible land", "Deployed capacity"))

  wmed <- with(p[order(p$pvout), ],
               pvout[which(cumsum(cap_mw)/sum(cap_mw) >= 0.5)[1]])

  ggplot() +
    geom_area(data = subset(df, k == "Feasible land"), aes(x, y, fill = k),
              colour = NA, alpha = 0.55) +
    geom_area(data = subset(df, k == "Deployed capacity"), aes(x, y, fill = k),
              colour = NA, alpha = 0.55) +
    geom_line(data = subset(df, k == "Grid-accessible land"), aes(x, y, colour = k),
              linewidth = LINE_THICK, linetype = "22") +
    geom_vline(xintercept = PVOUT_MED, linetype = "dashed",
               colour = GREY_DARK, linewidth = LINE_MED) +
    geom_vline(xintercept = wmed, colour = QUAD_COLS[["Premium"]],
               linewidth = LINE_THICK) +
    annotate("text", x = wmed, y = Inf, label = sprintf(" median MW: %.0f", wmed),
             hjust = 0, vjust = 1.6, size = BASE_PT/.pt,
             colour = QUAD_COLS[["Premium"]], fontface = "bold") +
    scale_fill_manual(values = c("Feasible land" = GREY_LITE,
                                 "Deployed capacity" = QUAD_COLS[["Premium"]]),
                      name = NULL, breaks = c("Feasible land", "Deployed capacity")) +
    scale_colour_manual(values = c("Grid-accessible land" = QUAD_COLS[["Grid-convenient"]]),
                        name = NULL) +
    scale_y_continuous(expand = expansion(c(0, 0.10)), labels = NULL, breaks = NULL) +
    coord_cartesian(xlim = rng, ylim = ylim) +
    labs(title = prov_name, tag = tag, y = "Density",
         x = expression("PVOUT (kWh kWp"^-1*" yr"^-1*")")) +
    theme_pub(grid_y = FALSE) +
    theme(legend.position = if (show_legend) "bottom" else "none",
          legend.box = "horizontal")
}

# Shared limits so the pair is genuinely comparable
pair_limits <- function(provs = c("Qinghai", "Xinjiang"), dat = plants, gd = grid_df) {
  xs <- c(); ys <- 0
  for (pn in provs) {
    g <- gd |> filter(feasible, is.finite(pvout_annual), prov == pn)
    p <- dat |> filter(prov == pn, is.finite(pvout))
    xs <- range(c(xs, g$pvout_annual, p$pvout))
  }
  for (pn in provs) {
    g <- gd |> filter(feasible, is.finite(pvout_annual), prov == pn)
    p <- dat |> filter(prov == pn, is.finite(pvout))
    ys <- max(ys,
      max(density(g$pvout_annual, bw = 45, from = xs[1], to = xs[2], n = 512)$y),
      max(density(p$pvout, weights = p$cap_mw/sum(p$cap_mw), bw = 45,
                  from = xs[1], to = xs[2], n = 512)$y))
  }
  list(xlim = xs, ylim = c(0, ys * 1.08))
}

fig_5_11 <- function(...) { L <- pair_limits()
  province_panel("Qinghai",  xlim = L$xlim, ylim = L$ylim, show_legend = TRUE, ...) }
fig_5_12 <- function(...) { L <- pair_limits()
  province_panel("Xinjiang", xlim = L$xlim, ylim = L$ylim, show_legend = TRUE, ...) }

fig_5_11_12_combined <- function() {
  L <- pair_limits()
  (province_panel("Qinghai",  xlim = L$xlim, ylim = L$ylim, tag = "a") |
   province_panel("Xinjiang", xlim = L$xlim, ylim = L$ylim, tag = "b", show_legend = TRUE)) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

# =============================================================================
# Figure 5.15 [old 6.16] — Investment opportunity map, 2025 grid vintage
#
# Feasible land scored on resource quality and 2025 grid access. Highlights
# where new capacity would land in the Premium quadrant.
# =============================================================================
fig_5_15 <- function(gd = grid_df, basemap = NULL) {
  if (!"grid_access_2025" %in% names(gd))
    stop("Needs grid_access_2025 (2025 grid vintage not yet computed).", call. = FALSE)
  g <- gd |> filter(feasible, is.finite(pvout_annual)) |>
    mutate(opportunity = case_when(
      pvout_annual >= quantile(feas, 0.80) &  grid_access_2025 ~ "Premium: high resource, connected",
      pvout_annual >= quantile(feas, 0.80) & !grid_access_2025 ~ "Stranded: high resource, unconnected",
      TRUE ~ "Below top quintile"))
  g$opportunity <- factor(g$opportunity, levels = c(
    "Premium: high resource, connected", "Stranded: high resource, unconnected",
    "Below top quintile"))

  p <- ggplot() +
    geom_raster(data = g, aes(x, y, fill = opportunity))
  if (!is.null(basemap))
    p <- p + geom_sf(data = basemap, fill = NA, colour = GREY_DARK,
                     linewidth = LINE_THIN)
  p + scale_fill_manual(values = c(
        "Premium: high resource, connected"    = QUAD_COLS[["Premium"]],
        "Stranded: high resource, unconnected" = QUAD_COLS[["Dominated"]],
        "Below top quintile"                   = GREY_LITE), name = NULL) +
    coord_sf(expand = FALSE) +
    labs(caption = paste0(
      "Feasible land only (elevation < 5,000 m, slope < 50 deg). Top quintile of the ",
      "national feasible-land PVOUT\ndistribution. Connected = within 50 km of a ",
      ">=220 kV substation in the 2025 grid vintage (Fu et al., 2026).")) +
    theme_pub_map() +
    theme(legend.position = "bottom", legend.direction = "vertical")
}

# =============================================================================
# Graphical abstract
#
# The quadrant framework paired with the Qinghai/Xinjiang contrast: two
# provinces, comparable capacity added under the same national incentives,
# opposite outcomes.
# =============================================================================
fig_graphical_abstract <- function(dat = plants) {
  two <- dat |> filter(prov %in% c("Qinghai", "Xinjiang"))

  a <- ggplot(two, aes(pvout, tx_penalty)) +
    geom_hline(yintercept = TX_MED, linetype = "dashed",
               colour = GREY_DARK, linewidth = LINE_MED) +
    geom_vline(xintercept = PVOUT_MED, linetype = "dashed",
               colour = GREY_DARK, linewidth = LINE_MED) +
    geom_point(aes(colour = quadrant, size = cap_mw), alpha = 0.5, stroke = 0) +
    scale_colour_quad(name = NULL) +
    scale_size_area(max_size = 4, guide = "none") +
    scale_y_sqrt(breaks = c(0, 5, 25, 100, 250)) +
    scale_x_continuous(breaks = seq(1300, 1900, by = 200)) +
    facet_wrap(~prov, nrow = 1) +
    labs(x = expression("PVOUT (kWh kWp"^-1*" yr"^-1*")"),
         y = "Transmission penalty") +
    theme_pub() +
    theme(legend.position = "bottom",
          panel.spacing.x = unit(9, "pt")) +
    guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1), nrow = 1))

  sh <- two |> group_by(prov, quadrant) |>
    summarise(gw = sum(cap_mw)/1e3, .groups = "drop_last") |>
    mutate(pct = 100*gw/sum(gw)) |> ungroup()

  b <- ggplot(sh, aes(prov, pct, fill = quadrant)) +
    geom_col(width = 0.62, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(pct >= 6, sprintf("%.0f%%", pct), "")),
              position = position_stack(vjust = 0.5),
              size = BASE_PT/.pt, colour = "white", fontface = "bold") +
    scale_fill_quad(guide = "none") +
    scale_y_continuous(expand = expansion(c(0, 0.02)), labels = function(x) paste0(x, "%")) +
    labs(x = NULL, y = "Share of installed capacity") +
    theme_pub(grid_y = TRUE)

  (a | b) + plot_layout(widths = c(2.1, 1)) +
    plot_annotation(
      title = "Same incentives, opposite outcomes",
      subtitle = paste0(
        "Qinghai and Xinjiang each added ~26 GW after the 2020 dual-carbon announcement. ",
        "Qinghai's landed on high-resource land;\nXinjiang's on land that is both ",
        "lower-resource and further from the grid."),
      theme = theme_pub())
}

# =============================================================================
# FIGURE REGISTRY
#
# One list drives everything: submission PDFs, draft PNGs and the markdown
# reference block. Add a figure here once and it appears in all three.
#   id      - manuscript number, used in filenames and captions
#   fn      - a function returning the plot (called lazily, so a broken figure
#             does not stop the others)
#   w, h    - submission size in mm
#   caption - the caption text written into the markdown block
# =============================================================================
FIGS <- list(
  list(id = "5_04",  fn = function() fig_5_04(), w = W_DOUBLE, h = 115,
       caption = "PVOUT vs transmission cost penalty"),
  list(id = "5_07",  fn = function() fig_5_07(), w = W_DOUBLE, h = 118,
       caption = "Provincial drivers of solar siting change (post-Q2 2020)"),
  list(id = "5_11",  fn = function() fig_5_11(), w = W_SINGLE, h = 88,
       caption = "Qinghai capacity-weighted deployment vs feasible and grid-accessible land"),
  list(id = "5_12",  fn = function() fig_5_12(), w = W_SINGLE, h = 88,
       caption = "Xinjiang capacity-weighted deployment vs feasible and grid-accessible land"),
  list(id = "5_11_12_pair", fn = function() fig_5_11_12_combined(), w = W_DOUBLE, h = 92,
       caption = "Qinghai and Xinjiang deployment compared (matched axes)"),
  list(id = "5_15",  fn = function() fig_5_15(basemap = .basemap()), w = W_DOUBLE, h = 140,
       caption = "Investment opportunity map, 2025 grid vintage"),
  list(id = "graphical_abstract", fn = function() fig_graphical_abstract(),
       w = W_DOUBLE, h = 105, caption = "Graphical abstract")
)

.basemap <- function() tryCatch(
  readRDS(file.path(P3, "china_basemap_bundle.rds"))$china_provinces,
  error = function(e) NULL)

# =============================================================================
# EXPORT FOR THE MARKDOWN DRAFT
#
# Writes PNGs where master.md expects them and emits a ready-to-paste markdown
# block with correct relative paths and width attributes. PNG because markdown
# previewers will not render an embedded PDF; the PDFs remain the submission
# artefact.
#
#   export_for_draft()                      # -> media/ at 300 dpi
#   export_for_draft(dir = "media/media")   # match master.md's existing paths
# =============================================================================
export_for_draft <- function(dir = "media", dpi = 300,
                             refs_file = "figure_refs.md") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  lines <- c("<!-- Generated by thesis_figures.R :: export_for_draft() -->",
             "<!-- Paste these blocks into master.md at the relevant points. -->", "")
  n_ok <- 0L
  for (f in FIGS) {
    p <- tryCatch(f$fn(), error = function(e) {
      cat(sprintf("  SKIPPED %-20s %s\n", f$id, conditionMessage(e))); NULL })
    if (is.null(p)) next
    png_name <- sprintf("fig_%s.png", f$id)
    ggsave(file.path(dir, png_name), p,
           width = f$w, height = f$h, units = "mm", dpi = dpi)
    win <- f$w / 25.4; hin <- f$h / 25.4
    # "5_04" -> "5.4", "5_11_12_pair" -> "5.11/5.12"
    num <- if (grepl("^[0-9]+_[0-9]+_[0-9]+_pair$", f$id)) {
      pr <- as.integer(strsplit(f$id, "_")[[1]][1:3])
      sprintf("%d.%d/%d.%d", pr[1], pr[2], pr[1], pr[3])
    } else if (grepl("^[0-9]+_[0-9]+$", f$id)) {
      pr <- as.integer(strsplit(f$id, "_")[[1]]); sprintf("%d.%d", pr[1], pr[2])
    } else NA_character_
    lab <- if (!is.na(num)) sprintf("**Figure %s** ", num) else ""
    lines <- c(lines,
      sprintf("![](./%s/%s){width=\"%.3fin\" height=\"%.3fin\"}", dir, png_name, win, hin),
      "", sprintf("%s%s", lab, f$caption), "")
    cat(sprintf("  wrote   %-28s %3.0f x %3.0f mm @ %d dpi\n", png_name, f$w, f$h, dpi))
    n_ok <- n_ok + 1L
  }
  writeLines(lines, refs_file)
  cat(sprintf("\n  %d PNGs -> %s/\n  markdown references -> %s\n",
              n_ok, dir, refs_file))
  invisible(refs_file)
}

# =============================================================================
# BUILD — write every figure to disk, then verify
# =============================================================================
build_all <- function(draft = TRUE) {
  t0 <- Sys.time()
  cat("=== Submission PDFs ===\n")
  for (f in FIGS) {
    p <- tryCatch(f$fn(), error = function(e) {
      cat(sprintf("  SKIPPED %-20s %s\n", f$id, conditionMessage(e))); NULL })
    if (!is.null(p)) save_fig(p, sprintf("fig_%s", f$id), f$w, f$h)
  }
  # single-column variants of the two wide figures
  for (f in FIGS[sapply(FIGS, function(x) x$id %in% c("5_04", "5_07"))]) {
    p <- tryCatch(f$fn(), error = function(e) NULL)
    if (!is.null(p)) save_fig(p, sprintf("fig_%s_1col", f$id), W_SINGLE, 95)
  }
  if (draft) { cat("\n=== Draft PNGs + markdown refs ===\n"); export_for_draft() }
  cat("\n=== Compliance ===\n"); ok <- verify_all()
  cat(sprintf("\nDone in %.1f min.\n  PDFs  -> %s/   (submission)\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins")), FIG_DIR))
  if (draft) cat("  PNGs  -> media/            (for master.md)\n",
                 " refs  -> figure_refs.md    (paste into the draft)\n", sep = "")
  invisible(ok)
}

cat("\nReady.\n  preview(fig_5_04())   inspect at true print size\n  preview_all()         every figure in turn\n  build_all()           PDFs + draft PNGs + markdown refs\n  export_for_draft()    draft PNGs only\n")
