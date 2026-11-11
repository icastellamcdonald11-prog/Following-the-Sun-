source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# 00_setup.R  —  Thesis Project 4: Solar Efficiency Evaluation
# Loads all processed datasets from Project 3 and rebuilds the core
# working objects (deploy_df, grid_df, basemap) from scratch.
# =============================================================================

library(sf)
library(terra)
library(dplyr)
library(lubridate)
library(ggplot2)

# Disable S2 spherical geometry engine — avoids strict validation errors on
# polygons with degenerate edges (duplicate vertices). Falls back to GEOS.
sf::sf_use_s2(FALSE)

P3        <- FTS_PROCESSED
P3_ckpt   <- FTS_PROCESSED

# =============================================================================
# 1.  China basemap  (load first — needed for province spatial join below)
# =============================================================================
cat("Loading China basemap bundle...\n")
basemap          <- readRDS(file.path(P3, "china_basemap_bundle.rds"))
china_nation     <- basemap$china_nation
china_provinces  <- basemap$china_provinces
nine_dotted_line <- basemap$nine_dotted_line
cat("  Province features:", nrow(china_provinces), "\n")

# =============================================================================
# 2.  GRW solar plant polygons  (UPV ≥ 1 MW, China, through Q2 2024)
#     Columns confirmed in data:
#       build_year, build_quarter, year_quarter, area_km2,
#       power_density_mw_km2, capacity_mw_used  (capacity already computed)
# =============================================================================
cat("Loading GRW UPV plant data...\n")
upv_sf <- readRDS(file.path(P3, "UPV_core_GRW_2024Q2_ge1MW.rds"))
cat("  Loaded:", nrow(upv_sf), "plant polygons\n")

# Ensure consistent CRS with basemap
if (!st_crs(upv_sf) == st_crs(china_provinces)) {
  upv_sf <- st_transform(upv_sf, st_crs(china_provinces))
}

# =============================================================================
# 3.  Assign province via spatial join
# =============================================================================
cat("Assigning provinces (spatial join)...\n")

# Extract the province name column from china_provinces
cat("  All columns in china_provinces:", paste(names(china_provinces), collapse = ", "), "\n")

# First try known candidates; if none match, fall back to the first
# non-geometry character/factor column
candidates <- c("NAME_1", "name_1", "name", "prov", "province", "VARNAME_1", "NL_NAME_1")
prov_name_col <- Filter(function(x) x %in% names(china_provinces), candidates)[1]

if (is.na(prov_name_col)) {
  non_geom <- names(china_provinces)[names(china_provinces) != attr(china_provinces, "sf_column")]
  char_cols <- non_geom[sapply(non_geom, function(cn)
                    is.character(china_provinces[[cn]]) | is.factor(china_provinces[[cn]]))]
  prov_name_col <- char_cols[1]
  cat("  No known column matched — using first character column:", prov_name_col, "\n")
}
stopifnot("Could not identify a province name column in china_provinces" = !is.na(prov_name_col))
cat("  Using province name column:", prov_name_col, "\n")

prov_sf <- china_provinces[, c(prov_name_col, "geometry")]
names(prov_sf)[1] <- "prov"

# Use centroids for join (faster, avoids boundary ambiguity)
plant_centroids <- st_centroid(upv_sf)
joined          <- st_join(plant_centroids, prov_sf, join = st_within, left = TRUE)

upv_sf$prov <- joined$prov

# Plants that fall outside province polygons (coastal/boundary edge effects):
# assign to nearest province
n_unmatched <- sum(is.na(upv_sf$prov))
if (n_unmatched > 0) {
  cat("  Assigning", n_unmatched, "unmatched plants to nearest province...\n")
  unmatched_idx  <- which(is.na(upv_sf$prov))
  nearest_idx    <- st_nearest_feature(plant_centroids[unmatched_idx, ], prov_sf)
  upv_sf$prov[unmatched_idx] <- prov_sf$prov[nearest_idx]
}
cat("  Province assignment complete. NA remaining:", sum(is.na(upv_sf$prov)), "\n")

# =============================================================================
# 4.  Extract PVOUT at each plant location from the 10 km raster
# =============================================================================
cat("Extracting PVOUT from 10 km raster...\n")
pvout_tif <- file.path(P3, "pvout_annual_10km.tif")
if (!file.exists(pvout_tif)) {
  pvout_tif <- file.path(P3_ckpt, "pvout_annual_10km.tif")
  cat("  (primary path missing — using checkpoint copy)\n")
}
stopifnot("pvout_annual_10km.tif not found in processed/ or checkpoint/" = file.exists(pvout_tif))
pvout_r   <- rast(pvout_tif)
pvout_ext <- terra::extract(pvout_r, vect(plant_centroids), method = "simple")
# pvout_ext: two columns — ID and the raster layer value
pvout_col_name <- names(pvout_ext)[2]
upv_sf$pvout   <- pvout_ext[[pvout_col_name]]
cat("  PVOUT extracted. NAs:", sum(is.na(upv_sf$pvout)), "/", nrow(upv_sf), "\n")

# =============================================================================
# 5.  Build deploy_df  —  flat data.frame, one row per plant
#
#     Capacity-weighting schedule (thesis §5.2):
#       build_year < 2019   →  51.50 MW/km²
#       build_year == 2019  →  55.83 MW/km²
#       build_year == 2020  →  60.16 MW/km²
#       build_year == 2021  →  64.49 MW/km²
#       build_year >= 2022  →  65.90 MW/km²   (constant)
#
#     NOTE: upv_sf already carries capacity_mw_used (pre-computed with this
#     schedule). We keep it as cap_mw for consistency with downstream code.
# =============================================================================
cat("Building deploy_df...\n")
deploy_df <- upv_sf %>%
  st_drop_geometry() %>%
  transmute(
    build_year    = as.integer(build_year),
    build_quarter = as.integer(build_quarter),
    year_quarter  = as.character(year_quarter),
    area_km2      = as.numeric(area_km2),
    pd_mw_km2     = as.numeric(power_density_mw_km2),
    cap_mw        = as.numeric(capacity_mw_used),
    pvout         = as.numeric(pvout),
    prov          = as.character(prov)
  ) %>%
  filter(
    is.finite(pvout),
    is.finite(cap_mw),
    cap_mw > 0,
    !is.na(prov)
  )

cat("  deploy_df rows:", nrow(deploy_df), "\n")
cat("  Total capacity (GW):", round(sum(deploy_df$cap_mw) / 1e3, 1), "\n")
cat("  Quarters covered:", paste(range(deploy_df$year_quarter), collapse = " – "), "\n")
cat("  Provinces:", length(unique(deploy_df$prov)), "\n")

# =============================================================================
# 6.  Solar resource matrix  (10 km grid, PVOUT + topography)
#     Columns: cell_id, x, y, pvout_annual, cell_area_km2,
#              is_allowed_topography, elevation_m, slope_deg
# =============================================================================
cat("Loading 10 km solar resource grid...\n")
grid_df <- readRDS(file.path(P3, "grid10km_solar_pvout_table_correct_topo.rds"))

grid_df <- grid_df %>%
  mutate(feasible = is_allowed_topography == 1 | is_allowed_topography == TRUE)

cat("  Grid cells total:", nrow(grid_df), "\n")
cat("  Feasible cells:  ", sum(grid_df$feasible, na.rm = TRUE), "\n")

# National PVOUT thresholds (feasible land only)
feasible_pvout        <- grid_df$pvout_annual[grid_df$feasible & is.finite(grid_df$pvout_annual)]
pvout_national_median <- median(feasible_pvout)
pvout_q75_national    <- quantile(feasible_pvout, 0.75)
pvout_q80_national    <- quantile(feasible_pvout, 0.80)
pvout_q90_national    <- quantile(feasible_pvout, 0.90)

cat(sprintf(
  "  National feasible PVOUT — median: %.1f  Q75: %.1f  Q80: %.1f  Q90: %.1f\n",
  pvout_national_median, pvout_q75_national, pvout_q80_national, pvout_q90_national
))

# =============================================================================
# 7.  Province-level resource summaries and pre-computed output tables
# =============================================================================
cat("Loading province-level summary tables...\n")

pvout_prov_feasible <- read.csv(
  file.path(P3, "pvout_resource_by_province_FEASIBLE.csv"),
  stringsAsFactors = FALSE
)

prov_gap <- read.csv(
  file.path(P3, "province_capacity_weighted_pvout_gap_DIRECT.csv"),
  stringsAsFactors = FALSE
)

misalloc_china_wide <- read.csv(
  file.path(P3, "misallocation_cost_top25_CHINA_WIDE.csv"),
  stringsAsFactors = FALSE
)

misalloc_within_prov <- read.csv(
  file.path(P3, "misallocation_cost_top25_within_province.csv"),
  stringsAsFactors = FALSE
)

prov_top25_shares <- read.csv(
  file.path(P3, "province_top25_feasible_pvout_shares.csv"),
  stringsAsFactors = FALSE
)

# =============================================================================
# 8.  Population density raster  (WorldPop 2020, 1 km → 10 km CEA)
#     Source: chn_pd_2020_1km_ASCII_XYZ.zip  (X=lon, Y=lat, Z=persons/km²)
#     Steps:  unzip → read XYZ → build 1km WGS84 raster → project to CEA
#             → resample (sum) to 10km → join to grid_df as pop_density
#     Result: grid_df$pop_density  (persons/km² at 10km CEA resolution)
# =============================================================================
P4 <- FTS_ROOT
pop_cache_path <- file.path(P4, "data/derived/pop_density_10km_cea.tif")

if (file.exists(pop_cache_path)) {
  cat("Loading cached population raster...\n")
  pop_10km <- rast(pop_cache_path)
} else {
  cat("Processing population density data (762 MB CSV — takes a few minutes)...\n")

  pop_zip <- file.path(P4, "data/raw/chn_pd_2020_1km_ASCII_XYZ.zip")
  stopifnot(file.exists(pop_zip))

  # Unzip to temp dir — take the largest file regardless of extension
  pop_dir <- tempfile()
  dir.create(pop_dir)
  unzip(pop_zip, exdir = pop_dir)
  all_unzipped <- list.files(pop_dir, full.names = TRUE, recursive = TRUE)
  pop_csv <- all_unzipped[which.max(file.size(all_unzipped))]
  cat("  Largest unzipped file:", pop_csv, "\n")
  cat("  File size:", round(file.size(pop_csv) / 1e6, 1), "MB\n")

  # Read XYZ — let fread auto-detect structure, then normalise column names
  if (requireNamespace("data.table", quietly = TRUE)) {
    pop_df <- data.table::fread(pop_csv, showProgress = TRUE)
    pop_df <- as.data.frame(pop_df)
  } else {
    cat("  (data.table not available; using read.csv — may be slow)\n")
    pop_df <- read.csv(pop_csv, stringsAsFactors = FALSE)
  }
  cat("  Columns detected:", paste(names(pop_df), collapse = ", "), "\n")
  cat("  Rows read:", nrow(pop_df), "\n")
  stopifnot("Need at least 3 columns (X, Y, Z)" = ncol(pop_df) >= 3)

  # Rename first three columns to x, y, z regardless of original names
  names(pop_df)[1:3] <- c("x", "y", "z")

  # Drop NODATA / negative rows and keep only the three coordinate columns
  pop_df <- pop_df[is.finite(pop_df$z) & pop_df$z >= 0, c("x", "y", "z")]
  cat("  Valid XYZ rows:", nrow(pop_df), "\n")

  # Build 1 km WGS84 raster from XYZ
  pop_r_wgs84 <- rast(pop_df[, c("x", "y", "z")], type = "xyz", crs = "EPSG:4326")

  # Project to CEA (same projection as pvout_r)
  cat("  Reprojecting to CEA...\n")
  pop_r_cea <- project(pop_r_wgs84, crs(pvout_r), method = "bilinear",
                        res = 1000)          # keep 1 km in CEA

  # Resample to match pvout_r exactly (10 km grid, same extent/origin)
  cat("  Resampling to 10 km grid...\n")
  pop_10km <- resample(pop_r_cea, pvout_r, method = "sum")

  # Cache to disk for future sessions
  writeRaster(pop_10km, pop_cache_path, overwrite = TRUE)
  cat("  Cached to:", pop_cache_path, "\n")

  # Clean up temp
  unlink(pop_dir, recursive = TRUE)
}

# Join population density to grid_df via point extraction
cat("Extracting population density to grid_df...\n")
grid_pts <- vect(
  data.frame(x = grid_df$x, y = grid_df$y),
  geom = c("x", "y"),
  crs  = crs(pvout_r)
)
pop_extracted       <- terra::extract(pop_10km, grid_pts)
grid_df$pop_density <- pop_extracted[[2]]

cat(sprintf(
  "  pop_density attached. Range: %.2f – %.2f  (NAs: %d)\n",
  min(grid_df$pop_density, na.rm = TRUE),
  max(grid_df$pop_density, na.rm = TRUE),
  sum(is.na(grid_df$pop_density))
))

# =============================================================================
# 9.  Environment summary
# =============================================================================
cat("\n========================================\n")
cat("Environment ready. Key objects:\n\n")
cat(sprintf("  upv_sf              sf polygons, UPV plants      (%d rows)\n",  nrow(upv_sf)))
cat(sprintf("  deploy_df           flat plant data.frame        (%d rows)\n",  nrow(deploy_df)))
cat(sprintf("  grid_df             10 km resource grid          (%d rows)\n",  nrow(grid_df)))
cat(sprintf("  pvout_prov_feasible province PVOUT summary       (%d rows)\n",  nrow(pvout_prov_feasible)))
cat(sprintf("  prov_gap            province PVOUT gap table     (%d rows)\n",  nrow(prov_gap)))
cat(sprintf("  pop_density         population on grid_df (cells w/ data: %d)\n",
            sum(!is.na(grid_df$pop_density))))
cat("\n  Scalars:\n")
cat(sprintf("    pvout_national_median  = %.2f\n", pvout_national_median))
cat(sprintf("    pvout_q75_national     = %.2f  (top-quartile threshold)\n", pvout_q75_national))
cat(sprintf("    pvout_q80_national     = %.2f  (top-quintile threshold)\n", pvout_q80_national))
cat(sprintf("    pvout_q90_national     = %.2f  (top-decile threshold)\n",  pvout_q90_national))
cat("\n  Map objects: china_nation, china_provinces, nine_dotted_line\n")
cat("========================================\n")
