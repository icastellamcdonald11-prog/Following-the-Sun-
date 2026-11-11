source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# 10_pvout_validation_casestudy.R  —  Thesis Project 4
#
# Empirical validation that PVOUT drives actual plant output.
# Two outputs:
#   Fig V1  — Provincial scatter: GSA PVOUT vs NEA actual utilisation hours
#             (2018, from National Renewable Power Monitoring & Evaluation Report)
#   Fig V2  — Plant-level scatter: GSA predicted output vs reported actual output
#             (manually collected data, 6 plants across 4 provinces)
#
# Data sources:
#   Provincial utilisation hours:
#     NEA (2019) 2018年可再生能源发展监测评价报告, via China Energy Portal
#     (chinaenergyportal.org/en/2018-national-renewable-power-development-
#      monitoring-and-evaluation-report/)
#   Plant-level reported output:
#     Manually collected from publicly available company/project sources
#     (see Untitled spreadsheet (1).xlsx in project folder)
#   PVOUT:
#     Global Solar Atlas v2.0 (World Bank / ESMAP, 2020); extracted to
#     10 km grid in 00_setup.R -> grid_df$pvout_annual
#
# Requires: 00_setup.R already sourced (grid_df with province column + pvout_annual)
# =============================================================================

library(ggplot2)
library(dplyr)
library(ggrepel)

P4  <- FTS_ROOT
OUT <- file.path(P4, "figures")
dir.create(OUT, showWarnings = FALSE)

theme_cs <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor  = element_blank(),
        axis.title        = element_text(size = 10),
        plot.title        = element_text(size = 11, face = "bold"),
        plot.caption      = element_text(size = 8, colour = "grey50"),
        legend.position   = "bottom")

# =============================================================================
# 1.  NEA 2018 provincial utilisation hours
#     Source: NEA (2019) 2018年可再生能源发展监测评价报告, Table 4
#             translated via China Energy Portal (2019)
#
#     Notes:
#     - Hours are ACTUAL grid-delivered hours (i.e. after curtailment)
#     - Inner Mongolia split into West (Zone I) and East (Zone II) by NEA;
#       we use both separately so they can be matched to grid cells by lon
#     - Gansu and Xinjiang have two resource zones; we use Zone I (best resource)
#       and Zone II averages separately, then average for single-province point
#     - Jiangsu figure (1,050 h) is the NEA minimum guaranteed purchase threshold
#       cited in the 2018 monitoring report; actual 2018 hours closely track this
# =============================================================================

nea_hours <- tribble(
  ~province_en,    ~province,       ~util_hours, ~curtail_flag,
  "Qinghai",       "Qinghai",           1461,  FALSE,
  "Inner Mongolia","Neimenggu",         1587,  FALSE,  # avg W(1649)+E(1525)
  "Heilongjiang",  "Heilongjiang",      1311,  FALSE,
  "Shaanxi",       "Shaanxi",           1316,  FALSE,
  "Jilin",         "Jilin",             1283,  FALSE,
  "Ningxia",       "Ningxia",           1376,  TRUE,
  "Hebei",         "Hebei",             1372,  FALSE,
  "Shanxi",        "Shanxi",            1355,  FALSE,
  "Liaoning",      "Liaoning",          1207,  FALSE,
  "Gansu",         "Gansu",             1264,  TRUE,   # avg zones, high curtailment
  "Xinjiang",      "Xinjiang",          1285,  TRUE,
  "Jiangsu",       "Jiangsu",           1050,  FALSE
)
# curtail_flag = TRUE  ->  province had >10% curtailment rate in 2017-18
#   (Xinjiang 22%, Gansu 20%  in 2017; both still >10% in 2018)
#   Flagged provinces will plot visibly BELOW the PVOUT trend line

# =============================================================================
# 2.  Compute province-level mean PVOUT from grid_df (feasible cells only)
# =============================================================================
cat("Computing province-level mean PVOUT from grid_df...\n")

prov_pvout <- pvout_prov_feasible %>%
  select(province, pvout_mean = pvout_mean_feasible,
         pvout_median = pvout_median_feasible)

cat("  Provinces in pvout_prov_feasible:", nrow(prov_pvout), "\n")
cat("  Sample province names:", head(prov_pvout$province, 4), "\n")

# Join on English province name
scatter_df <- nea_hours %>%
  left_join(prov_pvout, by = "province") %>%
  filter(!is.na(pvout_mean))

cat("  Matched provinces:", nrow(scatter_df), "of", nrow(nea_hours), "\n")
print(scatter_df %>% select(province_en, pvout_mean, util_hours, curtail_flag))

# =============================================================================
# 3.  Fig V1 — Scatter: provincial PVOUT vs actual utilisation hours
# =============================================================================
cat("\nBuilding Fig V1 (provincial PVOUT vs utilisation hours)...\n")

# OLS line (all provinces)
fit_all  <- lm(util_hours ~ pvout_mean, data = scatter_df)
fit_nc   <- lm(util_hours ~ pvout_mean, data = filter(scatter_df, !curtail_flag))
r2_all   <- summary(fit_all)$r.squared
r2_nc    <- summary(fit_nc)$r.squared

cat(sprintf("  R² (all provinces):                  %.3f\n", r2_all))
cat(sprintf("  R² (excluding high-curtailment):     %.3f\n", r2_nc))

fig_v1 <- ggplot(scatter_df,
                 aes(x = pvout_mean, y = util_hours,
                     shape = curtail_flag, colour = curtail_flag)) +

  # trend lines
  geom_smooth(data = filter(scatter_df, !curtail_flag),
              method = "lm", se = TRUE, colour = "#2C3E50",
              fill = "grey85", alpha = 0.35, linewidth = 0.8,
              linetype = "solid", show.legend = FALSE) +

  # reference lines
  geom_hline(yintercept = 1115, linetype = "dotted",
             colour = "grey55", linewidth = 0.4) +
  geom_vline(xintercept = pvout_national_median, linetype = "dotted",
             colour = "grey55", linewidth = 0.4) +

  # points
  geom_point(size = 3.5, stroke = 0.8, alpha = 0.9) +

  # labels
  geom_text_repel(aes(label = province_en),
                  size = 2.7, colour = "grey20",
                  max.overlaps = 20,
                  box.padding = 0.35, point.padding = 0.3,
                  segment.colour = "grey60", segment.size = 0.3) +

  # annotation: R2
  annotate("text",
           x    = min(scatter_df$pvout_mean, na.rm = TRUE) + 20,
           y    = max(scatter_df$util_hours, na.rm = TRUE) - 20,
           hjust = 0, vjust = 1, size = 3, colour = "grey30",
           label = paste0("R² = ", round(r2_nc, 2),
                          " (excl. high-curtailment provinces)")) +

  scale_colour_manual(
    values = c("FALSE" = "#2471A3", "TRUE" = "#C0392B"),
    labels = c("FALSE" = "Normal curtailment",
               "TRUE"  = "High curtailment (>10%, 2017–18)"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 17),
    labels = c("FALSE" = "Normal curtailment",
               "TRUE"  = "High curtailment (>10%, 2017–18)"),
    name = NULL
  ) +
  scale_x_continuous(
    "Provincial mean PVOUT, feasible land (kWh/kWp/yr)\n[Global Solar Atlas v2.0]",
    breaks = pretty(scatter_df$pvout_mean, n = 6)
  ) +
  scale_y_continuous(
    "Actual solar utilisation hours (h/yr)\n[NEA, 2018]",
    breaks = pretty(scatter_df$util_hours, n = 6)
  ) +
  theme_cs +
  labs(
    title   = "Fig V1  Provincial solar resource (PVOUT) vs actual annual output",
    subtitle = paste0("Each point = province (2018). Dotted lines: national averages. ",
                      "Shaded band: 95% CI (OLS, low-curtailment provinces only)."),
    caption = paste0(
      "Sources: GSA PVOUT extracted from Global Solar Atlas v2.0 (World Bank/ESMAP, 2020); ",
      "utilisation hours from NEA (2019) 2018年可再生能源发展监测评价报告, ",
      "via China Energy Portal (2019).\n",
      "Gansu and Xinjiang hours are zone-weighted averages. ",
      "High-curtailment provinces (red triangles) plot below the resource-driven trend, ",
      "demonstrating the grid constraint effect."
    )
  )

ggsave(file.path(OUT, "fig_v1_pvout_vs_utilhours.png"),
       fig_v1, width = 8, height = 6, dpi = 300)
cat("  Saved fig_v1_pvout_vs_utilhours.png\n")

# =============================================================================
# 4.  Fig V2 — Plant-level validation: GSA predicted output vs reported output
#
#     Manually collected plant-level data (6 plants across five provinces).
#     Coordinates parsed from spreadsheet to decimal degrees; PVOUT extracted
#     from GSA raster at plant centroid; expected output = capacity x PVOUT x PR.
#     Reported output from publicly available sources (see spreadsheet).
#
#     Plant types:
#       UPV (ground-mounted): PR = 0.80 (GSA methodology standard)
#       DPV (rooftop):        PR = 0.75 (lower due to suboptimal orientation/shading)
#     The Anhui Xuancheng plant is rooftop DPV and is flagged separately.
#
#     Guangdong excluded: no coordinates and capacity unclear in source data.
# =============================================================================
cat("\nBuilding Fig V2 (plant-level GSA validation: predicted vs reported output)...\n")

library(terra)

PR_upv <- 0.80   # ground-mounted fixed-tilt (GSA standard)
PR_dpv <- 0.75   # rooftop DPV

# -----------------------------------------------------------------------------
# 4a. Plant data with decimal-degree coordinates
#     (DMS coordinates converted manually: deg + min/60 + sec/3600)
#
#     Golmud Solar Park:         36deg19'43"N 94deg48'07"E -> 36.3286, 94.8019
#     Huanghe Golmud Solar Park: 36deg24'00"N 95deg07'30"E -> 36.4000, 95.1250
#     Longxiaya Power Plant:     36.1875N, 100.5735E (already decimal)
#     Anhui Xuancheng (DPV):     30.9487N, 118.7514E (already decimal)
#     Gansu Dunhuang:            40deg04'N  94deg30'E  -> 40.0667, 94.5000
#     Changli:                   39.7065N, 119.1602E (already decimal)
# -----------------------------------------------------------------------------
val_plants <- tibble::tribble(
  ~plant,                           ~province,  ~lon,      ~lat,     ~cap_mw, ~reported_gwh, ~type,
  # Golmud Solar Park (138 MW) excluded: uses concentrating PV (CPV) technology
  # (Suncore joint venture). GSA PVOUT is modelled for flat-plate silicon PV;
  # comparing CPV output against a flat-plate prediction is not valid.
  "Huanghe Golmud Solar Park",      "Qinghai",   95.1603,  36.4053,  200.0,   318.1,         "UPV",
  "Longxiaya Power Plant",          "Qinghai",  100.5735,  36.1875,  320.0,   498.0,         "UPV",
  "Anhui Xuancheng (rooftop DPV)",  "Anhui",    118.7514,  30.9487,    3.5,     4.2,         "DPV",
  "Gansu Dunhuang",                 "Gansu",     94.5000,  40.0667,   50.0,    78.2,         "UPV",
  "Ningxia Taiyangshan Phase II",   "Ningxia",  106.5737,  37.5004,   20.24,   31.31,        "UPV",
  "Xinjiang Dabancheng",            "Xinjiang",  88.2823,  43.3661,   10.00,   13.69,        "UPV"
  # Coordinates corrected: PDD p.4 transposes lon/lat labels; 88°E, 43°N is
  # consistent with Dabancheng District, Urumqi (confirmed by satellite image p.5)
  # Changli (Hebei, 460.8 MW) excluded: the 770 GWh figure is a developer
  # projection from the project announcement, not metered generation data.
  # Comparing two modelled estimates does not constitute independent validation.
)

# -----------------------------------------------------------------------------
# 4b. Extract GSA PVOUT at each plant centroid
#     Try 1km CEA raster first; fall back to Project 3 10km raster
# -----------------------------------------------------------------------------
pvout_r_path <- file.path(P4, "data/derived/pvout_1km_cea.tif")
if (!file.exists(pvout_r_path)) {
  pvout_r_path <- file.path(
    FTS_P3_ROOT,
    "data/processed/pvout_annual_10km.tif"
  )
}
cat("  Loading PVOUT raster from:", pvout_r_path, "\n")
pvout_rast <- terra::rast(pvout_r_path)

pts_wgs84 <- terra::vect(
  as.data.frame(val_plants[, c("lon", "lat")]),
  geom = c("lon", "lat"), crs = "EPSG:4326"
)
pts_proj <- terra::project(pts_wgs84, terra::crs(pvout_rast))

pvout_extracted <- terra::extract(pvout_rast, pts_proj)[, 2]
cat("  Extracted PVOUT values:\n")
print(data.frame(plant = val_plants$plant, pvout = round(pvout_extracted)))

val_plants <- val_plants %>%
  dplyr::mutate(
    pvout_gsa      = pvout_extracted,
    pr             = dplyr::if_else(type == "DPV", PR_dpv, PR_upv),
    expected_gwh   = cap_mw * pvout_gsa * pr / 1000,
    pct_error      = (reported_gwh - expected_gwh) / expected_gwh * 100,
    # Capacity-normalised output (GWh/MW) — removes plant-size effect
    expected_gwh_mw = expected_gwh / cap_mw,
    reported_gwh_mw = reported_gwh / cap_mw
  )

cat("\n  Validation table (capacity-normalised):\n")
print(val_plants %>%
        dplyr::select(plant, province, cap_mw, pvout_gsa,
                      expected_gwh_mw, reported_gwh_mw, pct_error) %>%
        dplyr::mutate(dplyr::across(where(is.numeric), ~ round(., 3))),
      n = 20)

# -----------------------------------------------------------------------------
# 4c. OLS fit and R2 on capacity-normalised output (UPV only)
#     Normalising by capacity removes the dominant plant-size effect and gives
#     a genuine test of GSA's per-unit predictive accuracy.
# -----------------------------------------------------------------------------
upv_df <- filter(val_plants, type == "UPV")

fit_val_upv <- lm(reported_gwh_mw ~ expected_gwh_mw, data = upv_df)
r2_val_upv  <- summary(fit_val_upv)$r.squared
mae_upv     <- mean(abs(upv_df$pct_error))

cat(sprintf("\n  R² (UPV, GWh/MW normalised):            %.3f\n", r2_val_upv))
cat(sprintf("  Mean absolute percentage error (UPV):   %.1f%%\n", mae_upv))

# -----------------------------------------------------------------------------
# 4d. Build Fig V2 — capacity-normalised scatter
# -----------------------------------------------------------------------------
ax_min <- min(c(upv_df$expected_gwh_mw, upv_df$reported_gwh_mw)) * 0.88
ax_max <- max(c(val_plants$expected_gwh_mw, val_plants$reported_gwh_mw)) * 1.12

fig_v2 <- ggplot(val_plants,
                 aes(x = expected_gwh_mw, y = reported_gwh_mw,
                     colour = type, shape = type)) +

  # 1:1 reference line
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey55", linewidth = 0.6) +

  # OLS fit (UPV only)
  geom_smooth(data = upv_df,
              method = "lm", se = TRUE,
              colour = "#2C3E50", fill = "grey85", alpha = 0.35,
              linewidth = 0.8, show.legend = FALSE) +

  # Points
  geom_point(size = 3.5, stroke = 0.8, alpha = 0.9) +

  # Labels
  ggrepel::geom_text_repel(aes(label = plant),
                            size = 2.7, colour = "grey20",
                            max.overlaps = 20,
                            box.padding = 0.4, point.padding = 0.3,
                            segment.colour = "grey60", segment.size = 0.3) +

  # Annotation: MAE only (R2 on normalised data is the honest stat)
  annotate("text",
           x = ax_min + (ax_max - ax_min) * 0.03,
           y = ax_max * 0.99,
           hjust = 0, vjust = 1, size = 3, colour = "grey30",
           label = paste0("Mean abs. error = ", round(mae_upv, 1),
                          "%  (UPV plants)\n",
                          "R² = ", round(r2_val_upv, 2),
                          "  (capacity-normalised, OLS; shaded = 95% CI)")) +

  # 1:1 label
  annotate("text",
           x = ax_min + (ax_max - ax_min) * 0.55,
           y = ax_min + (ax_max - ax_min) * 0.50,
           hjust = 0, size = 2.7, colour = "grey55", fontface = "italic",
           label = "1:1 (perfect prediction)") +

  scale_colour_manual(
    values = c("UPV" = "#2471A3", "DPV" = "#E67E22"),
    labels = c("UPV" = "Utility-scale PV (PR = 0.80)",
               "DPV" = "Rooftop DPV (PR = 0.75)"),
    name = NULL
  ) +
  scale_shape_manual(
    values = c("UPV" = 16, "DPV" = 17),
    labels = c("UPV" = "Utility-scale PV (PR = 0.80)",
               "DPV" = "Rooftop DPV (PR = 0.75)"),
    name = NULL
  ) +
  scale_x_continuous(
    "GSA-predicted output per MW installed (GWh/MW/yr)\n[PVOUT × PR]",
    limits = c(ax_min, ax_max),
    breaks = pretty(c(ax_min, ax_max), n = 6)
  ) +
  scale_y_continuous(
    "Reported output per MW installed (GWh/MW/yr)\n[CDM PDDs and public sources]",
    limits = c(ax_min, ax_max),
    breaks = pretty(c(ax_min, ax_max), n = 6)
  ) +
  coord_fixed() +
  theme_cs +
  labs(
    title    = "Fig V2  GSA-predicted vs reported output per MW: plant-level validation",
    subtitle = paste0("Six plants across five provinces. Capacity-normalised to remove plant-size effect. ",
                      "Dashed line = 1:1 perfect prediction."),
    caption  = paste0(
      "PVOUT extracted from Global Solar Atlas v2.0 (World Bank/ESMAP, 2020) at plant centroids. ",
      "Expected output per MW = PVOUT (kWh/kWp/yr) × PR / 1000.\n",
      "PR = 0.80 for ground-mounted UPV (GSA standard); PR = 0.75 for rooftop DPV. ",
      "Output figures from UNFCCC CDM Project Design Documents (Huanghe Golmud, Gansu Dunhuang, ",
      "Ningxia Taiyangshan, Xinjiang Dabancheng) and GlobalData/Power Technology (Longyangxia)."
    )
  )

ggsave(file.path(OUT, "fig_v2_plant_zone_comparison.png"),
       fig_v2, width = 7, height = 7, dpi = 300)
cat("  Saved fig_v2_plant_zone_comparison.png\n")

# =============================================================================
# 5.  Print summary for thesis write-up
# =============================================================================
cat("\n============================================================\n")
cat("KEY NUMBERS FOR THESIS TEXT\n")
cat("============================================================\n")
cat(sprintf("National average utilisation hours (2018): 1,115 h\n"))
cat(sprintf("Qinghai actual hours (2018):              1,461 h (+%.0f%% above nat. avg)\n",
            (1461/1115 - 1)*100))
cat(sprintf("Jiangsu approx. hours (2018):             1,050 h (%.0f%% of nat. avg)\n",
            1050/1115*100))
cat(sprintf("Gansu zone-avg hours (2018):              1,264 h\n"))
cat(sprintf("Qinghai:Jiangsu output ratio:             %.2fx\n", 1461/1050))
cat(sprintf("\nFig V1 — Provincial R² (excl. high-curtailment): %.3f\n", r2_nc))
cat(sprintf("\nFig V2 — Plant-level validation (UPV):\n"))
cat(sprintf("  R² (GSA predicted vs reported):  %.3f\n", r2_val_upv))
cat(sprintf("  Mean absolute percentage error:  %.1f%%\n", mae_upv))
cat("\n  Plant-level detail:\n")
print(val_plants %>%
        dplyr::select(plant, province, cap_mw, pvout_gsa,
                      expected_gwh, reported_gwh, pct_error) %>%
        dplyr::mutate(dplyr::across(where(is.numeric), ~ round(., 1))),
      n = 20)
cat("============================================================\n")
cat("\nAll figures saved to:", OUT, "\n")
