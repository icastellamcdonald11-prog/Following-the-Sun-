source("config.R")  # repository paths; run scripts from the repository root
# 05_provincial_quadrant_map.R
# Two figures:
#   fig_05a — choropleth: dominant quadrant / Premium share by province
#   fig_05b — scatter: provincial capacity vs Premium share, labelled
#
# Data: plants_ngv_national.csv (from 03_national_analysis.R)
# Province geometries: china_basemap_bundle.rds (Project 3)

suppressPackageStartupMessages({
  library(sf); library(dplyr); library(ggplot2); library(scales); library(ggrepel)
})
sf_use_s2(FALSE)

P3  <- FTS_PROCESSED
P4  <- FTS_ROOT
OUT <- file.path(P4, "transmission")

# ── 1. Load data ───────────────────────────────────────────────────────────────
plants_all <- read.csv(file.path(OUT, "plants_ngv_national.csv"))

# Cumulative fleet — all years to Q2 2024
plants <- plants_all

cat(sprintf("Total plants (cumulative): %d\n", nrow(plants)))

basemap <- readRDS(file.path(P3, "china_basemap_bundle.rds"))
china_provinces  <- basemap$china_provinces
china_nation     <- basemap$china_nation
nine_dotted_line <- basemap$nine_dotted_line

# ── 2. Provincial quadrant summary ────────────────────────────────────────────
quad_order <- c("Premium", "Stranded", "Grid-convenient", "Dominated")
quad_cols  <- c(
  "Premium"         = "#2e9e56",
  "Stranded"        = "#e8936a",
  "Grid-convenient" = "#5b9bd5",
  "Dominated"       = "#c0392b"
)

prov_quad <- plants %>%
  mutate(quadrant = factor(quadrant, levels = quad_order)) %>%
  group_by(prov_display, quadrant) %>%
  summarise(cap_mw = sum(cap_mw, na.rm = TRUE), .groups = "drop") %>%
  group_by(prov_display) %>%
  mutate(
    cap_total   = sum(cap_mw),
    share       = cap_mw / cap_total * 100
  ) %>%
  ungroup()

# Wide: one row per province
prov_wide <- prov_quad %>%
  select(prov_display, quadrant, share, cap_total) %>%
  tidyr::pivot_wider(names_from = quadrant, values_from = share,
                     values_fill = 0) %>%
  mutate(
    cap_gw           = cap_total / 1e3,
    premium_gw       = cap_gw * Premium / 100,   # absolute GW in Premium quadrant
    stranded_gw      = cap_gw * Stranded / 100,
    dominated_gw     = cap_gw * Dominated / 100,
    # "above-median resource" = Premium + Stranded
    pct_above_median = Premium + Stranded,
    # dominant quadrant for choropleth
    dominant = quad_order[apply(
      cbind(Premium, Stranded, `Grid-convenient`, Dominated), 1, which.max
    )]
  )

cat("Provincial summary (top 10 by capacity):\n")
print(prov_wide %>% arrange(desc(cap_gw)) %>%
  select(prov_display, cap_gw, Premium, Stranded, `Grid-convenient`, Dominated) %>%
  head(10), digits = 1)

# ── 3. Join to province geometry ──────────────────────────────────────────────
# Standardise province name column
prov_sf <- china_provinces %>%
  rename(prov_display = name) %>%
  left_join(prov_wide, by = "prov_display") %>%
  st_make_valid()

cat(sprintf("\nProvinces matched: %d / %d\n",
            sum(!is.na(prov_sf$cap_gw)), nrow(prov_sf)))

# ── 4. Map projection ─────────────────────────────────────────────────────────
# Albers Equal Area for China (EPSG:102025)
albers <- "+proj=aea +lat_1=25 +lat_2=47 +lat_0=35 +lon_0=105 +datum=WGS84"

prov_sf_a   <- st_transform(prov_sf,          albers)
nation_a    <- st_transform(china_nation,      albers)
ndl_a       <- st_transform(nine_dotted_line,  albers)

# ── 5a. Choropleth: Premium share ─────────────────────────────────────────────
fig_choro <- ggplot() +
  geom_sf(data = prov_sf_a,
          aes(fill = Premium),
          colour = "white", linewidth = 0.25) +
  geom_sf(data = nation_a,
          fill = NA, colour = "grey30", linewidth = 0.4) +
  geom_sf(data = ndl_a,
          colour = "grey50", linewidth = 0.3, linetype = "dotted") +
  scale_fill_gradientn(
    colours  = c("#c0392b", "#e8c77a", "#2e9e56"),
    values   = scales::rescale(c(0, 30, 100)),
    limits   = c(0, 100),
    breaks   = c(0, 25, 50, 75, 100),
    labels   = paste0(c(0, 25, 50, 75, 100), "%"),
    name     = "Premium-sited\ncapacity (%)",
    na.value = "grey88",
    guide    = guide_colorbar(barwidth = 1, barheight = 8,
                               ticks.colour = "white")
  ) +
  labs(
    title    = "Share of Provincial UPV Capacity in Premium Quadrant",
    subtitle = "Premium = above-median PVOUT (>1,470 kWh/kWp/yr) × below-median transmission cost (<11.5 kWh/kWp/yr), cumulative to Q2 2024",
    caption  = paste0(
      "Post-shuangtan capacity (Q2 2020 – Q2 2024).  ",
      "Cumulative installed capacity to Q2 2024.  ",
      "Grid nodes: Fu et al. (2026), ≥220 kV substations, 2020 vintage.  ",
      "TX cost: Wang et al. (2022).  Plants: GRW Q2 2024 (Robinson et al., 2025)."
    )
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle   = element_text(size = 8.5, colour = "grey35", hjust = 0.5,
                                   margin = margin(b = 6)),
    plot.caption    = element_text(size = 7.5, colour = "grey45", hjust = 0),
    legend.position = c(0.95, 0.35),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin     = margin(8, 8, 8, 8)
  )

ggsave(file.path(OUT, "fig_05a_choropleth_premium_share.pdf"),
       fig_choro, width = 10, height = 7)
ggsave(file.path(OUT, "fig_05a_choropleth_premium_share.png"),
       fig_choro, width = 10, height = 7, dpi = 250)
cat("Saved fig_05a_choropleth_premium_share\n")

# ── 5b. Choropleth: dominant quadrant ─────────────────────────────────────────
prov_sf_a$dominant <- factor(prov_sf_a$dominant, levels = quad_order)

fig_dominant <- ggplot() +
  geom_sf(data = prov_sf_a,
          aes(fill = dominant),
          colour = "white", linewidth = 0.25) +
  geom_sf(data = nation_a,
          fill = NA, colour = "grey30", linewidth = 0.4) +
  geom_sf(data = ndl_a,
          colour = "grey60", linewidth = 0.3, linetype = "dotted") +
  scale_fill_manual(
    values   = quad_cols,
    na.value = "grey88",
    name     = "Dominant quadrant",
    drop     = FALSE
  ) +
  labs(
    title    = "Dominant Siting Quadrant by Province",
    subtitle = "Quadrant with the largest share of provincial installed capacity",
    caption  = paste0(
      "Cumulative installed capacity to Q2 2024.  ",
      "Grid nodes: Fu et al. (2026), ≥220 kV substations, 2020 vintage.  ",
      "TX cost: Wang et al. (2022).  Plants: GRW Q2 2024 (Robinson et al., 2025)."
    )
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle   = element_text(size = 8.5, colour = "grey35", hjust = 0.5,
                                   margin = margin(b = 6)),
    plot.caption    = element_text(size = 7.5, colour = "grey45", hjust = 0),
    legend.position = c(0.95, 0.35),
    legend.background = element_rect(fill = "white", colour = NA),
    plot.background = element_rect(fill = "white", colour = NA),
    plot.margin     = margin(8, 8, 8, 8)
  )

ggsave(file.path(OUT, "fig_05b_choropleth_dominant_quadrant.pdf"),
       fig_dominant, width = 10, height = 7)
ggsave(file.path(OUT, "fig_05b_choropleth_dominant_quadrant.png"),
       fig_dominant, width = 10, height = 7, dpi = 250)
cat("Saved fig_05b_choropleth_dominant_quadrant\n")

# ── 6. Scatter: provincial capacity vs Premium share ──────────────────────────
# Label all provinces with > 5 GW or unusual quadrant composition
label_df <- prov_wide %>%
  filter(cap_gw > 5 |
           Premium > 70 |
           Dominated > 50 |
           pct_above_median < 20)

fig_scatter <- ggplot(prov_wide, aes(x = Premium, y = cap_gw)) +
  # Quadrant background shading
  annotate("rect", xmin = 50, xmax = 100, ymin = -Inf, ymax = Inf,
           fill = "#2e9e56", alpha = 0.05) +
  annotate("rect", xmin = -Inf, xmax = 50, ymin = -Inf, ymax = Inf,
           fill = "#c0392b", alpha = 0.05) +
  geom_vline(xintercept = 50, linetype = "dashed", colour = "grey60",
             linewidth = 0.4) +
  # Points sized by capacity, coloured by above-median resource share
  geom_point(aes(colour = pct_above_median, size = cap_gw),
             alpha = 0.85, shape = 16) +
  geom_text_repel(data = label_df,
                  aes(label = prov_display),
                  size = 2.8, colour = "grey20",
                  segment.colour = "grey60", segment.size = 0.3,
                  box.padding = 0.4, max.overlaps = 25) +
  scale_colour_gradientn(
    colours = c("#c0392b", "#e8c77a", "#2e9e56"),
    values  = scales::rescale(c(0, 40, 100)),
    limits  = c(0, 100),
    name    = "Above-median\nresource (%)\n(Premium + Stranded)",
    guide   = guide_colorbar(barwidth = 1, barheight = 6)
  ) +
  scale_size_continuous(range = c(2, 12), guide = "none") +
  scale_x_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
  scale_y_continuous(labels = function(y) paste0(y, " GW")) +
  labs(
    title    = "Provincial UPV Capacity vs. Premium-Sited Share",
    subtitle = paste0(
      "Point size proportional to installed capacity.  ",
      "Vertical dashed line: 50% Premium threshold.  ",
      "Colour: combined above-median resource share (Premium + Stranded)."
    ),
    x       = "Share of provincial capacity in Premium quadrant (%)",
    y       = "Total installed capacity (GW)",
    caption = paste0(
      "Cumulative installed capacity to Q2 2024.  ",
      "Grid nodes: Fu et al. (2026), ≥220 kV substations, 2020 vintage.  ",
      "TX cost: Wang et al. (2022).  Plants: GRW Q2 2024 (Robinson et al., 2025)."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title      = element_text(face = "bold", size = 12),
    plot.subtitle   = element_text(size = 8.5, colour = "grey35"),
    plot.caption    = element_text(size = 7.5, colour = "grey45", hjust = 0),
    legend.position = "right",
    plot.background = element_rect(fill = "white", colour = NA),
    panel.grid.minor = element_blank(),
    plot.margin     = margin(8, 8, 8, 8)
  )

ggsave(file.path(OUT, "fig_05c_scatter_capacity_vs_premium.pdf"),
       fig_scatter, width = 11, height = 7)
ggsave(file.path(OUT, "fig_05c_scatter_capacity_vs_premium.png"),
       fig_scatter, width = 11, height = 7, dpi = 250)
cat("Saved fig_05c_scatter_capacity_vs_premium\n")

# ── 7. Print full province table ──────────────────────────────────────────────
cat("\n=== Full provincial quadrant table ===\n")
prov_wide %>%
  arrange(desc(cap_gw)) %>%
  select(prov_display, cap_gw, Premium, Stranded, `Grid-convenient`, Dominated,
         pct_above_median) %>%
  print(digits = 1, n = 35)

cat("\nDone.\n")
