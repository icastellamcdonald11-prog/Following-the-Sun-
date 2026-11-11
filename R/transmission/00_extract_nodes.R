source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# 00_extract_nodes.R
# Extract >=220kV substation centroids from the Figshare spatial dataset
# (Chen et al., 2026 — A Spatial Dataset of China's Power Transmission Network
#  Across 2015, 2020 and 2025, Figshare doi: 10.6084/m9.figshare.31569094)
#
# Outputs (to this folder):
#   grid_nodes_2020.csv  — 2020 substations >=220kV  (used for plant-level TX costs)
#   grid_nodes_2025.csv  — 2025 substations >=220kV  (used for investment opportunity map)
# =============================================================================

suppressPackageStartupMessages(library(sf))
suppressPackageStartupMessages(library(dplyr))

TX_OLD <- file.path(FTS_ROOT, "transmission")
OUT    <- file.path(FTS_ROOT, "transmission")

VOLT_MIN <- 220000   # 220 kV — minimum for utility-scale solar connection in China

extract_nodes <- function(shp_path, year, volt_min = VOLT_MIN) {
  cat(sprintf("Reading Substation_%d.shp ...\n", year))
  sub <- st_read(shp_path, quiet = TRUE)

  # voltage_1 is the parsed primary voltage (numeric, in Volts)
  sub_hi <- sub %>%
    filter(!is.na(voltage_1), voltage_1 >= volt_min)

  # Extract WGS84 coordinates
  coords <- sub_hi %>%
    st_transform("EPSG:4326") %>%
    st_coordinates() %>%
    as.data.frame() %>%
    rename(lng = X, lat = Y)

  nodes <- bind_cols(
    data.frame(
      id       = seq_len(nrow(sub_hi)),
      voltage  = sub_hi$voltage_1,
      status   = sub_hi$status
    ),
    coords
  ) %>%
    filter(is.finite(lat), is.finite(lng))

  cat(sprintf("  Year %d: %d substations >= %g kV\n",
              year, nrow(nodes), volt_min / 1000))

  # Province-level summary for QH and XJ
  xj <- nodes %>% filter(lng >= 73, lng <= 97, lat >= 34, lat <= 50)
  qh <- nodes %>% filter(lng >= 89, lng <= 103, lat >= 31, lat <= 40)
  cat(sprintf("    Xinjiang: %d  |  Qinghai: %d\n", nrow(xj), nrow(qh)))
  cat("    Voltage breakdown:\n")
  print(round(prop.table(table(nodes$voltage / 1000)) * 100, 1))

  nodes
}

nodes_2020 <- extract_nodes(
  file.path(TX_OLD, "network_2020/Transmission_network_dataset_2020/Substation_2020.shp"),
  2020
)

nodes_2025 <- extract_nodes(
  file.path(TX_OLD, "network_2025/Transmission_network_dataset_2025/Substation_2025.shp"),
  2025
)

write.csv(nodes_2020, file.path(OUT, "grid_nodes_2020.csv"), row.names = FALSE)
write.csv(nodes_2025, file.path(OUT, "grid_nodes_2025.csv"), row.names = FALSE)

cat("\nSaved grid_nodes_2020.csv and grid_nodes_2025.csv\n")
cat("Done.\n")
