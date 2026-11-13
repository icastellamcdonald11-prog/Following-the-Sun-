# Data

No third-party data are redistributed in this repository except the small derived tables listed under **Included**. Check each provider's terms before redistributing anything below.

## Included (derived)

| Path | Contents | Derived from |
|---|---|---|
| `fig_cache/plants_quad.rds` | 16,713 plants: centroid lon/lat, capacity (MW), PVOUT, province, build quarter, distance to nearest ≥220 kV substation (2020), transmission penalty, quadrant | GRW, Global Solar Atlas, Fu et al. |
| `fig_cache/grid_access.rds` | 97,457 cells of the 10 km equal-area grid: PVOUT, elevation, slope, feasibility, province, 2020 grid access | Global Solar Atlas, SRTM, Fu et al., provincial boundaries |
| `transmission/grid_nodes_2020.csv`, `grid_nodes_2025.csv` | Substations of 220 kV and above (3,450 and 6,539) | Fu et al. (2026) |
| `robustness_outputs/*.csv` | Saved outputs of the thesis robustness scripts | as above, plus WorldPop and Xin et al. |

Note: `robustness_outputs/table_R1_ks_test.csv` labels both Kolmogorov–Smirnov tests as "post < pre"; R's `alternative = "less"` tests the opposite direction. The direction-explicit statistics reported in the paper (Supplementary Table S14) are produced by `R/paper/si_distribution_tests.R`.

## Processed inputs (`data/processed/`, or set `FTS_PROCESSED`)

Needed by `run_paper.R` and the full pipeline. They are not in the repository; they will be added to the Zenodo record (https://doi.org/10.5281/zenodo.24186119) or are available from the corresponding author.

| File | Used by |
|---|---|
| `china_basemap_bundle.rds` (provinces including Taiwan, Hong Kong and Macao; national outline; South China Sea dashed line; inset settings) | Figs. 4–5, SI figures |
| `pvout_annual_10km.tif` (PVOUT, 10 km equal-area grid, EPSG:6933) | caches, Fig. 5 |
| `elevation_10km.tif`, `slope_deg_10km.tif`, `topography_allowed_mask_10km.tif` | feasibility filter, Fig. 4 contours |
| `grid10km_solar_pvout_table_correct_topo.rds` | land grid |
| `UPV_core_GRW_2024Q2_ge1MW.rds` (GRW polygons in China, ≥ 1 MW, with capacity) | caches, Fig. 4 footprints |
| `pvout_annual.tif`, `pvout_resource_by_province_FEASIBLE.csv`, `province_*.csv`, `misallocation_cost_*.csv` | `00_setup.R` |

## Derived rasters (`data/derived/`)

Built by `R/thesis/02_hires_1km_kde_figs.R` and `R/thesis/00_setup.R`, or taken from the archived release:
`pvout_1km_cea.tif`, `topo_mask_1km_cea.tif` (Fig. 4c,d), `pop_density_10km_cea.tif`, `pop_density_1km_cea.tif`, `cpvpd_processed.rds`.

## Raw sources (`data/raw/` and `transmission/`)

| Dataset | Source | Place at |
|---|---|---|
| Global Renewables Watch, solar, Q4 2017–Q2 2024 (Robinson et al. 2025) | arXiv:2503.14860 and the GRW data release | used to build `UPV_core_GRW_2024Q2_ge1MW.rds` (not scripted) |
| Global Solar Atlas v2.0 PVOUT (ESMAP 2019) | https://globalsolaratlas.info/download | used to build `pvout_annual*.tif` (not scripted) |
| SRTM elevation (Farr et al. 2007) | NASA SRTM | used to build the elevation and slope rasters (not scripted) |
| China power transmission network 2015/2020/2025 (Fu et al. 2026) | https://doi.org/10.6084/m9.figshare.31569094 | `transmission/network_2020/…/Substation_2020.shp`, `transmission/network_2025/…/Substation_2025.shp` |
| WorldPop China 2020, 1 km, persons km⁻² | https://doi.org/10.5258/SOTON/WP00674 | `data/raw/chn_pd_2020_1km_ASCII_XYZ.zip` |
| Provincial fixed-asset investment 1951–2020 (Xin et al. 2024) | https://zenodo.org/records/12706096 | `data/raw/1951-2020_provincial_fixed_assets.xlsx` |
| CPVPD-2024 (Yang et al. 2025) | https://doi.org/10.5281/zenodo.15618227 | `data/raw/CPVPD-2024_4326/` |
| Provincial boundaries (UC Berkeley Library 2022) | https://geodata.lib.berkeley.edu/catalog/berkeley-s7b118 | used to build `china_basemap_bundle.rds` (not scripted) |
| NEA 2018 provincial utilization hours | China Energy Portal (2019) | entered in `R/thesis/10_pvout_validation_casestudy.R` |
