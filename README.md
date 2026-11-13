# Following the Sun? Code and derived data

Code for **Castella McDonald, I. & Pan, Y. *Following the Sun? Utility-scale solar siting efficiency and the resource frontier in China, 2017–2024*** (manuscript, 2026). <!-- TODO: journal reference and DOI on publication -->

The analysis measures the resource quality of China's utility-scale photovoltaic fleet (16,713 plants, 356.3 GW, Q4 2017–Q2 2024). It combines a satellite plant inventory with photovoltaic output potential (PVOUT) and the locations of high-voltage substations, and classifies each plant by resource quality and transmission cost.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24186119.svg)](https://doi.org/10.5281/zenodo.24186119)

Archived version: https://doi.org/10.5281/zenodo.24186119 (v1.0.0)

## Quick start: reproduce the paper's figures and recomputed numbers

```bash
git clone https://github.com/icastellamcdonald11-prog/Following-the-Sun-.git following-the-sun-code
cd following-the-sun-code
# place the processed inputs (data/README.md) in data/processed/ and data/derived/
Rscript run_paper.R
```

This runs, in order:

| Script | Produces | Paper item |
|---|---|---|
| `R/paper/make_figures.R` | `outputs/figures/fig1.pdf`–`fig5.pdf` | Figs. 1–5 |
| `R/paper/foregone_generation.R` | `outputs/analysis/foregone_generation*.csv` | Methods eq. 15; Supplementary Note 9, Table S10 |
| `R/paper/deterioration_shares.R` | `outputs/analysis/deterioration_*.csv` | Methods eq. 12; Supplementary Note 10, Table S11 |
| `R/paper/si_distribution_tests.R` | `outputs/analysis/si_distribution_tests.csv` | Supplementary Note 17, Table S14 |
| `R/paper/si_robustness_figs.R` | `outputs/figures/si/*.pdf` | Supplementary Figs. S20–S22 |

They read the derived plant and land tables in `fig_cache/`, the substation lists in `transmission/` and the saved test results in `robustness_outputs/`, all included here.

## Repository layout

```
config.R                 paths (FTS_ROOT, FTS_PROCESSED); sourced by every script
run_paper.R              runs the five paper scripts above
R/paper/                 scripts written for the manuscript (figures, recomputed numbers)
R/thesis/                analysis scripts from the underlying thesis
R/transmission/          substation extraction and transmission-cost scripts
fig_cache/               derived: plant table (plants_quad.rds), 10 km land grid (grid_access.rds)
transmission/            derived: substations >= 220 kV, 2020 and 2025 (from Fu et al. 2026)
robustness_outputs/      derived: saved outputs of the robustness scripts (CSV)
data/                    inputs not included; see data/README.md
outputs/                 written by the scripts
```

## Full pipeline (from processed inputs)

Run from the repository root. Each step writes files that later steps read.

1. `R/transmission/00_extract_nodes.R`: substations of 220 kV and above from the Fu et al. (2026) transmission dataset → `transmission/grid_nodes_2020.csv`, `grid_nodes_2025.csv`.
2. `R/thesis/00_setup.R`: loads the processed inputs and builds `deploy_df` (plants) and `grid_df` (10 km land grid), including population density. Source it in the R session before steps 4–5.
3. `R/thesis/build_caches.R`: builds `fig_cache/plants_quad.rds` (centroid PVOUT, capacity, province, distance to the nearest 2020 substation, transmission penalty, quadrant) and `fig_cache/grid_access.rds` (land grid with 2020 grid access).
4. `R/thesis/02_hires_1km_kde_figs.R`: builds the 1 km PVOUT and terrain rasters in `data/derived/` used by the Fig. 4 maps.
5. Statistics reported in the paper, each run after `00_setup.R`:
   - `11_siting_efficiency_regression.R`: within-province rank, fixed-effects model (Methods eq. 13), speed regression (eq. 14).
   - `13_model_spec_robustness.R`: alternative standard errors (`robustness_outputs/se_robustness_table.csv`).
   - `12_variance_decomposition_figs.R`: quarterly dispersion and quartile shares (Supplementary Note 17).
   - `14_land_opportunity_cost.R`, `14b_land_opportunity_cost_FAI.R`: land-cost tests (Supplementary Note 18).
   - `robustness_checks.R`: leave-one-province-out and threshold sensitivity (Supplementary Notes 6 and 8).
   - `10_pvout_validation_casestudy.R`: PVOUT against utilization hours (Supplementary Note 5).
   - `08_cpvpd_validation.R`: comparison with CPVPD-2024 (Supplementary Note 3).
   - `03_subgroup_kde_figs.R`, `06_nea_capacity_figs.R`, `R/transmission/01_*`, `05_*`: Supplementary figures carried over from the thesis.
6. `run_paper.R`.

**Not scripted:** the preparation of the processed inputs from the raw downloads (clipping the GRW polygons to China and applying the 1 MW filter and capacity schedule; aggregating PVOUT, elevation and slope to the 10 km grid; the basemap bundle) was carried out interactively and is not included. The processed inputs themselves are deposited with the archived release (see `data/README.md`).

## Key definitions (as in the paper's Methods)

- **Capacity**: footprint area × power density of 51.50 MW km⁻² before 2019, rising by 4.33 MW km⁻² a year to 65.90 MW km⁻² from 2022.
- **Plant PVOUT**: value of the 10 km cell containing the polygon centroid.
- **Feasible land**: elevation ≤ 5,000 m and slope ≤ 50°; national median PVOUT 1,470 kWh kWp⁻¹ yr⁻¹.
- **Transmission penalty**: (0.0002 CNY kWh⁻¹ km⁻¹ × distance, capped at 600 km) / 0.37 CNY kWh⁻¹ × PVOUT; median across plants 11.4 kWh kWp⁻¹ yr⁻¹.
- **Post-2020**: commissioned in Q2 2020 or later (quadrant shares of post-2020 additions count from Q3 2020).

## Requirements

Tested with R 4.5.2: sf 1.0.23, terra 1.8.86, dplyr 1.1.4, ggplot2 4.0.1, patchwork 1.3.2, ggrepel 0.9.6, scales 1.4.0, fixest 0.14.0, data.table 1.18.2.1, tidyr 1.3.2.

```r
install.packages(c("sf", "terra", "dplyr", "ggplot2", "patchwork", "ggrepel",
                   "scales", "fixest", "data.table", "tidyr"))
```

Figures use R's base `pdf()` device with Helvetica. Fonts are referenced, not embedded; run Ghostscript `embedFonts` on the PDFs if a publisher requires embedding.

## Data sources

See `data/README.md` for every input, its source, licence and expected location.

## Licence

Code: MIT (see `LICENSE`). Derived data in `fig_cache/`, `transmission/` and `robustness_outputs/` are derived from the datasets listed in `data/README.md` and remain subject to their licences.

## Citation

See `CITATION.cff`.
