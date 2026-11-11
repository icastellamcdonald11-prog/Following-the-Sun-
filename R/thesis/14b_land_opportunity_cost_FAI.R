source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# 14b_land_opportunity_cost_FAI.R  —  Thesis Project 4
# Replaces population density with provincial Fixed Asset Investment (FAI)
# as the land opportunity cost proxy.
#
# Data: 1951–2020_provincial_fixed_assets.xlsx
#   — Annual FAI per province, 100m RMB at current prices, 1951–2020
#   — Source: National Bureau of Statistics / compiled dataset
#
# Rationale: FAI per km² captures how intensively a province is investing in
# physical capital, which drives land demand and therefore opportunity cost.
# Unlike population density (static stock), FAI reflects the contemporaneous
# economic pressure on land at the time solar deployment decisions were made.
#
# Three metrics constructed:
#   fai_avg:        Mean FAI 2018–2020 (100m RMB/yr, current prices)
#   fai_per_km2:    fai_avg / province area (km²) — land-intensity proxy
#   fai_cagr_5yr:   CAGR 2015–2020 — captures recent investment acceleration
#
# Tests (parallel to 14_land_opportunity_cost.R):
#   A. Province-level FAI intensity vs mean deployment PVOUT
#   B. FAI intensity vs post-2020 siting rank change (Fig R8, R9)
#   C. OLS: does FAI explain province-level siting patterns?
#   D. Plant-level interaction: post2020 × log(FAI/km²)
#
# Requires: 00_setup.R + 11_siting_efficiency_regression.R sourced
#   (needs: china_provinces sf, prov_summary, deploy_en, pvout_national_median)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readxl)
  library(sf)
  library(ggplot2)
  library(ggrepel)
  library(fixest)
})

P4    <- FTS_ROOT
FAI_PATH <- file.path(FTS_ROOT, "data", "raw", "1951-2020_provincial_fixed_assets.xlsx")
ROB   <- file.path(P4, "robustness_outputs")
OUT   <- file.path(P4, "figures")
dir.create(ROB, showWarnings = FALSE)
dir.create(OUT, showWarnings = FALSE)

col_pre  <- "#2471A3"
col_post <- "#C0392B"

theme_rob <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        axis.title        = element_text(size = 10),
        plot.title        = element_text(size = 11, face = "bold"),
        plot.subtitle     = element_text(size = 9, colour = "grey40"),
        plot.caption      = element_text(size = 8, colour = "grey50"),
        legend.position   = "bottom")

# =============================================================================
# §1  LOAD AND CLEAN FAI DATA
# =============================================================================
cat("Loading FAI data...\n")

fai_raw <- read_xlsx(path.expand(FAI_PATH),
                     sheet = "current price",
                     col_names = TRUE)

# Drop footnote row (last row, all-NA on year columns)
fai_raw <- fai_raw[!is.na(fai_raw[[3]]), ]

# Drop national aggregate
fai_raw <- fai_raw[fai_raw[[1]] != "the whoe nation", ]

# Rename first two columns
names(fai_raw)[1:2] <- c("prov_fai", "prov_zh")

# Year columns are stored as numeric column names; convert
year_cols <- as.character(1951:2020)
names(fai_raw)[3:ncol(fai_raw)] <- year_cols

cat(sprintf("  Provinces in FAI data: %d\n", nrow(fai_raw)))
cat(sprintf("  Years: %s – %s\n", year_cols[1], year_cols[length(year_cols)]))

# Name harmonisation: FAI names → analysis prov_en names
fai_name_map <- c(
  "Inner Mongolia" = "Neimenggu",
  "Tibet"          = "Xizang"
  # All others match directly
)

fai_raw <- fai_raw %>%
  mutate(prov_en = dplyr::recode(prov_fai, !!!fai_name_map))

# Spot-check
cat("\nName mapping check:\n")
print(fai_raw %>% dplyr::select(prov_fai, prov_en) %>%
        dplyr::filter(prov_fai != prov_en))

# =============================================================================
# §2  COMPUTE FAI METRICS
# =============================================================================
cat("\nComputing FAI metrics...\n")

# Average FAI 2018–2020 (period bracketing the shuangtan transition)
fai_metrics <- fai_raw %>%
  mutate(
    fai_2018 = as.numeric(.data[["2018"]]),
    fai_2019 = as.numeric(.data[["2019"]]),
    fai_2020 = as.numeric(.data[["2020"]]),
    fai_2015 = as.numeric(.data[["2015"]]),
    fai_avg   = (fai_2018 + fai_2019 + fai_2020) / 3,  # 100m RMB/yr
    # 5-yr CAGR: captures how fast investment was accelerating pre-shuangtan
    fai_cagr_5yr = (fai_2020 / fai_2015)^(1/5) - 1
  ) %>%
  dplyr::select(prov_en, prov_fai, fai_avg, fai_cagr_5yr,
                fai_2015, fai_2018, fai_2019, fai_2020)

# =============================================================================
# §3  PROVINCE AREA FROM SHAPEFILE
# =============================================================================
cat("Computing province areas from shapefile...\n")

# Reproject to equal-area (CEA) for accurate area calculation
prov_cea <- st_transform(china_provinces, "EPSG:6933")

# Province name column (same logic as 11_siting_efficiency_regression.R)
prov_name_col <- "name"   # confirmed English name column from 11_*.R

prov_areas <- data.frame(
  prov_en      = st_drop_geometry(prov_cea)[[prov_name_col]],
  area_km2     = as.numeric(st_area(prov_cea)) / 1e6
) %>%
  group_by(prov_en) %>%
  summarise(area_km2 = sum(area_km2), .groups = "drop")

cat(sprintf("  Provinces with area: %d\n", nrow(prov_areas)))

# =============================================================================
# §4  MERGE AND BUILD PROVINCE-LEVEL ANALYSIS DATASET
# =============================================================================
prov_fai <- fai_metrics %>%
  left_join(prov_areas, by = "prov_en") %>%
  mutate(
    fai_per_km2     = fai_avg / area_km2,        # 100m RMB per km²
    log_fai_per_km2 = log(fai_per_km2),
    log_fai_avg     = log(fai_avg)
  ) %>%
  dplyr::filter(!is.na(area_km2), is.finite(fai_per_km2))

cat("\nFAI per km² summary (100m RMB/km²/yr):\n")
print(prov_fai %>%
        dplyr::select(prov_en, fai_avg, area_km2, fai_per_km2, fai_cagr_5yr) %>%
        mutate(across(where(is.numeric), ~round(., 3))) %>%
        arrange(desc(fai_per_km2)))

write.csv(prov_fai, file.path(ROB, "province_fai_metrics.csv"), row.names = FALSE)
cat("Saved province_fai_metrics.csv\n")

# =============================================================================
# §5  MERGE WITH PROVINCE SITING SUMMARY
# =============================================================================
if (!exists("prov_summary")) stop("Run 11_siting_efficiency_regression.R first.")

prov_joint <- prov_summary %>%
  left_join(prov_fai %>%
              dplyr::select(prov_en, fai_avg, fai_per_km2, log_fai_per_km2,
                            log_fai_avg, fai_cagr_5yr, area_km2),
            by = "prov_en") %>%
  dplyr::filter(!is.na(fai_per_km2)) %>%
  mutate(
    # Province-level pre-2020 mean deployment PVOUT (from prov_summary)
    # rank_pre / rank_post already in prov_summary from 11_*.R
    rank_decline = delta_rank < 0
  )

cat(sprintf("\nProvinces in joint analysis: %d\n", nrow(prov_joint)))

# =============================================================================
# §6  SECTION A — FAI INTENSITY VS MEAN DEPLOYMENT PVOUT
# =============================================================================
cat("\n===== A: FAI intensity vs deployment PVOUT =====\n")

cor_fai_rank_pre  <- cor(prov_joint$log_fai_per_km2, prov_joint$rank_pre,
                         use = "complete.obs")
cor_fai_rank_post <- cor(prov_joint$log_fai_per_km2, prov_joint$rank_post,
                         use = "complete.obs")

cat(sprintf("  Correlation: log(FAI/km²) vs pre-2020 mean rank:  r = %.3f\n",
            cor_fai_rank_pre))
cat(sprintf("  Correlation: log(FAI/km²) vs post-2020 mean rank: r = %.3f\n",
            cor_fai_rank_post))

# OLS: rank_pre ~ log(FAI/km²)
m_fai_rank <- lm(rank_pre ~ log_fai_per_km2, data = prov_joint,
                 weights = cap_total)
cat(sprintf("\n  OLS: rank_pre ~ log(FAI/km²)\n"))
cat(sprintf("    beta = %.3f  (SE %.3f, p = %.4f)  R² = %.3f\n",
            coef(m_fai_rank)[2],
            sqrt(diag(vcov(m_fai_rank)))[2],
            summary(m_fai_rank)$coefficients[2, 4],
            summary(m_fai_rank)$r.squared))

# =============================================================================
# §7  SECTION B — FAI INTENSITY VS POST-2020 RANK CHANGE (Figs R8, R9)
# =============================================================================
cat("\n===== B: FAI intensity vs post-2020 siting change =====\n")

cor_fai_delta <- cor(prov_joint$log_fai_per_km2, prov_joint$delta_rank,
                     use = "complete.obs")
cat(sprintf("  Correlation: log(FAI/km²) vs Δrank: r = %.3f\n", cor_fai_delta))

m_fai_delta1 <- lm(delta_rank ~ log_fai_per_km2,
                   data = prov_joint, weights = cap_total)
m_fai_delta2 <- lm(delta_rank ~ log_fai_per_km2 + rank_pre + log(cap_total),
                   data = prov_joint, weights = cap_total)
m_fai_cagr   <- lm(delta_rank ~ fai_cagr_5yr + rank_pre + log(cap_total),
                   data = prov_joint, weights = cap_total)

cat(sprintf("\n  OLS: delta_rank ~ log(FAI/km²) [bivariate]\n"))
cat(sprintf("    beta = %.3f  (SE %.3f, p = %.4f)  R² = %.3f\n",
            coef(m_fai_delta1)[2],
            sqrt(diag(vcov(m_fai_delta1)))[2],
            summary(m_fai_delta1)$coefficients[2, 4],
            summary(m_fai_delta1)$r.squared))

cat(sprintf("\n  OLS: delta_rank ~ log(FAI/km²) + rank_pre + log(cap) [controlled]\n"))
cat(sprintf("    beta = %.3f  (SE %.3f, p = %.4f)  R² = %.3f\n",
            coef(m_fai_delta2)[2],
            sqrt(diag(vcov(m_fai_delta2)))[2],
            summary(m_fai_delta2)$coefficients[2, 4],
            summary(m_fai_delta2)$r.squared))

cat(sprintf("\n  OLS: delta_rank ~ FAI CAGR 2015-2020 + rank_pre + log(cap)\n"))
cat(sprintf("    beta = %.3f  (SE %.3f, p = %.4f)  R² = %.3f\n",
            coef(m_fai_cagr)[2],
            sqrt(diag(vcov(m_fai_cagr)))[2],
            summary(m_fai_cagr)$coefficients[2, 4],
            summary(m_fai_cagr)$r.squared))

# --- Figure R8: FAI/km² vs mean deployment rank (pre and post) ---
cat("Building Figure R8...\n")

# Long format for pre/post
rank_long <- prov_joint %>%
  dplyr::select(prov_en, log_fai_per_km2, cap_total,
                rank_pre, rank_post, fai_cagr_5yr) %>%
  pivot_longer(cols = c(rank_pre, rank_post),
               names_to = "period", values_to = "mean_rank") %>%
  mutate(period = ifelse(period == "rank_pre",
                         "Pre-Q2 2020", "Post-Q2 2020"),
         period = factor(period, levels = c("Pre-Q2 2020", "Post-Q2 2020")))

fig_r8 <- ggplot(rank_long,
                 aes(x = log_fai_per_km2, y = mean_rank,
                     colour = period)) +
  geom_hline(yintercept = 50, linetype = "dotted",
             colour = "grey55", linewidth = 0.4) +
  geom_smooth(method = "lm", aes(weight = cap_total),
              se = TRUE, alpha = 0.15, linewidth = 0.9) +
  geom_point(aes(size = cap_total / 1e3), alpha = 0.8) +
  geom_text_repel(
    data = dplyr::filter(rank_long,
                         cap_total > 25 | mean_rank > 85 | mean_rank < 35),
    aes(label = prov_en),
    size = 2.4, max.overlaps = 15,
    box.padding = 0.3, segment.colour = "grey70",
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c("Pre-Q2 2020"  = col_pre,
               "Post-Q2 2020" = col_post),
    name = NULL
  ) +
  scale_size_continuous("Total cap (GW)", range = c(1.5, 7), guide = "none") +
  scale_x_continuous("log(Fixed Asset Investment per km²)  [100m RMB/km²/yr, 2018–2020 avg]") +
  scale_y_continuous(
    "Cap-weighted mean within-province PVOUT rank",
    breaks = seq(20, 100, 20)
  ) +
  labs(
    title    = "Figure R8: Provincial economic intensity vs siting quality",
    subtitle = paste0(
      "Each point = province × period. ",
      "If land costs explain siting quality, higher-FAI provinces should \n",
      "cluster lower on the y-axis. Dotted line = random-siting baseline (50th pct)."
    ),
    caption  = paste0(
      "FAI: fixed asset investment 2018–2020 average at current prices ",
      "(NBS; 100m RMB). Province area from GADM. ",
      "Rank: within-province PVOUT percentile (10km feasible-land grid)."
    )
  ) +
  theme_rob

ggsave(file.path(OUT, "fig_R8_fai_vs_rank.png"),
       fig_r8, width = 8, height = 6, dpi = 300)
cat("  Saved fig_R8_fai_vs_rank.png\n")

# --- Figure R9: FAI/km² vs delta_rank ---
cat("Building Figure R9...\n")

fig_r9 <- ggplot(prov_joint,
                 aes(x = log_fai_per_km2, y = delta_rank)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey60", linewidth = 0.4) +
  geom_smooth(method = "lm", aes(weight = cap_total),
              se = TRUE, colour = "#2C3E50",
              fill = "grey80", alpha = 0.3, linewidth = 0.9) +
  geom_point(aes(size = cap_total / 1e3,
                 colour = delta_rank < 0),
             alpha = 0.85) +
  geom_text_repel(
    data = dplyr::filter(prov_joint,
                         abs(delta_rank) > 8 |
                         log_fai_per_km2 > quantile(log_fai_per_km2, 0.85) |
                         cap_total > 30),
    aes(label = prov_en),
    size = 2.5, max.overlaps = 20,
    box.padding = 0.4,
    segment.colour = "grey60", segment.size = 0.3
  ) +
  # Quadrant annotations
  annotate("text",
           x = max(prov_joint$log_fai_per_km2, na.rm=TRUE) * 0.92,
           y = max(prov_joint$delta_rank, na.rm=TRUE) * 0.85,
           label = "High FAI\nBetter siting",
           hjust = 1, size = 2.5, colour = "grey55", fontface = "italic") +
  annotate("text",
           x = max(prov_joint$log_fai_per_km2, na.rm=TRUE) * 0.92,
           y = min(prov_joint$delta_rank, na.rm=TRUE) * 0.85,
           label = "High FAI\nWorse siting",
           hjust = 1, size = 2.5, colour = "#C0392B", fontface = "italic") +
  scale_colour_manual(
    values = c("FALSE" = col_pre, "TRUE" = col_post),
    labels = c("FALSE" = "Rank improved post-2020",
               "TRUE"  = "Rank declined post-2020"),
    name = NULL
  ) +
  scale_size_continuous("Total cap (GW)", range = c(2, 9)) +
  scale_x_continuous("log(Fixed Asset Investment per km²)  [100m RMB/km²/yr, 2018–2020 avg]") +
  scale_y_continuous(
    "Change in mean within-province rank (post minus pre-Q2 2020)\n← worse siting  |  better siting →"
  ) +
  labs(
    title    = "Figure R9: FAI intensity vs post-2020 siting rank change",
    subtitle = sprintf(
      "r = %.3f between log(FAI/km²) and Δrank (capacity-weighted OLS). %s",
      cor_fai_delta,
      if (abs(cor_fai_delta) < 0.2)
        "Near-zero correlation: land cost intensity does not predict post-2020 siting deterioration."
      else
        "See regression results for controlled estimates."
    ),
    caption  = paste0(
      "Province-level. FAI: NBS fixed asset investment 2018–2020 avg (100m RMB, current prices). ",
      sprintf("OLS (cap-weighted): β = %.3f (SE %.3f, p = %.3f).",
              coef(m_fai_delta1)[2],
              sqrt(diag(vcov(m_fai_delta1)))[2],
              summary(m_fai_delta1)$coefficients[2,4])
    )
  ) +
  theme_rob

ggsave(file.path(OUT, "fig_R9_fai_vs_rank_change.png"),
       fig_r9, width = 8.5, height = 6.5, dpi = 300)
cat("  Saved fig_R9_fai_vs_rank_change.png\n")


# =============================================================================
# §8  SECTION C — Plant-level interaction: post2020 × log(FAI/km²)
# =============================================================================
cat("\n===== C: Plant-level interaction =====\n")

deploy_fai <- deploy_en %>%
  dplyr::filter(is.finite(within_pct_rank), is.finite(cap_mw)) %>%
  left_join(prov_fai %>%
              dplyr::select(prov_en, log_fai_per_km2, fai_cagr_5yr),
            by = "prov_en") %>%
  dplyr::filter(is.finite(log_fai_per_km2)) %>%
  mutate(year_c = as.numeric(build_year) - 2020)

cat(sprintf("  Plants with FAI data: %d\n", nrow(deploy_fai)))

# Main spec (no interaction)
m_fai_main <- feols(
  within_pct_rank ~ post2020 + year_c | prov_en,
  data = deploy_fai, weights = ~cap_mw, vcov = ~prov_en
)

# Interaction: post2020 × log(FAI/km²)
# log_fai_per_km2 is province-level → absorbed by FE as main effect;
# the interaction post2020 × log_fai_per_km2 is identified because
# post2020 varies within province
m_fai_int <- feols(
  within_pct_rank ~ post2020 * log_fai_per_km2 + year_c | prov_en,
  data = deploy_fai, weights = ~cap_mw, vcov = ~prov_en
)

# Interaction with FAI CAGR
m_fai_cagr_int <- feols(
  within_pct_rank ~ post2020 * fai_cagr_5yr + year_c | prov_en,
  data = deploy_fai, weights = ~cap_mw, vcov = ~prov_en
)

cat("\nMain spec (reproduced):\n")
cat(sprintf("  beta post2020 = %.3f  (SE %.3f, p = %.4f)\n",
            coef(m_fai_main)["post2020TRUE"],
            se(m_fai_main)["post2020TRUE"],
            pvalue(m_fai_main)["post2020TRUE"]))

cat("\nInteraction: post2020 × log(FAI/km²):\n")
etable(m_fai_int, se.below = TRUE, coefstat = "se",
       digits = 3, fitstat = c("n", "r2"))

cat(sprintf("\n  beta post2020                       = %.3f (SE %.3f, p = %.4f)\n",
            coef(m_fai_int)["post2020TRUE"],
            se(m_fai_int)["post2020TRUE"],
            pvalue(m_fai_int)["post2020TRUE"]))
cat(sprintf("  beta post2020 × log(FAI/km²)        = %.3f (SE %.3f, p = %.4f)\n",
            coef(m_fai_int)["post2020TRUE:log_fai_per_km2"],
            se(m_fai_int)["post2020TRUE:log_fai_per_km2"],
            pvalue(m_fai_int)["post2020TRUE:log_fai_per_km2"]))

cat(sprintf("\n  beta post2020 × FAI CAGR 2015-2020  = %.3f (SE %.3f, p = %.4f)\n",
            coef(m_fai_cagr_int)["post2020TRUE:fai_cagr_5yr"],
            se(m_fai_cagr_int)["post2020TRUE:fai_cagr_5yr"],
            pvalue(m_fai_cagr_int)["post2020TRUE:fai_cagr_5yr"]))

# =============================================================================
# §9  SAVE FULL RESULTS TABLE
# =============================================================================
results_fai <- data.frame(
  Test = c(
    "A: log(FAI/km²) vs pre-2020 mean rank (Pearson r)",
    "A: log(FAI/km²) vs post-2020 mean rank (Pearson r)",
    "B: OLS delta_rank ~ log(FAI/km²) [bivariate] — beta",
    "B: OLS delta_rank ~ log(FAI/km²) [bivariate] — p",
    "B: OLS delta_rank ~ log(FAI/km²) [controlled] — beta",
    "B: OLS delta_rank ~ log(FAI/km²) [controlled] — p",
    "B: OLS delta_rank ~ FAI CAGR [controlled] — beta",
    "B: OLS delta_rank ~ FAI CAGR [controlled] — p",
    "C: post2020 × log(FAI/km²) interaction — beta",
    "C: post2020 × log(FAI/km²) interaction — p",
    "C: post2020 × FAI CAGR interaction — beta",
    "C: post2020 × FAI CAGR interaction — p"
  ),
  Value = round(c(
    cor_fai_rank_pre,
    cor_fai_rank_post,
    coef(m_fai_delta1)[2],
    summary(m_fai_delta1)$coefficients[2, 4],
    coef(m_fai_delta2)[2],
    summary(m_fai_delta2)$coefficients[2, 4],
    coef(m_fai_cagr)[2],
    summary(m_fai_cagr)$coefficients[2, 4],
    coef(m_fai_int)["post2020TRUE:log_fai_per_km2"],
    pvalue(m_fai_int)["post2020TRUE:log_fai_per_km2"],
    coef(m_fai_cagr_int)["post2020TRUE:fai_cagr_5yr"],
    pvalue(m_fai_cagr_int)["post2020TRUE:fai_cagr_5yr"]
  ), 4)
)
write.csv(results_fai, file.path(ROB, "fai_opportunity_cost_results.csv"),
          row.names = FALSE)

# =============================================================================
# §10  SUMMARY
# =============================================================================
cat("\n================================================================\n")
cat("SUMMARY: FAI-based land opportunity cost analysis\n")
cat("================================================================\n")
cat(sprintf("\nA. FAI intensity vs siting quality:\n"))
cat(sprintf("   log(FAI/km²) vs pre-2020 rank:  r = %.3f\n", cor_fai_rank_pre))
cat(sprintf("   log(FAI/km²) vs post-2020 rank: r = %.3f\n", cor_fai_rank_post))
cat(sprintf("\nB. FAI intensity vs post-2020 rank change:\n"))
cat(sprintf("   Bivariate:  beta = %.3f, p = %.4f, R² = %.3f\n",
            coef(m_fai_delta1)[2],
            summary(m_fai_delta1)$coefficients[2,4],
            summary(m_fai_delta1)$r.squared))
cat(sprintf("   Controlled: beta = %.3f, p = %.4f, R² = %.3f\n",
            coef(m_fai_delta2)[2],
            summary(m_fai_delta2)$coefficients[2,4],
            summary(m_fai_delta2)$r.squared))
cat(sprintf("\nC. Plant-level interaction:\n"))
cat(sprintf("   post2020 × log(FAI/km²): beta = %.3f, p = %.4f\n",
            coef(m_fai_int)["post2020TRUE:log_fai_per_km2"],
            pvalue(m_fai_int)["post2020TRUE:log_fai_per_km2"]))
cat(sprintf("   post2020 × FAI CAGR:     beta = %.3f, p = %.4f\n",
            coef(m_fai_cagr_int)["post2020TRUE:fai_cagr_5yr"],
            pvalue(m_fai_cagr_int)["post2020TRUE:fai_cagr_5yr"]))
int_sig <- pvalue(m_fai_int)["post2020TRUE:log_fai_per_km2"] < 0.10
cat(sprintf("\n   Interaction significant at 10%%? %s\n",
            ifelse(int_sig,
                   "YES — high-FAI provinces show differential post-2020 siting change",
                   "NO  — land opportunity costs (FAI) do not explain post-2020 siting shift")))
cat("================================================================\n")
cat("\nFigures: fig_R8_fai_vs_rank.png, fig_R9_fai_vs_rank_change.png\n")
cat("Tables:  province_fai_metrics.csv, fai_opportunity_cost_results.csv\n")
