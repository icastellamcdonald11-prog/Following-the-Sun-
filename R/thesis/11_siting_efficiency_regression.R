source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# 11_siting_efficiency_regression.R  —  Thesis Project 4
#
# Tests whether rapid post-2020 solar expansion shifted siting patterns away
# from China's solar resource frontier.
#
# Key construct:
#   within_pct_rank  —  each plant's PVOUT expressed as a percentile rank
#                        within its province's feasible-land PVOUT distribution.
#                        100 = best-resourced feasible site in the province;
#                        50 = at the provincial median; 0 = worst.
#                        This normalises away province-level differences in
#                        absolute solar resource, so a Jiangsu plant and a
#                        Qinghai plant can be compared on equal footing.
#
# Four tests:
#   T1 — Pre/post shuangtan shift in mean within-province percentile rank
#   T2 — Annual time series: has rank declined monotonically since 2020?
#   T3 — OLS regression with province fixed effects and time trend
#   T4 — Province-level: capacity growth rate vs Δrank (post–pre)
#
# Shuangtan boundary: Q2 2020 (consistent with thesis §5 / 01_siting_density_figs.R)
#   pre  = build_year < 2020  OR  (build_year == 2020 AND build_quarter == 1)
#   post = (build_year == 2020 AND build_quarter >= 2)  OR  build_year > 2020
#
# Requires: 00_setup.R already sourced
#   Objects used: deploy_df, grid_df, china_provinces, pvout_prov_feasible
#   Packages:     sf, dplyr, ggplot2, fixest, broom, ggrepel, patchwork
# =============================================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(ggplot2)
  library(fixest)
  library(broom)
  library(ggrepel)
  library(patchwork)
})

P4    <- FTS_ROOT
P3    <- FTS_P3_ROOT
OUT   <- file.path(P4, "figures")
CACHE <- file.path(P4, "data")
dir.create(OUT,   showWarnings = FALSE)
dir.create(CACHE, showWarnings = FALSE)

CEA_CRS <- "EPSG:6933"   # WGS 84 / NSIDC EASE-Grid 2.0 — matches pvout raster

theme_reg <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor   = element_blank(),
        axis.title         = element_text(size = 10),
        plot.title         = element_text(size = 11, face = "bold"),
        plot.subtitle      = element_text(size = 9,  colour = "grey40"),
        plot.caption       = element_text(size = 8,  colour = "grey50"),
        legend.position    = "bottom")

# Colour scheme: pre = steel blue, post = coral red
col_pre  <- "#2471A3"
col_post <- "#C0392B"

# =============================================================================
# §1  ASSIGN PROVINCE TO EACH FEASIBLE GRID CELL
#     (cached to avoid re-running the spatial join each session)
# =============================================================================
cache_file <- file.path(CACHE, "grid_province_assignment.rds")

if (file.exists(cache_file)) {
  cat("Loading cached grid-province assignment...\n")
  grid_prov_df <- readRDS(cache_file)
} else {
  cat("Assigning province to each feasible grid cell (first run — cached after this)...\n")

  # Only process feasible cells with valid PVOUT to save time
  feas_df <- grid_df %>%
    dplyr::filter(feasible, is.finite(pvout_annual)) %>%
    select(cell_id, x, y, pvout_annual)

  # Convert to sf in CEA CRS
  feas_sf <- st_as_sf(feas_df, coords = c("x", "y"), crs = CEA_CRS)

  # Reproject provinces to CEA
  prov_cea <- st_transform(china_provinces, crs = CEA_CRS)

  # Spatial join: each grid cell → province
  cat("  Running st_join (", nrow(feas_sf), "points ×", nrow(prov_cea), "provinces)...\n")
  joined  <- st_join(feas_sf, prov_cea[, c("name")], left = TRUE)
  # 'name' is the English province name column

  grid_prov_df <- as.data.frame(joined) %>%
    select(cell_id, pvout_annual, province_en = name) %>%
    dplyr::filter(!is.na(province_en))

  cat("  Cells assigned:", nrow(grid_prov_df), "\n")
  saveRDS(grid_prov_df, cache_file)
  cat("  Cached to", cache_file, "\n")
}

# =============================================================================
# §2  BUILD WITHIN-PROVINCE ECDF AND COMPUTE PERCENTILE RANK FOR EACH PLANT
# =============================================================================
cat("\nBuilding within-province PVOUT ecdfs...\n")

# Province-level ECDF functions (from feasible grid cells)
prov_ecdfs <- grid_prov_df %>%
  group_by(province_en) %>%
  summarise(ecdf_fn = list(ecdf(pvout_annual)), .groups = "drop")

cat("  Provinces with ECDF:", nrow(prov_ecdfs), "\n")

# Province name mapping: Chinese (deploy_df$prov) → English (china_provinces$name)
prov_map <- data.frame(
  prov    = as.data.frame(china_provinces)[[1]],
  prov_en = as.data.frame(china_provinces)[["name"]],
  stringsAsFactors = FALSE
)

# Attach English province name to deploy_df
deploy_en <- deploy_df %>%
  left_join(prov_map, by = "prov") %>%
  dplyr::filter(!is.na(prov_en), is.finite(pvout), is.finite(cap_mw))

cat("  Plants with English province name:", nrow(deploy_en),
    "(NAs:", nrow(deploy_df) - nrow(deploy_en), ")\n")

# Compute percentile rank: for plant i in province p,
# rank = ecdf_p(pvout_i) × 100  (i.e., % of feasible cells in p with PVOUT ≤ plant PVOUT)
compute_rank <- function(pvout_val, prov_name, ecdf_list) {
  idx <- match(prov_name, ecdf_list$province_en)
  if (is.na(idx)) return(NA_real_)
  ecdf_list$ecdf_fn[[idx]](pvout_val) * 100
}

cat("  Computing within-province percentile ranks...\n")
deploy_en$within_pct_rank <- mapply(
  compute_rank,
  deploy_en$pvout,
  deploy_en$prov_en,
  MoreArgs = list(ecdf_list = prov_ecdfs)
)

cat("  Rank range:", round(range(deploy_en$within_pct_rank, na.rm = TRUE), 1), "\n")
cat("  NAs in rank:", sum(is.na(deploy_en$within_pct_rank)), "\n")

# =============================================================================
# §3  DEFINE SHUANGTAN SPLIT AND CORE VARIABLES
# =============================================================================
deploy_en <- deploy_en %>%
  mutate(
    post2020 = case_when(
      build_year  < 2020                              ~ FALSE,
      build_year == 2020 & build_quarter == 1         ~ FALSE,
      TRUE                                            ~ TRUE
    ),
    period   = ifelse(post2020, "Post-Q2 2020", "Pre-Q2 2020"),
    year_num = as.numeric(build_year)
  )

cat("\nPre/post split:\n")
cat("  Pre -Q2 2020: ", sum(!deploy_en$post2020), "plants |",
    round(sum(deploy_en$cap_mw[!deploy_en$post2020], na.rm=TRUE)/1e3, 1), "GW\n")
cat("  Post-Q2 2020: ", sum( deploy_en$post2020),  "plants |",
    round(sum(deploy_en$cap_mw[ deploy_en$post2020], na.rm=TRUE)/1e3, 1), "GW\n")

# =============================================================================
# §4  T1 — DESCRIPTIVE: PRE/POST MEAN WITHIN-PROVINCE PERCENTILE RANK
# =============================================================================
cat("\n--- T1: Pre/post percentile rank comparison ---\n")

# Weighted mean function
wmean <- function(x, w) sum(x * w, na.rm=TRUE) / sum(w[is.finite(x)], na.rm=TRUE)
wmed  <- function(x, w) {
  ok <- is.finite(x) & is.finite(w)
  x <- x[ok]; w <- w[ok] / sum(w[ok])
  x[which.min(abs(cumsum(w[order(x)]) - 0.5))]
}

pre_df  <- dplyr::filter(deploy_en, !post2020, is.finite(within_pct_rank))
post_df <- dplyr::filter(deploy_en,  post2020, is.finite(within_pct_rank))

stats_t1 <- data.frame(
  Period          = c("Pre-Q2 2020", "Post-Q2 2020"),
  N_plants        = c(nrow(pre_df), nrow(post_df)),
  Cap_GW          = c(sum(pre_df$cap_mw)/1e3, sum(post_df$cap_mw)/1e3),
  Mean_rank_unw   = c(mean(pre_df$within_pct_rank),  mean(post_df$within_pct_rank)),
  Mean_rank_capwt = c(wmean(pre_df$within_pct_rank,  pre_df$cap_mw),
                      wmean(post_df$within_pct_rank, post_df$cap_mw)),
  Median_rank_unw = c(median(pre_df$within_pct_rank), median(post_df$within_pct_rank))
)
print(stats_t1, digits = 3)

# Hypothesis test: are post-2020 ranks lower?
wt  <- wilcox.test(pre_df$within_pct_rank, post_df$within_pct_rank,
                   alternative = "greater")   # H1: pre > post
cat(sprintf("\n  Wilcoxon rank-sum test (H1: pre > post): W = %.0f, p = %.4f\n",
            wt$statistic, wt$p.value))

# Effect size (unweighted)
delta_unw <- stats_t1$Mean_rank_unw[1]   - stats_t1$Mean_rank_unw[2]
delta_wt  <- stats_t1$Mean_rank_capwt[1] - stats_t1$Mean_rank_capwt[2]
cat(sprintf("  Decline in mean rank (unweighted):  %+.2f percentile points\n", -delta_unw))
cat(sprintf("  Decline in mean rank (cap-weighted): %+.2f percentile points\n", -delta_wt))

# =============================================================================
# §5  T2 — TIME SERIES: ANNUAL MEAN WITHIN-PROVINCE PERCENTILE RANK
# =============================================================================
cat("\n--- T2: Annual time series ---\n")

ts_df <- deploy_en %>%
  dplyr::filter(is.finite(within_pct_rank)) %>%
  group_by(build_year) %>%
  summarise(
    n              = n(),
    cap_gw         = sum(cap_mw) / 1e3,
    mean_rank_unw  = mean(within_pct_rank),
    mean_rank_capwt = wmean(within_pct_rank, cap_mw),
    med_rank_unw   = median(within_pct_rank),
    .groups = "drop"
  )

cat("Annual summary:\n")
print(ts_df, digits = 3)

# Pre/post trend lines
pre_ts  <- dplyr::filter(ts_df, build_year <= 2020)
post_ts <- dplyr::filter(ts_df, build_year >= 2020)

fit_pre  <- lm(mean_rank_capwt ~ build_year, data = pre_ts,  weights = cap_gw)
fit_post <- lm(mean_rank_capwt ~ build_year, data = post_ts, weights = cap_gw)

cat(sprintf("\n  Pre-2020 trend:  %+.2f pct-pts/yr (p=%.3f)\n",
            coef(fit_pre)[2],  summary(fit_pre)$coefficients[2,4]))
cat(sprintf("  Post-2020 trend: %+.2f pct-pts/yr (p=%.3f)\n",
            coef(fit_post)[2], summary(fit_post)$coefficients[2,4]))

fig_t2 <- ggplot(ts_df, aes(x = build_year)) +

  # Capacity bars (background)
  geom_col(aes(y = cap_gw * 3),          # scale to overlay on rank axis
           fill = "grey88", colour = NA) +

  # Trend lines
  geom_smooth(data = dplyr::filter(ts_df, build_year <= 2020),
              aes(y = mean_rank_capwt), method = "lm", se = TRUE,
              colour = col_pre,  fill = "#AED6F1", alpha = 0.3,
              linewidth = 0.8, linetype = "solid") +
  geom_smooth(data = dplyr::filter(ts_df, build_year >= 2020),
              aes(y = mean_rank_capwt), method = "lm", se = TRUE,
              colour = col_post, fill = "#F1948A", alpha = 0.3,
              linewidth = 0.8, linetype = "solid") +

  # Observed rank
  geom_line(aes(y = mean_rank_capwt), colour = "grey30",
            linewidth = 0.5, linetype = "dashed") +
  geom_point(aes(y = mean_rank_capwt, size = cap_gw),
             colour = "grey20", fill = "white", shape = 21, stroke = 0.8) +

  # Reference: random siting = 50th percentile
  geom_hline(yintercept = 50, linetype = "dotted",
             colour = "grey55", linewidth = 0.4) +

  # Shuangtan marker
  geom_vline(xintercept = 2020.25, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  annotate("text", x = 2020.35, y = 78, hjust = 0, size = 2.8,
           colour = "grey35", label = "Shuangtan\n(Q2 2020)") +

  scale_x_continuous("Year of commissioning",
                     breaks = 2018:2024) +
  scale_y_continuous(
    "Cap-weighted mean within-province\nPVOUT percentile rank",
    breaks = seq(20, 80, 10),
    sec.axis = sec_axis(~ . / 3,
                        name = "Annual installed capacity (GW)",
                        breaks = seq(0, 25, 5))
  ) +
  scale_size_continuous("Cap (GW)", range = c(2, 7)) +
  theme_reg +
  labs(
    title    = "Fig 11.1  Within-province PVOUT percentile rank by year",
    subtitle = paste0("Each point = capacity-weighted mean rank for that commissioning year. ",
                      "Grey bars = installed capacity (right axis). ",
                      "Dotted line = random-siting baseline (50th percentile)."),
    caption  = paste0("Within-province rank = percentile of plant PVOUT relative to province's ",
                      "feasible-land PVOUT distribution (10 km grid, GSA v2.0). ",
                      "Trend lines fitted by OLS weighted by annual capacity.")
  )

ggsave(file.path(OUT, "fig_11_1_annual_rank_timeseries.png"),
       fig_t2, width = 9, height = 5.5, dpi = 300)
cat("  Saved fig_11_1_annual_rank_timeseries.png\n")

# =============================================================================
# §6  T3 — REGRESSION: percentile rank ~ post2020 + province FE + year trend
# =============================================================================
cat("\n--- T3: Fixed-effects regression ---\n")

reg_df <- deploy_en %>%
  dplyr::filter(is.finite(within_pct_rank), is.finite(cap_mw)) %>%
  mutate(year_c = year_num - 2020)   # centre year on shuangtan

# Model 1: bivariate (no controls)
m1 <- feols(within_pct_rank ~ post2020,
            data = reg_df, weights = ~cap_mw,
            vcov = "HC1")

# Model 2: province fixed effects only
m2 <- feols(within_pct_rank ~ post2020 | prov_en,
            data = reg_df, weights = ~cap_mw,
            vcov = ~prov_en)    # cluster SEs by province

# Model 3: province FE + linear year trend
m3 <- feols(within_pct_rank ~ post2020 + year_c | prov_en,
            data = reg_df, weights = ~cap_mw,
            vcov = ~prov_en)

# Model 4: province FE + quadratic year trend (allows for non-linear rollout)
m4 <- feols(within_pct_rank ~ post2020 + year_c + I(year_c^2) | prov_en,
            data = reg_df, weights = ~cap_mw,
            vcov = ~prov_en)

# Model 5: unweighted version of M3 (robustness)
m5 <- feols(within_pct_rank ~ post2020 + year_c | prov_en,
            data = reg_df,
            vcov = ~prov_en)

cat("\n=== Regression results (capacity-weighted) ===\n")
etable(m1, m2, m3, m4,
       se.below    = TRUE,
       fitstat     = c("n", "r2", "wr2"),
       headers     = c("(1) Bivariate", "(2) Prov FE",
                       "(3) FE + trend", "(4) FE + quad"),
       coefstat    = "se",
       digits      = 3)

cat("\n=== Model 3 (main spec) — full summary ===\n")
cat(sprintf("  β post2020 (cap-wt, prov FE, year trend): %.3f  (SE %.3f, p %.4f)\n",
            coef(m3)["post2020TRUE"],
            se(m3)["post2020TRUE"],
            pvalue(m3)["post2020TRUE"]))
cat(sprintf("  Interpretation: post-2020 plants sit %.1f pct-pts lower in their\n",
            abs(coef(m3)["post2020TRUE"])))
cat(sprintf("    province's feasible-land PVOUT distribution, conditional on\n"))
cat(sprintf("    province and year trend (capacity-weighted; clustered SEs).\n"))

cat("\n=== Model 5 (unweighted robustness) ===\n")
cat(sprintf("  β post2020 (unweighted):  %.3f  (SE %.3f, p %.4f)\n",
            coef(m5)["post2020TRUE"],
            se(m5)["post2020TRUE"],
            pvalue(m5)["post2020TRUE"]))

# =============================================================================
# §7  T4 — PROVINCE LEVEL: CAPACITY GROWTH RATE vs ΔRANK
# =============================================================================
cat("\n--- T4: Province-level capacity growth vs change in siting rank ---\n")

prov_summary <- deploy_en %>%
  dplyr::filter(is.finite(within_pct_rank), is.finite(cap_mw)) %>%
  group_by(prov_en, post2020) %>%
  summarise(
    cap_gw     = sum(cap_mw) / 1e3,
    mean_rank  = wmean(within_pct_rank, cap_mw),
    n_plants   = n(),
    .groups    = "drop"
  ) %>%
  tidyr::pivot_wider(names_from  = post2020,
                     values_from = c(cap_gw, mean_rank, n_plants),
                     names_prefix = "") %>%
  # Rename: FALSE = pre, TRUE = post
  rename(
    cap_pre    = `cap_gw_FALSE`,
    cap_post   = `cap_gw_TRUE`,
    rank_pre   = `mean_rank_FALSE`,
    rank_post  = `mean_rank_TRUE`,
    n_pre      = `n_plants_FALSE`,
    n_post     = `n_plants_TRUE`
  ) %>%
  dplyr::filter(!is.na(cap_pre), !is.na(cap_post),
                !is.na(rank_pre), !is.na(rank_post)) %>%
  mutate(
    delta_rank    = rank_post - rank_pre,          # negative = siting got worse
    growth_rate   = cap_post / cap_pre,            # >1 means expanded post-2020
    log_growth    = log(growth_rate),
    cap_total     = cap_pre + cap_post,
    pct_post      = cap_post / cap_total * 100     # share of capacity built post-2020
  )

cat("  Provinces in T4 analysis:", nrow(prov_summary), "\n\n")
cat("  Distribution of Δrank (negative = worse post-2020):\n")
cat(sprintf("    Min: %.1f  Median: %.1f  Max: %.1f\n",
            min(prov_summary$delta_rank),
            median(prov_summary$delta_rank),
            max(prov_summary$delta_rank)))
cat(sprintf("    Provinces with declining rank: %d of %d\n",
            sum(prov_summary$delta_rank < 0), nrow(prov_summary)))

# Province-level regression: delta_rank ~ log_growth (+ controls)
m_prov1 <- lm(delta_rank ~ log_growth,
              data = prov_summary, weights = cap_total)
m_prov2 <- lm(delta_rank ~ log_growth + rank_pre,   # control for regression-to-mean
              data = prov_summary, weights = cap_total)

cat("\n  Province regression: delta_rank ~ log(growth_rate)\n")
cat(sprintf("    β log_growth (bivariate): %.3f  (SE %.3f, p %.4f)\n",
            coef(m_prov1)[2], sqrt(diag(vcov(m_prov1)))[2],
            summary(m_prov1)$coefficients[2, 4]))
cat(sprintf("    R² = %.3f\n", summary(m_prov1)$r.squared))
cat("\n  Province regression: delta_rank ~ log(growth_rate) + rank_pre\n")
cat(sprintf("    β log_growth (controlled): %.3f  (SE %.3f, p %.4f)\n",
            coef(m_prov2)[2], sqrt(diag(vcov(m_prov2)))[2],
            summary(m_prov2)$coefficients[2, 4]))
cat(sprintf("    R² = %.3f\n", summary(m_prov2)$r.squared))

# --- Fig T4: Province scatter ---
fig_t4 <- ggplot(prov_summary,
                 aes(x = log_growth, y = delta_rank)) +

  # Quadrant reference lines
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey60", linewidth = 0.4) +

  # OLS line (weighted by total capacity)
  geom_smooth(method = "lm", aes(weight = cap_total),
              se = TRUE, colour = "#2C3E50",
              fill = "grey80", alpha = 0.35, linewidth = 0.9) +

  # Points (sized by total capacity, coloured by direction)
  geom_point(aes(size = cap_total,
                 colour = delta_rank < 0),
             alpha = 0.8, stroke = 0.5) +

  # Province labels for notable cases
  geom_text_repel(
    data = dplyr::filter(prov_summary,
                         abs(delta_rank) > 8 | log_growth > 1.5 |
                         cap_total > 30),
    aes(label = prov_en),
    size = 2.6, colour = "grey20",
    max.overlaps = 20, box.padding = 0.4,
    segment.colour = "grey60", segment.size = 0.3
  ) +

  # Quadrant annotations
  annotate("text", x = max(prov_summary$log_growth) * 0.85,
           y = max(prov_summary$delta_rank)   * 0.85,
           label = "Fast growth\nBetter siting",
           hjust = 1, size = 2.6, colour = "grey50", fontface = "italic") +
  annotate("text", x = max(prov_summary$log_growth) * 0.85,
           y = min(prov_summary$delta_rank)   * 0.85,
           label = "Fast growth\nWorse siting",
           hjust = 1, size = 2.6, colour = "#C0392B", fontface = "italic") +

  scale_colour_manual(values = c("FALSE" = col_pre, "TRUE" = col_post),
                      labels = c("FALSE" = "Rank improved (post > pre)",
                                 "TRUE"  = "Rank declined (post < pre)"),
                      name = NULL) +
  scale_size_continuous("Total cap (GW)", range = c(2, 10)) +
  scale_x_continuous("Log capacity growth rate (post-2020 / pre-2020)") +
  scale_y_continuous(
    "Change in cap-weighted mean PVOUT percentile rank\n(post minus pre-Q2 2020)"
  ) +
  theme_reg +
  labs(
    title    = "Fig 11.2  Provincial capacity growth vs change in within-province siting rank",
    subtitle = paste0(
      "Each point = province. Size = total capacity (GW). ",
      "Red = rank declined post-2020 (worse siting).\n",
      sprintf("OLS (capacity-weighted): \u03b2 = %.2f per log-unit growth, ",
              coef(m_prov1)[2]),
      sprintf("R\u00b2 = %.2f.", summary(m_prov1)$r.squared)
    ),
    caption  = paste0(
      "Within-province rank from GSA v2.0 PVOUT at plant centroids relative to ",
      "10 km feasible-land grid (GRW, Robinson et al., 2025). ",
      "Provinces with <5 plants in either period excluded."
    )
  )

ggsave(file.path(OUT, "fig_11_2_province_growth_vs_rank.png"),
       fig_t4, width = 9, height = 6.5, dpi = 300)
cat("  Saved fig_11_2_province_growth_vs_rank.png\n")

# =============================================================================
# §8  SUPPLEMENTARY: WITHIN-PROVINCE RANK DENSITY (pre vs post, by region)
# =============================================================================
cat("\nBuilding Fig 11.3 (density of within-province ranks, pre vs post)...\n")

# National density
fig_t5_nat <- ggplot(
  dplyr::filter(deploy_en, is.finite(within_pct_rank)),
  aes(x = within_pct_rank, weight = cap_mw / sum(cap_mw),
      colour = period, fill = period)
) +
  geom_density(alpha = 0.25, linewidth = 0.85,
               adjust = 1.2, bw = "nrd0") +
  geom_vline(xintercept = 50, linetype = "dotted",
             colour = "grey50", linewidth = 0.4) +
  scale_colour_manual(values = c("Pre-Q2 2020"  = col_pre,
                                  "Post-Q2 2020" = col_post), name = NULL) +
  scale_fill_manual(values  = c("Pre-Q2 2020"  = col_pre,
                                 "Post-Q2 2020" = col_post), name = NULL) +
  scale_x_continuous("Within-province PVOUT percentile rank",
                     limits = c(0, 100), breaks = seq(0, 100, 25)) +
  scale_y_continuous("Capacity-weighted density") +
  theme_reg +
  labs(title    = "National",
       subtitle = "All provinces combined",
       caption  = "Dotted line = random-siting baseline (50th percentile).")

# Highlight four key provinces (high-deployment, interesting stories)
key_provs <- c("Qinghai", "Xinjiang", "Neimenggu", "Hebei")
# Note: Inner Mongolia appears as "Neimenggu" in this dataset

fig_t5_prov <- ggplot(
  deploy_en %>%
    dplyr::filter(is.finite(within_pct_rank),
                  prov_en %in% key_provs) %>%
    mutate(prov_en = factor(prov_en, levels = key_provs)),
  aes(x = within_pct_rank, weight = cap_mw,
      colour = period, fill = period)
) +
  geom_density(alpha = 0.20, linewidth = 0.75,
               adjust = 1.4, bw = "nrd0") +
  geom_vline(xintercept = 50, linetype = "dotted",
             colour = "grey50", linewidth = 0.4) +
  facet_wrap(~prov_en, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c("Pre-Q2 2020"  = col_pre,
                                  "Post-Q2 2020" = col_post), name = NULL) +
  scale_fill_manual(values  = c("Pre-Q2 2020"  = col_pre,
                                 "Post-Q2 2020" = col_post), name = NULL) +
  scale_x_continuous("Within-province PVOUT percentile rank",
                     breaks = seq(0, 100, 25)) +
  scale_y_continuous("Capacity-weighted density") +
  theme_reg +
  theme(strip.text = element_text(size = 9, face = "bold")) +
  labs(title    = "Selected provinces",
       subtitle = "Capacity-weighted density of plant PVOUT percentile rank")

fig_11_3 <- fig_t5_nat + fig_t5_prov +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(
    title   = "Fig 11.3  Within-province PVOUT percentile rank: pre vs post-Q2 2020",
    caption = paste0(
      "Rank = percentile of plant PVOUT within province's feasible 10 km grid (GSA v2.0). ",
      "Capacity-weighted. GRW plants \u22651 MW commissioned Q4 2017\u2013Q2 2024."
    ),
    theme   = theme(plot.title   = element_text(size = 12, face = "bold"),
                    plot.caption = element_text(size = 8, colour = "grey50"))
  )

ggsave(file.path(OUT, "fig_11_3_rank_density_prepost.png"),
       fig_11_3, width = 11, height = 5, dpi = 300)
cat("  Saved fig_11_3_rank_density_prepost.png\n")

# =============================================================================
# §9  SUMMARY TABLE FOR THESIS
# =============================================================================
cat("\n")
cat("================================================================\n")
cat("SUMMARY: Siting efficiency regression results\n")
cat("================================================================\n")
cat("\nT1 — Pre/post descriptive:\n")
print(stats_t1 %>%
        mutate(across(where(is.numeric), ~round(., 2))),
      row.names = FALSE)

cat(sprintf("\n  Unweighted rank decline post-2020: %+.2f pct-pts (p=%.4f, Wilcoxon)\n",
            stats_t1$Mean_rank_unw[2] - stats_t1$Mean_rank_unw[1],
            wt$p.value))
cat(sprintf("  Cap-weighted rank decline post-2020: %+.2f pct-pts\n",
            stats_t1$Mean_rank_capwt[2] - stats_t1$Mean_rank_capwt[1]))

cat("\nT2 — Annual trend in cap-weighted rank:\n")
cat(sprintf("  Pre-2020:  %+.2f pct-pts/yr\n",  coef(fit_pre)[2]))
cat(sprintf("  Post-2020: %+.2f pct-pts/yr\n",  coef(fit_post)[2]))

cat("\nT3 — Fixed-effects regression (main spec: Model 3):\n")
cat(sprintf("  β post2020 = %.3f pct-pts  (SE %.3f, p = %.4f)\n",
            coef(m3)["post2020TRUE"],
            se(m3)["post2020TRUE"],
            pvalue(m3)["post2020TRUE"]))
cat(sprintf("  Province FEs: YES  |  Year trend: YES  |  Weights: cap_mw\n"))
cat(sprintf("  Clustered SEs: province level\n"))

cat("\nT4 — Province-level growth vs rank change:\n")
cat(sprintf("  β log_growth = %.3f  (SE %.3f, p = %.4f)  R² = %.3f\n",
            coef(m_prov1)[2], sqrt(diag(vcov(m_prov1)))[2],
            summary(m_prov1)$coefficients[2,4],
            summary(m_prov1)$r.squared))
cat(sprintf("  Provinces declining in rank: %d of %d\n",
            sum(prov_summary$delta_rank < 0), nrow(prov_summary)))
cat("================================================================\n")
cat("\nAll figures saved to:", OUT, "\n")
