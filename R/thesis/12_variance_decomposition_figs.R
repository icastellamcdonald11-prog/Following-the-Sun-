source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# 12_variance_decomposition_figs.R  —  Thesis Project 4
# Robustness Tests 1–3: Variance decomposition, quartile concentration,
# and stochastic dominance.
#
# Produces:
#   Figure R1  — Quarterly CV of new-build PVOUT (Q1 2018 – Q2 2024)
#   Figure R2  — Quarterly share of new capacity by PVOUT quartile
#   Table R1   — KS test (pre vs. post Q2 2020 PVOUT distributions)
#
# Requires: source("00_setup.R")  (deploy_df, pvout_national_median,
#           pvout_q75_national, grid_df)
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

P4  <- FTS_ROOT
OUT <- file.path(P4, "figures")
dir.create(OUT, showWarnings = FALSE)

ROB <- file.path(P4, "robustness_outputs")
dir.create(ROB, showWarnings = FALSE)

# =============================================================================
# Shared theme (matches 01_siting_density_figs.R)
# =============================================================================
theme_rob <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.line.x  = element_line(colour = "grey30", linewidth = 0.3),
    axis.line.y  = element_line(colour = "grey30", linewidth = 0.3),
    axis.ticks   = element_line(colour = "grey30", linewidth = 0.3),
    axis.title   = element_text(size = 10),
    legend.position   = "bottom",
    legend.key.width  = unit(0.8, "cm"),
    legend.key.height = unit(0.35, "cm"),
    legend.text       = element_text(size = 8.5),
    plot.title   = element_text(size = 11, face = "bold"),
    plot.caption = element_text(size = 8, colour = "grey50")
  )

# Quarter ordering: sorted lexicographically (format is "YYYYQn", e.g. "2020Q2")
all_quarters <- sort(unique(deploy_df$year_quarter))

# Pre / post shuangtan boundary: post = Q2 2020 inclusive
# Use build_year / build_quarter (integers) to avoid string-format ambiguity.
# Called on deploy_df rows; returns logical vector.
is_post_row <- function(yr, qn) yr > 2020 | (yr == 2020 & qn >= 2)

# =============================================================================
# SECTION 1 — Capacity-weighted CV per quarter
# =============================================================================
cat("\n===== R1: Quarterly CV (capacity-weighted) =====\n")

# Capacity-weighted mean and variance
wt_stats <- function(x, w) {
  w  <- w / sum(w)
  mu <- sum(w * x)
  v  <- sum(w * (x - mu)^2)
  list(mean = mu, sd = sqrt(v))
}

quarterly_cv <- deploy_df %>%
  filter(is.finite(pvout), is.finite(cap_mw), cap_mw > 0) %>%
  group_by(year_quarter, build_year, build_quarter) %>%
  summarise(
    n        = n(),
    cap_total = sum(cap_mw),
    cw_mean  = {s <- wt_stats(pvout, cap_mw); s$mean},
    cw_sd    = {s <- wt_stats(pvout, cap_mw); s$sd},
    .groups  = "drop"
  ) %>%
  mutate(
    cv     = cw_sd / cw_mean,
    period = ifelse(is_post_row(build_year, build_quarter), "Post Q2 2020", "Pre Q2 2020"),
    qfac   = factor(year_quarter, levels = all_quarters)
  ) %>%
  # Drop quarters with fewer than 5 plants (too few to compute meaningful CV)
  filter(n >= 5) %>%
  arrange(qfac)

cat("Quarterly CV summary:\n")
print(quarterly_cv %>% select(year_quarter, n, cw_mean, cw_sd, cv, period))

# Pre vs. post mean CV
cv_period_means <- quarterly_cv %>%
  group_by(period) %>%
  summarise(mean_cv = round(mean(cv), 4), .groups = "drop")
cat("\nMean CV by period:\n")
print(cv_period_means)

# --- Figure R1 ---
cat("Building Figure R1...\n")

# x-axis: numeric index for spacing; labels placed at every 4th quarter
n_q    <- nrow(quarterly_cv)
x_idx  <- seq_len(n_q)
x_brk  <- x_idx[seq(1, n_q, by = 4)]
x_lbl  <- quarterly_cv$year_quarter[x_brk]

# Vertical line position: first post-shuangtan quarter (Q2 2020)
vline_idx <- which(quarterly_cv$period == "Post Q2 2020")[1]

quarterly_cv$x <- x_idx

# Smoothed trend (LOESS) per period for visual clarity
pre_df  <- quarterly_cv %>% filter(period == "Pre Q2 2020")
post_df <- quarterly_cv %>% filter(period == "Post Q2 2020")

fig_r1 <- ggplot(quarterly_cv, aes(x = x, y = cv)) +
  # Background shading for post period
  annotate("rect",
           xmin = vline_idx - 0.5, xmax = max(x_idx) + 0.5,
           ymin = -Inf, ymax = Inf,
           fill = "#F5E6E6", alpha = 0.5) +
  # Raw points
  geom_point(aes(colour = period), size = 1.6, alpha = 0.8) +
  # LOESS trend lines per period
  geom_smooth(data = pre_df,  aes(colour = "Pre Q2 2020"),
              method = "loess", span = 0.6, se = FALSE,
              linewidth = 0.9, linetype = "solid") +
  geom_smooth(data = post_df, aes(colour = "Post Q2 2020"),
              method = "loess", span = 0.6, se = FALSE,
              linewidth = 0.9, linetype = "solid") +
  # Shuangtan boundary
  geom_vline(xintercept = vline_idx - 0.5,
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  annotate("text", x = vline_idx - 0.3, y = Inf,
           label = "Post-shuangtan", hjust = 0, vjust = 1.5,
           size = 2.8, colour = "grey35") +
  scale_colour_manual(
    name   = NULL,
    values = c("Pre Q2 2020" = "#2C6FAC", "Post Q2 2020" = "#C0392B")
  ) +
  scale_x_continuous(
    breaks = x_brk, labels = x_lbl,
    expand = c(0.02, 0)
  ) +
  scale_y_continuous(
    name   = "Coefficient of variation (CV = σ/μ)",
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  labs(
    x       = NULL,
    title   = "Figure R1: Dispersion of new-build PVOUT over time",
    caption = "Capacity-weighted CV per quarter. LOESS trend lines shown separately for pre- and post-shuangtan periods."
  ) +
  theme_rob +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5))

ggsave(file.path(OUT, "fig_R1_cv_quarterly.png"),
       fig_r1, width = 7, height = 4, dpi = 300)
cat("  Saved fig_R1_cv_quarterly.png\n")


# =============================================================================
# SECTION 2 — Quartile concentration analysis
# =============================================================================
cat("\n===== R2: Quartile concentration =====\n")

# National feasible-land PVOUT quartile thresholds
feasible_pvout_vec <- grid_df$pvout_annual[
  grid_df$feasible & is.finite(grid_df$pvout_annual) & grid_df$pvout_annual > 0
]
pvout_q25_national <- quantile(feasible_pvout_vec, 0.25)

cat(sprintf("National feasible PVOUT quartile thresholds:\n"))
cat(sprintf("  Q25: %.2f  |  Q50 (median): %.2f  |  Q75: %.2f\n",
            pvout_q25_national, pvout_national_median, pvout_q75_national))

# Assign each plant to a national PVOUT quartile
deploy_df_q <- deploy_df %>%
  filter(is.finite(pvout), is.finite(cap_mw), cap_mw > 0) %>%
  mutate(
    pvout_quartile = case_when(
      pvout <  pvout_q25_national    ~ "Q1 (lowest)",
      pvout <  pvout_national_median ~ "Q2",
      pvout <  pvout_q75_national    ~ "Q3",
      TRUE                            ~ "Q4 (highest)"
    ),
    pvout_quartile = factor(pvout_quartile,
                            levels = c("Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"))
  )

# Quarterly share of capacity by quartile
quartile_shares <- deploy_df_q %>%
  group_by(year_quarter, build_year, build_quarter, pvout_quartile) %>%
  summarise(cap_q = sum(cap_mw), .groups = "drop") %>%
  group_by(year_quarter) %>%
  mutate(
    cap_total = sum(cap_q),
    share     = cap_q / cap_total * 100
  ) %>%
  ungroup() %>%
  mutate(
    qfac   = factor(year_quarter, levels = all_quarters),
    period = ifelse(is_post_row(build_year, build_quarter), "Post Q2 2020", "Pre Q2 2020")
  )

# Period-mean shares for annotation (pre vs. post)
period_shares <- quartile_shares %>%
  group_by(period, pvout_quartile) %>%
  summarise(mean_share = round(mean(share), 1), .groups = "drop")

cat("\nMean capacity share by PVOUT quartile and period:\n")
print(period_shares %>% pivot_wider(names_from = period, values_from = mean_share))

# --- Figure R2 ---
cat("Building Figure R2...\n")

# Colour palette for quartiles: Q1 (light/warm) → Q4 (deep blue)
q_cols <- c(
  "Q1 (lowest)"  = "#D9534F",   # warm red
  "Q2"           = "#F0AD4E",   # amber
  "Q3"           = "#5BC0DE",   # light blue
  "Q4 (highest)" = "#1A3A5C"    # deep navy
)

# x-axis: use numeric index again to ensure even spacing
all_q_levels <- levels(factor(quartile_shares$qfac))
q_xmap <- data.frame(
  year_quarter = all_quarters,
  x            = seq_along(all_quarters)
)
quartile_shares <- quartile_shares %>%
  left_join(q_xmap, by = "year_quarter")

x_brk2 <- q_xmap$x[seq(1, nrow(q_xmap), by = 4)]
x_lbl2 <- q_xmap$year_quarter[seq(1, nrow(q_xmap), by = 4)]
# First post-shuangtan quarter in the x-map
post_quarters2 <- quartile_shares$year_quarter[quartile_shares$period == "Post Q2 2020"]
first_post_q2  <- sort(unique(post_quarters2))[1]
vline_idx2     <- q_xmap$x[q_xmap$year_quarter == first_post_q2]

fig_r2 <- ggplot(quartile_shares,
                 aes(x = x, y = share, fill = pvout_quartile)) +
  # Background shading post period
  annotate("rect",
           xmin = vline_idx2 - 0.5, xmax = max(quartile_shares$x) + 0.5,
           ymin = -Inf, ymax = Inf,
           fill = "#F5F5F5", alpha = 0.6) +
  geom_bar(stat = "identity", width = 0.85, colour = NA) +
  geom_vline(xintercept = vline_idx2 - 0.5,
             linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  annotate("text", x = vline_idx2 - 0.3, y = 101,
           label = "Post-shuangtan", hjust = 0, vjust = 1,
           size = 2.8, colour = "grey35") +
  scale_fill_manual(
    name   = "PVOUT quartile\n(national feasible land)",
    values = q_cols
  ) +
  scale_x_continuous(
    breaks = x_brk2, labels = x_lbl2,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    name   = "Share of new quarterly capacity (%)",
    limits = c(0, 105),
    breaks = seq(0, 100, 25),
    expand = c(0, 0)
  ) +
  labs(
    x       = NULL,
    title   = "Figure R2: PVOUT quartile concentration of new-build capacity",
    caption = paste0(
      "Capacity-weighted quartile shares per quarter. Quartile thresholds based on feasible national land distribution\n",
      sprintf("(Q25=%.0f, Median=%.0f, Q75=%.0f kWh/kWp/yr).",
              pvout_q25_national, pvout_national_median, pvout_q75_national)
    )
  ) +
  theme_rob +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5),
    legend.position = "right"
  )

ggsave(file.path(OUT, "fig_R2_quartile_shares.png"),
       fig_r2, width = 8, height = 4.5, dpi = 300)
cat("  Saved fig_R2_quartile_shares.png\n")


# =============================================================================
# SECTION 3 — KS test for first-order stochastic dominance
# =============================================================================
cat("\n===== R3: KS test (FOSD) =====\n")

# Split capacity-weighted distributions
# For the KS test we replicate each plant observation by (cap_mw / min_cap) to
# produce an integer-weighted empirical distribution.
# More directly: use ks.test() on the plant-level PVOUT vectors, then note that
# the two-sample KS test tests whether one ECDF is ≥ the other everywhere.
# FOSD: pre dominates post if F_pre(x) ≤ F_post(x) for all x (fewer mass below).
# One-sided KS: alternative = "less" tests H₀: F_pre ≥ F_post vs Ha: F_pre < F_post
# i.e., post puts more mass below every x → post is stochastically dominated by pre.

ok      <- is.finite(deploy_df$pvout) & is.finite(deploy_df$cap_mw) & deploy_df$cap_mw > 0
post_lgl <- is_post_row(deploy_df$build_year, deploy_df$build_quarter)

pre_pvout  <- deploy_df$pvout[ok & !post_lgl]
post_pvout <- deploy_df$pvout[ok &  post_lgl]
pre_cap    <- deploy_df$cap_mw[ok & !post_lgl]
post_cap   <- deploy_df$cap_mw[ok &  post_lgl]

cat(sprintf("Pre-Q2 2020 plants:  %d (%.1f GW)\n",
            length(pre_pvout), sum(pre_cap) / 1e3))
cat(sprintf("Post-Q2 2020 plants: %d (%.1f GW)\n",
            length(post_pvout), sum(post_cap) / 1e3))

# --- Unweighted KS test (plant count) ---
ks_uw <- ks.test(post_pvout, pre_pvout, alternative = "less")
cat(sprintf("\nUnweighted KS test (post vs. pre, alternative = 'less'):\n"))
cat(sprintf("  D = %.4f,  p = %.2e\n", ks_uw$statistic, ks_uw$p.value))

# --- Capacity-weighted KS test ---
# Expand to pseudo-sample by rounding cap to nearest 0.1 MW and replicating
# (avoids distortion from very large plants while retaining weighting)
expand_weighted <- function(pvout_vec, cap_vec, unit = 0.1) {
  n_rep <- pmax(round(cap_vec / unit), 1L)
  rep(pvout_vec, times = n_rep)
}
pre_exp  <- expand_weighted(pre_pvout,  pre_cap)
post_exp <- expand_weighted(post_pvout, post_cap)

cat(sprintf("\nCapacity-weighted expanded samples: pre=%d, post=%d\n",
            length(pre_exp), length(post_exp)))

ks_cw <- ks.test(post_exp, pre_exp, alternative = "less")
cat(sprintf("\nCapacity-weighted KS test (post vs. pre, alternative = 'less'):\n"))
cat(sprintf("  D = %.4f,  p = %.2e\n", ks_cw$statistic, ks_cw$p.value))

# Two-sided KS (for completeness)
ks_cw_2s <- ks.test(post_exp, pre_exp)
cat(sprintf("\nCapacity-weighted KS test (two-sided):\n"))
cat(sprintf("  D = %.4f,  p = %.2e\n", ks_cw_2s$statistic, ks_cw_2s$p.value))

# --- Table R1 ---
table_r1 <- data.frame(
  Test          = c("Unweighted (plant count)",
                    "Capacity-weighted",
                    "Capacity-weighted (two-sided)"),
  Pre_N         = c(length(pre_pvout),  length(pre_exp),  length(pre_exp)),
  Post_N        = c(length(post_pvout), length(post_exp), length(post_exp)),
  Alternative   = c("Post < Pre (FOSD)", "Post < Pre (FOSD)", "Two-sided"),
  D_statistic   = round(c(ks_uw$statistic, ks_cw$statistic, ks_cw_2s$statistic), 4),
  p_value       = formatC(c(ks_uw$p.value, ks_cw$p.value, ks_cw_2s$p.value),
                          format = "e", digits = 2),
  Interpretation = c(
    ifelse(ks_uw$p.value   < 0.01, "Reject H0 at 1%", "Fail to reject H0"),
    ifelse(ks_cw$p.value   < 0.01, "Reject H0 at 1%", "Fail to reject H0"),
    ifelse(ks_cw_2s$p.value < 0.01, "Reject H0 at 1%", "Fail to reject H0")
  )
)

cat("\nTable R1:\n")
print(table_r1)
write.csv(table_r1, file.path(ROB, "table_R1_ks_test.csv"), row.names = FALSE)
cat("  Saved table_R1_ks_test.csv\n")

# Also save quarterly CV and quartile share tables
write.csv(quarterly_cv %>% select(-x),
          file.path(ROB, "fig_R1_quarterly_cv.csv"), row.names = FALSE)
write.csv(quartile_shares %>% select(-x),
          file.path(ROB, "fig_R2_quartile_shares.csv"), row.names = FALSE)

cat("\nAll done. Figures in figures/, tables in robustness_outputs/\n")
