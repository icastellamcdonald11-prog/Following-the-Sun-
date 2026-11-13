source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# threshold_sensitivity.R: Supplementary Note 19 (choice of thresholds)
#
# Run from the project root:  Rscript paper/analysis/threshold_sensitivity.R
#
# A. Resource threshold: share of plants and of capacity below quantiles of
#    feasible-land PVOUT (median, top 25%, 20%, 10%); recomputes Table S8.
# B. Grid threshold: quadrant shares (fleet, Qinghai, Xinjiang) with the
#    penalty threshold at the 25th, 50th (main) and 75th percentiles across
#    plants, and with a distance rule of 50 km (the grid-accessible definition
#    used for land in Figs 1 and 4).
# C. Deterioration benchmark B in equation 12: feasible-land median (1,470),
#    capacity-weighted mean PVOUT of capacity installed before Q2 2020, and
#    4.22 kWh kWp-1 d-1 (1,540, thesis value). Main text: pre-2020 fleet mean (D9).
# Outputs: paper/analysis/threshold_sensitivity_{A,B,C}.csv
# =============================================================================

suppressPackageStartupMessages(library(dplyr))
p <- readRDS("fig_cache/plants_quad.rds")
g <- readRDS("fig_cache/grid_access.rds")
feas <- g$pvout_annual[g$feasible & is.finite(g$pvout_annual)]
m <- median(feas)
qn <- function(x) as.integer(substr(x, 1, 4)) * 4 + as.integer(substr(x, 6, 6))
p$q <- qn(p$year_quarter)
en <- function(x) dplyr::recode(x, Neimenggu = "Inner Mongolia", Xizang = "Tibet")

# ---- A. resource threshold ---------------------------------------------------
A <- do.call(rbind, lapply(c(0.50, 0.75, 0.80, 0.90), function(pr) {
  t <- unname(quantile(feas, pr))
  data.frame(threshold = sprintf("Q%d", round(100 * pr)), cutoff = round(t),
             plants_below = round(100 * mean(p$pvout < t), 1),
             capacity_below = round(100 * sum(p$cap_mw[p$pvout < t]) / sum(p$cap_mw), 1))
}))
write.csv(A, "outputs/analysis/threshold_sensitivity_A.csv", row.names = FALSE)
print(A, row.names = FALSE)

# ---- B. grid threshold -------------------------------------------------------
quad <- function(hi_res, near) ifelse(hi_res & near, "Premium", ifelse(hi_res, "Stranded",
                                ifelse(near, "Grid-convenient", "Dominated")))
shares <- function(qd, sel) {
  s <- tapply(p$cap_mw[sel], factor(qd[sel], c("Premium", "Stranded", "Grid-convenient", "Dominated")), sum)
  s[is.na(s)] <- 0; round(100 * s / sum(s), 1)
}
rules <- list(
  "Penalty, 25th percentile" = p$tx_penalty < quantile(p$tx_penalty, 0.25),
  "Penalty, median (main)"   = p$tx_penalty < median(p$tx_penalty),
  "Penalty, 75th percentile" = p$tx_penalty < quantile(p$tx_penalty, 0.75),
  "Distance <= 50 km"        = p$dist_km <= 50)
thr <- c(quantile(p$tx_penalty, c(0.25, 0.5, 0.75)), NA)
B <- do.call(rbind, Map(function(nm, near, th) {
  qd <- quad(p$pvout >= m, near)
  sets <- list(National = rep(TRUE, nrow(p)), Qinghai = p$prov == "Qinghai",
               Xinjiang = p$prov == "Xinjiang",
               `Xinjiang post-2020` = p$prov == "Xinjiang" & p$q > qn("2020Q2"))
  do.call(rbind, lapply(names(sets), function(k)
    data.frame(rule = nm, threshold = round(th, 1), group = k, t(shares(qd, sets[[k]])),
               check.names = FALSE)))
}, names(rules), rules, thr))
write.csv(B, "outputs/analysis/threshold_sensitivity_B.csv", row.names = FALSE)
print(B, row.names = FALSE)

# ---- C. deterioration benchmark ----------------------------------------------
post <- p[p$q >= qn("2020Q2"), ]
pre_mean <- weighted.mean(p$pvout[p$q < qn("2020Q2")], p$cap_mw[p$q < qn("2020Q2")])
bench <- c("Feasible-land median" = m, "Pre-2020 fleet mean (main)" = pre_mean, "4.22 kWh/kWp/d (thesis)" = 4.22 * 365)
C <- do.call(rbind, lapply(names(bench), function(nm) {
  B0 <- bench[[nm]]
  s <- post |> group_by(prov) |>
    summarise(k = sum(cap_mw) / 1e3, pv = weighted.mean(pvout, cap_mw), .groups = "drop") |>
    mutate(delta = k * (pv - B0))
  neg <- sum(s$delta[s$delta < 0])
  s <- s |> mutate(share = ifelse(delta < 0, 100 * delta / neg, NA)) |> arrange(desc(share))
  top <- head(s[!is.na(s$share), ], 4)
  data.frame(benchmark = nm, B = round(B0, 1),
             rank1 = sprintf("%s %.1f", en(top$prov[1]), top$share[1]),
             rank2 = sprintf("%s %.1f", en(top$prov[2]), top$share[2]),
             rank3 = sprintf("%s %.1f", en(top$prov[3]), top$share[3]),
             rank4 = sprintf("%s %.1f", en(top$prov[4]), top$share[4]),
             top4 = round(sum(top$share), 1),
             xinjiang = round(s$share[s$prov == "Xinjiang"], 1),
             n_below = sum(!is.na(s$share)))
}))
write.csv(C, "outputs/analysis/threshold_sensitivity_C.csv", row.names = FALSE)
print(C, row.names = FALSE)

# exact decomposition check: with B = pre-2020 fleet mean, sum_p Delta_p / K_total
# equals the change in the fleet's capacity-weighted PVOUT
d_fleet <- weighted.mean(p$pvout, p$cap_mw) - pre_mean
s <- post |> group_by(prov) |> summarise(k = sum(cap_mw), pv = weighted.mean(pvout, cap_mw))
cat(sprintf("Change in fleet mean %.3f; sum of contributions / total capacity %.3f\n",
            d_fleet, sum(s$k * (s$pv - pre_mean)) / sum(p$cap_mw)))
