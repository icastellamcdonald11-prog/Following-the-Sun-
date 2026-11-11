source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# deterioration_shares.R: provincial contributions to the post-2020 change
# (Methods eq. 12) at benchmark B = 4.22 kWh kWp-1 d-1 = 1,540.3 kWh kWp-1 yr-1
# (DECISIONS D6, 11 Nov 2026), with the bootstrap CI on the same benchmark.
#
# Run from the project root:  Rscript paper/analysis/deterioration_shares.R
#
#   Delta_p = K_p_post * (PVOUT_p_post_capweighted - B)
#   sigma_p = Delta_p / sum_{q: Delta_q < 0} Delta_q * 100      (Delta_p < 0)
# Post-2020 = commissioned in Q2 2020 or later (as in Methods).
# Bootstrap: resample plants with replacement within each province, 1,000
# times, recompute sigma_p; 2.5th and 97.5th percentiles. Seed 42.
# Outputs: paper/analysis/deterioration_shares.csv (all provinces),
#          paper/analysis/deterioration_bootstrap.csv (CIs, deteriorating provinces)
# =============================================================================

suppressPackageStartupMessages(library(dplyr))
B <- 4.22 * 365
p <- readRDS("fig_cache/plants_quad.rds")
qn <- function(x) as.integer(substr(x, 1, 4)) * 4 + as.integer(substr(x, 6, 6))
post <- p[qn(p$year_quarter) >= qn("2020Q2"), ]
en <- function(x) dplyr::recode(x, Neimenggu = "Inner Mongolia", Xizang = "Tibet")

shares <- function(d) {
  s <- d |> group_by(prov) |>
    summarise(cap_gw = sum(cap_mw) / 1e3, pv = weighted.mean(pvout, cap_mw), .groups = "drop") |>
    mutate(gap = pv - B, delta = cap_gw * gap)
  neg <- sum(s$delta[s$delta < 0])
  s |> mutate(share = ifelse(delta < 0, 100 * delta / neg, NA_real_))
}
s <- shares(post) |> arrange(desc(share)) |> mutate(province = en(prov))
write.csv(s |> select(province, cap_gw, pv, gap, delta, share) |>
            mutate(across(c(cap_gw, pv, gap, delta, share), \(x) round(x, 2))),
          "outputs/analysis/deterioration_shares.csv", row.names = FALSE)

set.seed(42)
boot <- replicate(1000, simplify = FALSE, {
  b <- post |> group_by(prov) |> slice_sample(prop = 1, replace = TRUE) |> ungroup()
  r <- shares(b); setNames(r$share, r$prov)
})
det <- s$prov[!is.na(s$share)]
ci <- t(sapply(det, function(pr) {
  v <- sapply(boot, function(x) if (pr %in% names(x)) x[[pr]] else NA_real_)
  v[is.na(v)] <- 0
  quantile(v, c(0.025, 0.975))
}))
out <- data.frame(province = en(det), share = round(s$share[match(det, s$prov)], 1),
                  ci_low = round(ci[, 1], 1), ci_high = round(ci[, 2], 1))
write.csv(out, "outputs/analysis/deterioration_bootstrap.csv", row.names = FALSE)

cat(sprintf("B = %.1f; post-2020 plants %d (%.1f GW)\n", B, nrow(post), sum(post$cap_mw) / 1e3))
print(head(out, 8), row.names = FALSE)
cat(sprintf("Top four share: %.1f%%; all others: %.1f%%\n",
            sum(head(out$share, 4)), 100 - sum(head(out$share, 4))))
cat("Net positive contributors:\n")
print(s |> filter(delta > 0) |> arrange(desc(delta)) |> select(province, cap_gw, pv, delta) |>
        mutate(across(where(is.numeric), \(x) round(x, 1))), row.names = FALSE)
