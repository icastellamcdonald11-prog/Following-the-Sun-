source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# si_distribution_tests.R: Supplementary Note 17 (distributional shift)
#
# Run from the project root:  Rscript paper/analysis/si_distribution_tests.R
#
# Compares PVOUT of capacity commissioned before Q2 2020 with capacity
# commissioned from Q2 2020 (the split used in 12_variance_decomposition_figs.R),
# unweighted and capacity-weighted, and reports the Kolmogorov-Smirnov
# statistic in both directions:
#   D_below = max_x [F_post(x) - F_pre(x)]  (post has more mass below x)
#   D_above = max_x [F_pre(x) - F_post(x)]  (post has more mass above x)
# Unweighted p-values from ks.test(); for the capacity-weighted case the
# statistic is computed exactly on weighted ECDFs. (12_variance_decomposition_
# figs.R obtained p-values by replicating plants in 0.1 MW units, which treats
# capacity as independent observations, so those p-values are not reported.)
# Output: paper/analysis/si_distribution_tests.csv
# =============================================================================

p <- readRDS("fig_cache/plants_quad.rds")
qn <- function(x) as.integer(substr(x, 1, 4)) * 4 + as.integer(substr(x, 6, 6))
post <- qn(p$year_quarter) >= qn("2020Q2"); pre <- !post

wecdf <- function(x, w) { o <- order(x); list(x = x[o], F = cumsum(w[o]) / sum(w)) }
Fat <- function(e, t) { i <- findInterval(t, e$x); ifelse(i == 0, 0, e$F[pmax(i, 1)]) }
grid_t <- sort(unique(p$pvout))

dstats <- function(w) {
  e_pre <- wecdf(p$pvout[pre], w[pre]); e_post <- wecdf(p$pvout[post], w[post])
  d <- Fat(e_post, grid_t) - Fat(e_pre, grid_t)
  c(D_below = max(d), x_below = grid_t[which.max(d)],
    D_above = max(-d), x_above = grid_t[which.min(d)])
}
uw <- dstats(rep(1, nrow(p))); cw <- dstats(p$cap_mw)
suppressWarnings({
  ks_two <- ks.test(p$pvout[post], p$pvout[pre])
  ks_gt  <- ks.test(p$pvout[post], p$pvout[pre], alternative = "greater")  # post CDF above pre
  ks_lt  <- ks.test(p$pvout[post], p$pvout[pre], alternative = "less")     # post CDF below pre
})
share_below <- function(v, w) c(pre = sum(w[pre & p$pvout < v]) / sum(w[pre]),
                                post = sum(w[post & p$pvout < v]) / sum(w[post])) * 100
m <- 1470.097

out <- data.frame(
  quantity = c("plants pre", "plants post", "capacity pre GW", "capacity post GW",
               "mean PVOUT pre (unweighted)", "mean PVOUT post (unweighted)",
               "mean PVOUT pre (capacity-weighted)", "mean PVOUT post (capacity-weighted)",
               "unweighted D_below", "unweighted D_below at", "unweighted p (one-sided, below)",
               "unweighted D_above", "unweighted p (one-sided, above)", "unweighted D two-sided", "unweighted p two-sided",
               "capacity-weighted D_below", "capacity-weighted D_below at",
               "capacity-weighted D_above", "capacity-weighted D_above at",
               "capacity share below national median pre %", "capacity share below national median post %"),
  value = c(sum(pre), sum(post), sum(p$cap_mw[pre]) / 1e3, sum(p$cap_mw[post]) / 1e3,
            mean(p$pvout[pre]), mean(p$pvout[post]),
            weighted.mean(p$pvout[pre], p$cap_mw[pre]), weighted.mean(p$pvout[post], p$cap_mw[post]),
            uw[["D_below"]], uw[["x_below"]], ks_gt$p.value,
            uw[["D_above"]], ks_lt$p.value, unname(ks_two$statistic), ks_two$p.value,
            cw[["D_below"]], cw[["x_below"]], cw[["D_above"]], cw[["x_above"]],
            share_below(m, p$cap_mw)[["pre"]], share_below(m, p$cap_mw)[["post"]]))
out$value <- signif(out$value, 4)
write.csv(out, "outputs/analysis/si_distribution_tests.csv", row.names = FALSE)
print(out, row.names = FALSE)
