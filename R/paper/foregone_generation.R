source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# foregone_generation.R: monetisation of the post-2020 decline (Methods eq. 15)
#
# Run from the project root:  Rscript paper/analysis/foregone_generation.R
#
# Definition (as in the thesis: "variation across the post-2020 period, applied
# to installed capacity"), now computed from the capacity-weighted PVOUT of the
# cumulative installed fleet in annual units (the Fig. 3a series):
#   lower bound  dP_lo = PVOUT(Q2 2020)        - min over Q2 2020..Q2 2024
#   upper bound  dP_hi = max over Q2 2020..Q2 2024 - min over Q2 2020..Q2 2024
#   dE = dP x K, with K = installed capacity at Q2 2024
#   V  = dE x P;  NPV = V x sum_{t=1..25} (1+r)^-t, r = 6% and 8%
#   USD at the 2024 annual average rate, 7.1957 CNY per US$
#   (Board of Governors of the Federal Reserve System, G.5A; FRED AEXCHUS)
# Writes paper/analysis/foregone_generation.csv and ..._sensitivity.csv
# =============================================================================

plants <- readRDS("fig_cache/plants_quad.rds")
qn <- function(x) as.integer(substr(x, 1, 4)) + (as.integer(substr(x, 6, 6)) - 1) / 4
plants$t <- qn(plants$year_quarter)

qs <- sort(unique(plants$t))
series <- data.frame(t = qs, pv = sapply(qs, function(q) {
  s <- plants[plants$t <= q, ]; weighted.mean(s$pvout, s$cap_mw) }))
post <- series[series$t >= 2020.25, ]
qlab <- function(t) sprintf("%dQ%d", floor(t + 1e-9), round((t - floor(t + 1e-9)) * 4) + 1)

pv_q2_2020 <- series$pv[abs(series$t - 2020.25) < 1e-9]
pv_min <- min(post$pv); pv_max <- max(post$pv)
dP <- c(lo = pv_q2_2020 - pv_min, hi = pv_max - pv_min)          # kWh kWp-1 yr-1
K_gw <- sum(plants$cap_mw) / 1e3                                    # GW at Q2 2024
dE_twh <- dP * K_gw * 1e6 / 1e9                                     # kWh/kWp * kW -> TWh

PRICE <- 0.37          # CNY/kWh, grid-parity benchmark (Methods)
FX    <- 7.1957        # CNY per US$, 2024 annual average
annuity <- function(r, T = 25) sum((1 + r)^-(1:T))
V_cny  <- dE_twh * 1e9 * PRICE / 1e9                                # bn CNY / yr
npv_cny <- c(r8_lo = V_cny[["lo"]] * annuity(0.08), r6_hi = V_cny[["hi"]] * annuity(0.06))

out <- data.frame(
  quantity = c("PVOUT at Q2 2020", "post-2020 maximum", "post-2020 minimum",
               "dPVOUT lower", "dPVOUT upper", "installed capacity Q2 2024",
               "dE lower", "dE upper", "annual value lower", "annual value upper",
               "NPV lower (8%, lower dE)", "NPV upper (6%, upper dE)",
               "annual value lower USD", "annual value upper USD",
               "NPV lower USD", "NPV upper USD"),
  value = c(pv_q2_2020, pv_max, pv_min, dP[["lo"]], dP[["hi"]], K_gw,
            dE_twh[["lo"]], dE_twh[["hi"]], V_cny[["lo"]], V_cny[["hi"]],
            npv_cny[["r8_lo"]], npv_cny[["r6_hi"]],
            V_cny[["lo"]] / FX, V_cny[["hi"]] / FX,
            npv_cny[["r8_lo"]] / FX, npv_cny[["r6_hi"]] / FX),
  unit = c(rep("kWh kWp-1 yr-1", 5), "GW", "TWh yr-1", "TWh yr-1",
           "bn CNY yr-1", "bn CNY yr-1", "bn CNY", "bn CNY",
           "bn US$ yr-1", "bn US$ yr-1", "bn US$", "bn US$"),
  note = c("", qlab(post$t[which.max(post$pv)]), qlab(post$t[which.min(post$pv)]),
           rep("", 13)))
out$value <- round(out$value, 2)
dir.create("outputs/analysis", showWarnings = FALSE)
write.csv(out, "outputs/analysis/foregone_generation.csv", row.names = FALSE)
print(out, row.names = FALSE)

# Sensitivity (Supplementary Table S10): price scenarios x discount rates
sens <- expand.grid(price = c(0.14, 0.20, 0.26, 0.37), r = c(0.06, 0.08))
sens$annual_lo_bn_cny <- dE_twh[["lo"]] * sens$price
sens$annual_hi_bn_cny <- dE_twh[["hi"]] * sens$price
sens$npv_lo_bn_cny <- sens$annual_lo_bn_cny * sapply(sens$r, annuity)
sens$npv_hi_bn_cny <- sens$annual_hi_bn_cny * sapply(sens$r, annuity)
sens$npv_lo_bn_usd <- sens$npv_lo_bn_cny / FX
sens$npv_hi_bn_usd <- sens$npv_hi_bn_cny / FX
sens[, -(1:2)] <- round(sens[, -(1:2)], 1)
write.csv(sens, "outputs/analysis/foregone_generation_sensitivity.csv", row.names = FALSE)
print(sens, row.names = FALSE)
