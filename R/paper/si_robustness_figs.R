source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# si_robustness_figs.R: Supplementary Notes 17-18 figures, redrawn in the
# manuscript style from saved outputs of the thesis robustness scripts
# (no analysis recomputed).
#   robustness_outputs/fig_R1_quarterly_cv.csv          (12_variance_decomposition_figs.R)
#   robustness_outputs/fig712_regression_data.csv       (11_siting_efficiency_regression.R)
#   robustness_outputs/province_fai_metrics.csv         (14b_land_opportunity_cost_FAI.R)
#   robustness_outputs/province_pvout_pop_correlation.csv (14_land_opportunity_cost.R)
# Run from the project root:  Rscript paper/analysis/si_robustness_figs.R
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(sf); library(scales) })
invisible(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))

P3  <- FTS_PROCESSED
OUT <- "outputs/figures/si"
OI_BLUE <- "#0072B2"; OI_ORANGE <- "#E69F00"; GREY <- "#7A7A7A"
PROV_EN <- c(Neimenggu = "Inner Mongolia", Xizang = "Tibet")
en <- function(x) ifelse(x %in% names(PROV_EN), PROV_EN[x], x)
th <- theme_classic(base_size = 7, base_family = "Helvetica") +
  theme(axis.title = element_text(size = 7.5), legend.position = "bottom",
        legend.title = element_blank())
save <- function(p, f, h) {
  ggsave(file.path(OUT, f), p, device = grDevices::pdf, width = 183, height = h, units = "mm")
  cat("wrote", f, "\n")
}

# --- CV of capacity-weighted PVOUT of new capacity, by quarter ---------------
cv <- read.csv("robustness_outputs/fig_R1_quarterly_cv.csv") |>
  filter(year_quarter != "2017Q4") |>
  mutate(t = build_year + (build_quarter - 1) / 4,
         period = ifelse(period == "Post Q2 2020", "From Q2 2020", "Before Q2 2020"))
p1 <- ggplot(cv, aes(t, cv, colour = period, shape = period)) +
  geom_vline(xintercept = 2020.25, linetype = "dotted", colour = GREY) +
  geom_line(linewidth = 0.4) + geom_point(size = 1.4) +
  scale_colour_manual(values = c("Before Q2 2020" = OI_BLUE, "From Q2 2020" = OI_ORANGE)) +
  scale_x_continuous(breaks = 2018:2024) +
  labs(x = "Quarter of commissioning",
       y = "Coefficient of variation of PVOUT\n(capacity-weighted, new capacity)") + th
save(p1, "figS_cv_quarterly.pdf", 70)

# --- FAI intensity vs change in within-province rank -------------------------
prov_cn <- st_drop_geometry(readRDS(file.path(P3, "china_basemap_bundle.rds"))$china_provinces)
names(prov_cn)[1] <- "cn"
rk  <- read.csv("robustness_outputs/fig712_regression_data.csv", encoding = "UTF-8")
rk$prov_py <- prov_cn$name[match(rk$prov, prov_cn$cn)]
fai <- read.csv("robustness_outputs/province_fai_metrics.csv")
d <- inner_join(rk, fai, by = c("prov_py" = "prov_en")) |>
  mutate(cap = cap_gw_post + cap_gw_pre, name = en(prov_py),
         dir = ifelse(rank_change >= 0, "Rank rose after Q2 2020", "Rank fell after Q2 2020"))
stopifnot(nrow(d) == 31)
p2 <- ggplot(d, aes(log_fai_per_km2, rank_change)) +
  geom_hline(yintercept = 0, colour = GREY, linewidth = 0.3) +
  geom_smooth(method = "lm", mapping = aes(weight = cap), formula = y ~ x,
              colour = "black", fill = "grey85", linewidth = 0.4) +
  geom_point(aes(size = cap, fill = dir), shape = 21, colour = "white", stroke = 0.2) +
  ggrepel::geom_text_repel(aes(label = name), size = 2.3, seed = 3, segment.size = 0.2,
                           max.overlaps = 40) +
  scale_fill_manual(values = c("Rank rose after Q2 2020" = OI_BLUE,
                               "Rank fell after Q2 2020" = OI_ORANGE)) +
  scale_size_area(max_size = 7, guide = "none") +
  guides(fill = guide_legend(override.aes = list(size = 3))) +
  labs(x = expression("log fixed-asset investment per km"^2*" (100 million CNY km"^-2*", 2018–2020 mean)"),
       y = "Change in capacity-weighted mean\nwithin-province rank (percentile points)") + th
save(p2, "figS_fai_rank.pdf", 95)

# --- PVOUT-population correlation by province --------------------------------
pc <- read.csv("robustness_outputs/province_pvout_pop_correlation.csv") |>
  mutate(name = en(province_en)) |> arrange(cor_pearson)
pc$name <- factor(pc$name, levels = pc$name)
p3 <- ggplot(pc, aes(cor_pearson, name, fill = cor_pearson < 0)) +
  geom_col(width = 0.7) + geom_vline(xintercept = 0, colour = GREY, linewidth = 0.3) +
  scale_fill_manual(values = c(`TRUE` = OI_BLUE, `FALSE` = OI_ORANGE),
                    labels = c(`TRUE` = "Higher PVOUT on less populated land",
                               `FALSE` = "Higher PVOUT on more populated land")) +
  labs(x = "Pearson correlation between PVOUT and log(1 + population density), feasible cells",
       y = NULL) + th
save(p3, "figS_pvout_pop_cor.pdf", 110)
cat("negative correlations:", sum(pc$cor_pearson < 0), "of", nrow(pc), "\n")
cat("capacity-weighted r(log FAI, rank change):",
    round(cov.wt(cbind(d$log_fai_per_km2, d$rank_change), wt = d$cap, cor = TRUE)$cor[1, 2], 3), "\n")
