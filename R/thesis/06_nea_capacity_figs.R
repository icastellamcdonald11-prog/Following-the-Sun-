source("config.R")  # repository paths; run scripts from the repository root
# =============================================================================
# 06_nea_capacity_figs.R  —  Thesis Project 4
#
# Produces Figures 2.2 and 2.3:
#   2.2  Annual Solar PV Additions (UPV, DPV, Total)
#   2.3  Cumulative Solar PV Capacity (UPV, DPV, Total)
#
# Source: National Energy Administration (NEA) press releases.
# ** 2018/2019 UPV/DPV split derived from most recently reported mix.
# ** 2019/2020 cumulative UPV/DPV derived from most recently reported mix.
#
# Only change from original: Fig 2.3 now carries the same
# "China's dual-carbon goals announced" annotation as Fig 2.2.
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

P4  <- FTS_ROOT
OUT <- file.path(P4, "figures")
dir.create(OUT, showWarnings = FALSE)

# =============================================================================
# DATA
# =============================================================================
nea <- tribble(
  ~year, ~upv_added, ~dpv_added, ~total_added, ~upv_cum, ~dpv_cum, ~total_cum,
  2017,      33.62,      19.44,        53.06,   100.59,    26.66,      130,
  2018,      23.32,      20.95,        44.26,   123.84,    50.61,      174,
  2019,      17.91,      12.20,        30.11,   141.75,    62.81,      204,
  2020,      32.86,      15.52,        48.20,   172.10,    81.17,      253,
  2021,      25.60,      29.28,        54.88,   203.00,   110.45,      306,
  2022,      36.30,      51.11,        87.41,   239.20,   154.10,      393,
  2023,     120.00,      96.28,       216.30,   354.48,   254.43,   608.91,
  2024,     159.00,     120.00,       279.00,   511.00,   370.00,      881,
  2025,     164.00,     153.00,       317.00,   670.00,   530.00,     1200
)

# Xi Jinping announced dual-carbon goals at UNGA, 22 Sep 2020
DUAL_CARBON_X <- 2020 + (9 / 12)

# =============================================================================
# SHARED THEME  (matches original: classic, minimal grid, bold title)
# =============================================================================
theme_nea <- theme_classic(base_size = 11) +
  theme(
    axis.line        = element_line(colour = "grey30", linewidth = 0.35),
    axis.ticks       = element_line(colour = "grey30", linewidth = 0.35),
    axis.text        = element_text(colour = "grey20"),
    panel.grid.major.y = element_line(colour = "grey88", linewidth = 0.3),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "top",
    legend.key.width   = unit(1.5, "cm"),
    legend.key.height  = unit(0.4, "cm"),
    legend.text        = element_text(size = 9),
    legend.title       = element_blank(),
    plot.title         = element_text(face = "bold", size = 12),
    plot.subtitle      = element_text(size = 9, colour = "grey40"),
    plot.caption       = element_text(size = 7.5, colour = "grey45",
                                      hjust = 0, lineheight = 1.3)
  )

# Line types to match original
lt <- c("DPV"   = "solid",
        "Total" = "dashed",
        "UPV"   = "dotdash")

# =============================================================================
# FIG 2.2  —  Annual additions
# =============================================================================
added_long <- nea %>%
  select(year, UPV = upv_added, DPV = dpv_added, Total = total_added) %>%
  pivot_longer(-year, names_to = "segment", values_to = "gw") %>%
  mutate(segment = factor(segment, levels = c("DPV", "Total", "UPV")))

# Legend labels to match original
lab_added <- c("DPV" = "DPV added", "Total" = "Total added", "UPV" = "UPV added")

fig_2_2 <- ggplot(added_long, aes(x = year, y = gw,
                                   linetype = segment, group = segment)) +
  geom_vline(xintercept = DUAL_CARBON_X,
             linetype = "dashed", colour = "grey50", linewidth = 0.45) +
  annotate("text", x = DUAL_CARBON_X + 0.08, y = max(added_long$gw) * 0.95,
           label = "China\u2019s dual-carbon goals announced",
           hjust = 0, size = 3, colour = "grey30") +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.8, fill = "black", shape = 21) +
  scale_linetype_manual(values = lt, labels = lab_added) +
  scale_x_continuous(breaks = nea$year, expand = expansion(add = 0.3)) +
  scale_y_continuous(name   = "Added capacity (GW)",
                     limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05)),
                     breaks = seq(0, 300, 50)) +
  labs(
    title    = "Annual Solar PV Additions",
    subtitle = "Yearly added capacity by segment (NEA press releases)",
    x        = NULL,
    caption  = paste0("Source: National Energy Administration (NEA).\n",
                      "** 2018/2019 UPV/DPV split derived from most recently ",
                      "reported segment mix.")
  ) +
  theme_nea

ggsave(file.path(OUT, "fig_2_2_annual_additions.png"),
       fig_2_2, width = 8, height = 5, dpi = 300)
cat("Saved fig_2_2_annual_additions.png\n")

# =============================================================================
# FIG 2.3  —  Cumulative capacity
# =============================================================================
cum_long <- nea %>%
  select(year, UPV = upv_cum, DPV = dpv_cum, Total = total_cum) %>%
  pivot_longer(-year, names_to = "segment", values_to = "gw") %>%
  mutate(segment = factor(segment, levels = c("DPV", "Total", "UPV")))

lab_cum <- c("DPV" = "DPV cumulative", "Total" = "Total cumulative",
             "UPV" = "UPV cumulative")

fig_2_3 <- ggplot(cum_long, aes(x = year, y = gw,
                                 linetype = segment, group = segment)) +
  geom_vline(xintercept = DUAL_CARBON_X,
             linetype = "dashed", colour = "grey50", linewidth = 0.45) +
  annotate("text", x = DUAL_CARBON_X + 0.08, y = max(cum_long$gw) * 0.95,
           label = "China\u2019s dual-carbon goals announced",
           hjust = 0, size = 3, colour = "grey30") +
  geom_line(linewidth = 0.65) +
  geom_point(size = 1.8, fill = "black", shape = 21) +
  scale_linetype_manual(values = lt, labels = lab_cum) +
  scale_x_continuous(breaks = nea$year, expand = expansion(add = 0.3)) +
  scale_y_continuous(name   = "Capacity (GW)",
                     limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.05)),
                     labels = scales::comma,
                     breaks = seq(0, 1250, 250)) +
  labs(
    title    = "Cumulative Solar PV Capacity",
    subtitle = "Utility-scale (UPV), distributed (DPV), and total (NEA press releases)",
    x        = NULL,
    caption  = paste0("Source: National Energy Administration (NEA).\n",
                      "** 2019/2020 cumulative UPV/DPV derived from most recently ",
                      "reported segment mix.")
  ) +
  theme_nea

ggsave(file.path(OUT, "fig_2_3_cumulative_capacity.png"),
       fig_2_3, width = 8, height = 5, dpi = 300)
cat("Saved fig_2_3_cumulative_capacity.png\n")

# =============================================================================
# FIG 2.4  —  Implied decommissioned capacity
#
# Method: decommissioned(t) = cum(t-1) + added(t) - cum(t)
#
# Several years produce small negatives for UPV/DPV segments — these are
# artefacts of:
#   (a) ** estimated UPV/DPV splits in 2018–2020
#   (b) rounding in NEA figures (e.g. 2025 total rounds to 1200 GW)
# All negatives are floored to zero. Total decommissioned is more reliable
# than the segment breakdown and is shown as the primary series.
# =============================================================================

decomm <- nea %>%
  arrange(year) %>%
  mutate(
    total_decomm = pmax(0, lag(total_cum) + total_added - total_cum),
    upv_decomm   = pmax(0, lag(upv_cum)   + upv_added   - upv_cum),
    dpv_decomm   = pmax(0, lag(dpv_cum)   + dpv_added   - dpv_cum)
  ) %>%
  filter(!is.na(total_decomm))   # drop 2017 (no prior year)

# Flag years where any segment was floored (for caption transparency)
raw <- nea %>%
  arrange(year) %>%
  mutate(
    total_raw = lag(total_cum) + total_added - total_cum,
    upv_raw   = lag(upv_cum)   + upv_added   - upv_cum,
    dpv_raw   = lag(dpv_cum)   + dpv_added   - dpv_cum
  ) %>%
  filter(!is.na(total_raw))

n_floored <- sum(raw$total_raw < -0.5 | raw$upv_raw < -0.5 | raw$dpv_raw < -0.5,
                 na.rm = TRUE)

fig_2_4 <- ggplot(decomm, aes(x = year, y = total_decomm)) +
  geom_col(width = 0.65, fill = "#555555", colour = NA) +
  scale_x_continuous(breaks = decomm$year, expand = expansion(add = 0.4)) +
  scale_y_continuous(name   = "Implied decommissioned capacity (GW)",
                     limits = c(0, NA),
                     expand = expansion(mult = c(0, 0.08)),
                     breaks = seq(0, 10, 2)) +
  labs(
    title    = "Implied Decommissioned Solar PV Capacity",
    subtitle = "Derived as: prior-year cumulative + annual additions \u2212 end-year cumulative (NEA totals only)",
    x        = NULL,
    caption  = paste0(
      "Source: National Energy Administration (NEA); author calculations.\n",
      "Based on independently-reported NEA total cumulative figures; ",
      "segment-level (UPV/DPV) breakdown omitted due to definitional reclassifications in source data.\n",
      "Small negative implied values (rounding artefacts) treated as zero.")
  ) +
  theme_nea

ggsave(file.path(OUT, "fig_2_4_decommissioned.png"),
       fig_2_4, width = 8, height = 5, dpi = 300)
cat("Saved fig_2_4_decommissioned.png\n")

cat("\nImplied decommissioned capacity (GW):\n")
print(decomm %>% select(year, upv_decomm, dpv_decomm, total_decomm) %>%
        mutate(across(-year, ~ round(.x, 2))))
