# run_paper.R: reproduce the main-text figures and the numbers recomputed for
# the manuscript from the derived data in this repository.
#   Rscript run_paper.R
# Needs: fig_cache/, transmission/, robustness_outputs/ (included), the
# processed basemap and rasters in FTS_PROCESSED, and data/derived/*.tif
# (see data/README.md). Outputs go to outputs/.
scripts <- c("R/paper/make_figures.R",          # Figs 1-5
             "R/paper/foregone_generation.R",   # Methods eq. 15, SI Table S10
             "R/paper/deterioration_shares.R",  # Methods eq. 12, SI Table S11
             "R/paper/si_distribution_tests.R", # SI Note 17, Table S14
             "R/paper/si_robustness_figs.R",    # SI Figs S20-S22
             "R/paper/threshold_sensitivity.R") # SI Note 19, Tables S16-S17
for (s in scripts) {
  cat("\n==>", s, "\n")
  status <- system2("Rscript", s)
  if (status != 0) stop("failed: ", s)
}
cat("\nDone. See outputs/.\n")
