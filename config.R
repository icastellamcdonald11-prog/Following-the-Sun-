# config.R: repository paths. Every script sources this file first, so run
# scripts from the repository root, e.g.  Rscript R/paper/make_figures.R
#
# FTS_PROCESSED points to the processed input layers (see data/README.md).
# Set it with an environment variable if they live elsewhere, e.g.
#   FTS_PROCESSED=/path/to/processed Rscript R/paper/make_figures.R
FTS_ROOT      <- normalizePath(Sys.getenv("FTS_ROOT", "."), mustWork = TRUE)
FTS_PROCESSED <- Sys.getenv("FTS_PROCESSED", file.path(FTS_ROOT, "data", "processed"))
FTS_P3_ROOT   <- dirname(dirname(FTS_PROCESSED))   # parent of data/processed
setwd(FTS_ROOT)
for (d in c("outputs/figures/si", "outputs/analysis", "figures", "fig_cache",
            "robustness_outputs", "data/derived"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
