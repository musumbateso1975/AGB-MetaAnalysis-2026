if (!file.exists("AGB_MetaAnalysis.Rproj")) {
  here <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NA_character_)
  if (!is.na(here) && file.exists(file.path(here, "AGB_MetaAnalysis.Rproj"))) setwd(here) else
    stop("Set the working directory to the project folder (the one containing AGB_MetaAnalysis.Rproj).")
}

required <- c("data/DATA_FisherZ.csv", "data/DATA_RMS.csv",
              "data/AGB_master_database.csv")
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing input files: ", paste(missing, collapse = ", "))

unlink(c("output/figures", "output/tables"), recursive = TRUE)

scripts <- c(
  "R/00_setup.R",
  "R/01_data_import.R",
  "R/02_descriptive_stats.R",
  "R/03_fisher_z_meta.R",
  "R/10_multilevel_reanalysis.R",
  "R/07_figures.R",
  "R/11_main_figures.R",
  "R/12_terrain_samplesize_gedi.R",
  "R/13_algorithm_table.R",
  "R/14_estimands_paired_conventional.R",
  "R/15_final_validation.R",
  "R/16_two_scales_and_pairs.R",
  "R/17_study_summary_figures.R",
  "R/18_robustness_analyses.R",
  "R/19_model_diagnostics.R",
  "R/08_tables.R"
)

errors <- character(0)
t_start <- proc.time()
for (s in scripts) {
  cat(strrep("=", 60), "\n", s, "\n", strrep("=", 60), "\n", sep = "")
  tryCatch(source(s, local = if (grepl("^R/1[0-3]_", s)) new.env() else FALSE),
           error = function(e) errors <<- c(errors, paste0(s, ": ", conditionMessage(e))))
}

cat(strrep("=", 60), "\n")
cat("Scripts completed:", length(scripts) - length(errors), "/", length(scripts),
    sprintf("(%.1f min)\n", (proc.time() - t_start)[["elapsed"]] / 60))
if (length(errors)) {
  cat("Errors:\n"); for (e in errors) cat(" -", e, "\n")
} else {
  cat("Figures: output/figures (PNG 1000 dpi and PDF)\nTables:  output/tables (CSV and All_tables.xlsx)\n")
}
