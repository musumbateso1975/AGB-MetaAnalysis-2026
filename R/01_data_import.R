if (!exists("SENSOR_COLORS")) source("R/00_setup.R")
source("R/00c_export_helpers.R")
set.seed(2024)

cat("=== 01_data_import.R: reading the databases ===\n\n")

required_files <- c(
  "data/DATA_FisherZ.csv",
  "data/DATA_RMS.csv",
  "data/AGB_master_database.csv"
)

for (f in required_files) {
  if (!file.exists(f)) {
    stop(paste0(
      "\n❌ Fichier manquant : '", f, "'\n",
      "   Assurez-vous d'avoir :\n",
      "   1. Ouvert AGB_MetaAnalysis.Rproj (double-clic)\n",
      "   2. Ou fait : setwd('chemin/vers/AGB_MetaAnalysis_R')\n",
      "   Working directory: ", getwd()
    ))
  }
}
cat("✓ All files found in:", getwd(), "\n\n")

cat("--- Chargement DATA_FisherZ.csv ---\n")

dat_fz_raw <- readr::read_csv(
  "data/DATA_FisherZ.csv",
  col_types = readr::cols(.default = "c"), na = c("", "NA")
)
dat_fz_raw$Reference <- sub("^(.*? [0-9]{4}[a-z]?)\\. .*$", "\\1", dat_fz_raw$Reference)

CLUS_MAP <- c(S1 = "C_AfriSAR_Betul", S58 = "C_AfriSAR_Betul", S5 = "C_AfriSAR_Betul", S57 = "C_AfriSAR_Betul",
              S6 = "C_Haldwani", S12 = "C_Haldwani", S41 = "C_Sonoma", S54 = "C_Sonoma", S22 = "C_Conghua", S52 = "C_Conghua",
              S24 = "C_Maine", S32 = "C_Maine")
N_CAP <- 1000

removed_summary <- tibble::tribble(
  ~Dataset, ~Reason, ~n_obs, ~Studies,
  "Fisher_Z", "Non-comparable metric or estimand", 25L, "S48, S54, S68, S71, S72",
  "Fisher_Z", "Predictive R² that can be negative", 9L, "S63",
  "Fisher_Z", "MODIS estimates (250–1000 m) evaluated against aggregated LiDAR-derived AGB", 8L, "S67",
  "Fisher_Z", "Partial duplicate of models reported elsewhere in the same article", 4L, "S67",
  "Fisher_Z", "Accuracy of the existing GEDI L4A product", 2L, "S2",
  "Fisher_Z", "Basal area or volume, not biomass", 2L, "S27",
  "Fisher_Z", "Residual variance reported as R²", 2L, "S59",
  "Fisher_Z", "Training-fit R² (excluded from the primary analysis; sensitivity analysis only)", 10L, "S2, S6, S17, S31",
  "Fisher_Z", "Subset nested in an aggregate retained from the same study", 7L, "S40, S50, S60",
  "Fisher_Z", "Computed against GEDI L4A, not field data", 1L, "S4",
  "RMSE", "Computed against GEDI L4B, not field data", 12L, "S71, S72",
  "RMSE", "Carbon units (t C/ha)", 1L, "S61",
  "RMSE", "Accuracy of an existing product", 5L, "S2, S63",
  "RMSE", "MODIS estimates (250–1000 m) evaluated against aggregated LiDAR-derived AGB", 4L, "S67",
  "RMSE", "Partial duplicate of models reported elsewhere in the same article", 4L, "S67",
  "RMSE", "Computed against GEDI L4A, not field data", 3L, "S68",
  "RMSE", "Contradictory values in the article", 3L, "S38",
  "RMSE", "Implausible value (< 2 Mg/ha) not confirmed in the article", 2L, "S9",
  "RMSE", "Unit not convertible to Mg/ha (kg per 100 m²)", 2L, "S17",
  "RMSE", "Variable other than AGB", 2L, "S27",
  "RMSE", "AGB change, not AGB", 1L, "S62",
  "RMSE", "Training-fit RMSE (excluded from the primary analysis; sensitivity analysis only)", 2L, "S2",
  "RMSE", "Subset nested in an aggregate retained from the same study", 7L, "S40, S50, S60",
  "RMSE", "Computed against GEDI L4A, not field data", 1L, "S4"
)
removed_summary$k_studies <- lengths(strsplit(removed_summary$Studies, ", "))
removed_by_study <- tidyr::separate_rows(removed_summary, Studies, sep = ", ") |> dplyr::rename(Study_ID = Studies)

dat_fz_all0 <- dat_fz_raw |>
  dplyr::mutate(dplyr::across(c(Fisher_Z, R2, r, Training_samples, Spatial_resolution_m, Year), ~ suppressWarnings(as.numeric(.x))))
dat_fz_all <- dat_fz_all0[primary_r2(dat_fz_all0), ]

dat_fz <- dat_fz_all |>
  dplyr::filter(!is.na(Fisher_Z)) |>
  dplyr::mutate(
    Sensor_combination = factor(Sensor_combination,
      levels = c("LiDAR+Optical","LiDAR-only","Optical+SAR",
                 "Optical-only","SAR+LiDAR","Optical+SAR+LiDAR","SAR-only")),
    Model_category = factor(Model_category,
      levels = c("Statistical","Machine-learning","Deep-learning")),
    Climate_zone = factor(Climate_zone,
      levels = c("Mediterranean","Subtropical","Temperate",
                 "Arid","Tropical","Boreal")),
    GEDI_used  = factor(GEDI_used,  levels = c("Yes","No")),
    Study_scale = factor(Study_scale,
      levels = c("Local","Regional","National","Multi-country")),
    Cluster = dplyr::if_else(Study_ID %in% names(CLUS_MAP), CLUS_MAP[Study_ID], Study_ID),
    n_c = pmin(Training_samples, N_CAP),
    vi  = dplyr::if_else(!is.na(n_c) & n_c > 3, 1 / (n_c - 3), NA_real_),
    log_Resolution = log(Spatial_resolution_m + 1),
    log_Training_n = log(Training_samples + 1),
    r_bt  = tanh(Fisher_Z),
    R2_bt = tanh(Fisher_Z)^2
  )

median_nc <- median(dat_fz$n_c[dat_fz$n_c > 3], na.rm = TRUE)
median_vi <- 1 / (median_nc - 3)
dat_fz <- dat_fz |>
  dplyr::mutate(
    vi_imp     = dplyr::if_else(is.na(vi), median_vi, vi),
    vi_imputed = is.na(vi)
  )

cat("✓ DATA_FisherZ.csv loaded\n")
cat("  n =", nrow(dat_fz),
    "| Études =", length(unique(dat_fz$Study_ID)), "\n")
cat("  Mean Z =", round(mean(dat_fz$Fisher_Z), 3),
    "| Median Z =", round(median(dat_fz$Fisher_Z), 3), "\n")
cat("  vi manquants =", sum(dat_fz$vi_imputed),
    "(imputed: median capped n =", median_nc, ")\n")
cat("  Removed as non-comparable =", sum(removed_summary$n_obs[removed_summary$Dataset == "Fisher_Z"]), "obs.\n\n")

cat("--- Chargement DATA_RMS.csv ---\n")

dat_rmse_raw <- readr::read_csv(
  "data/DATA_RMS.csv",
  col_types = readr::cols(.default = "c"), na = c("", "NA")
)
dat_rmse_raw$Reference <- sub("^(.*? [0-9]{4}[a-z]?)\\. .*$", "\\1", dat_rmse_raw$Reference)

dat_rmse_all0 <- dat_rmse_raw |>
  dplyr::mutate(dplyr::across(c(RMSE, Training_samples, Spatial_resolution_m, Year), ~ suppressWarnings(as.numeric(.x))))
dat_rmse_all <- dat_rmse_all0[primary_rmse(dat_rmse_all0), ]

dat_rmse <- dat_rmse_all |>
  dplyr::filter(RMSE_unit == "Mg_ha") |>
  dplyr::mutate(
    Sensor_combination = factor(Sensor_combination,
      levels = c("LiDAR+Optical","LiDAR-only","Optical+SAR",
                 "Optical-only","SAR+LiDAR","Optical+SAR+LiDAR","SAR-only")),
    Model_category = factor(Model_category,
      levels = c("Statistical","Machine-learning","Deep-learning")),
    Climate_zone   = factor(Climate_zone,
      levels = c("Mediterranean","Subtropical","Temperate",
                 "Arid","Tropical","Boreal")),
    GEDI_used      = factor(GEDI_used, levels = c("Yes","No")),
    Study_scale    = factor(Study_scale,
      levels = c("Local","Regional","National","Multi-country")),
    Cluster        = dplyr::if_else(Study_ID %in% names(CLUS_MAP), CLUS_MAP[Study_ID], Study_ID),
    log_RMSE       = log(RMSE)
  ) |>
  dplyr::filter(RMSE > 0)

cat("✓ DATA_RMS.csv loaded\n")
cat("  n =", nrow(dat_rmse),
    "| Études =", length(unique(dat_rmse$Study_ID)), "\n")
cat("  Median RMSE =", round(median(dat_rmse$RMSE), 1), "Mg/ha\n")
cat("  Mean RMSE   =", round(mean(dat_rmse$RMSE), 1), "Mg/ha\n\n")

cat("--- Reading AGB_master_database.csv ---\n")

dat <- readr::read_csv(
  "data/AGB_master_database.csv",
  show_col_types = FALSE
)
dat$Reference <- sub("^(.*? [0-9]{4}[a-z]?)\\. .*$", "\\1", dat$Reference)

cat("✓ AGB_master_database.csv loaded\n")
cat("  n total =", nrow(dat),
    "| Études =", length(unique(dat$Study_ID)), "\n\n")

cat("--- Cross-checks ---\n")

studies_fz   <- unique(dat_fz$Study_ID)
studies_rmse <- unique(dat_rmse$Study_ID)
both         <- intersect(studies_fz, studies_rmse)

cat("Études FisherZ uniquement :", length(setdiff(studies_fz, studies_rmse)), "\n")
cat("Études RMSE uniquement    :", length(setdiff(studies_rmse, studies_fz)), "\n")
cat("Études dans les deux bases:", length(both), "\n\n")

cat("GEDI dans FisherZ  :", table(dat_fz$GEDI_used)["Yes"], "/ ",
    nrow(dat_fz), "(", round(mean(dat_fz$GEDI_used == "Yes")*100, 1), "%)\n")
cat("GEDI dans RMSE     :", table(dat_rmse$GEDI_used)["Yes"], "/ ",
    nrow(dat_rmse), "(", round(mean(dat_rmse$GEDI_used == "Yes")*100, 1), "%)\n\n")

cat("=== Data loaded ===\n")
cat("Objets disponibles : dat_fz, dat_rmse, dat\n")
cat("  dat_fz  → ", nrow(dat_fz),  "obs. Fisher's Z (", length(unique(dat_fz$Study_ID)),  "studies)\n")
cat("  dat_rmse →", nrow(dat_rmse), "obs. RMSE      (", length(unique(dat_rmse$Study_ID)), "studies)\n")
cat("  dat     → ", nrow(dat),     "obs. (master database)\n")

dat_rmse$Spatial_resolution_m <- as.numeric(dat_rmse$Spatial_resolution_m)
dat_fz$Spatial_resolution_m   <- as.numeric(dat_fz$Spatial_resolution_m)
dat_rmse$Training_samples     <- as.numeric(dat_rmse$Training_samples)
dat_fz$Training_samples       <- as.numeric(dat_fz$Training_samples)
