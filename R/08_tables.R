if (!exists("dat_fz")) { source("R/00_setup.R"); source("R/01_data_import.R") }
source("R/00b_figure_style.R"); source("R/00c_export_helpers.R")
set.seed(2024)
cat("=== 08_tables.R — Article and supplementary tables ===\n\n")

ML <- if (file.exists(file.path(CACHE_DIR, "reanalyse.rds"))) readRDS(file.path(CACHE_DIR, "reanalyse.rds")) else NULL
nrem_fz <- sum(removed_summary$n_obs[removed_summary$Dataset == "Fisher_Z"])
nrem_rm <- sum(removed_summary$n_obs[removed_summary$Dataset == "RMSE"])
fmt_ci <- function(x, lo, hi, d = 2) sprintf(paste0("%.", d, "f (%.", d, "f–%.", d, "f)"), x, lo, hi)

P1 <- if (!is.null(ML)) ML$P1 else NULL
t1 <- data.frame(
  Characteristic = c("Observations extracted", "Observations not analysed (non-comparable metric, estimand, reference, unit or scale; training fit)",
                     "Observations analysed (n)", "Studies analysed (k)", "Independent clusters",
                     "Observations using GEDI", "Publication years", "Countries",
                     "Mean Fisher's Z (SD)", "Median Fisher's Z (IQR)",
                     "Pooled R² (95% CI) — 3-level model, CR2", "95% prediction interval, untransformed R² (primary)", "95% prediction interval, Fisher's Z back-transformed (complementary)", "I² total (%), Fisher's Z",
                     "Median RMSE (IQR), Mg ha⁻¹", "Observations with imputed n"),
  FishersZ_dataset = c(nrow(dat_fz_all) + nrem_fz, nrem_fz, nrow(dat_fz), dplyr::n_distinct(dat_fz$Study_ID),
                       dplyr::n_distinct(dat_fz$Cluster),
                       sprintf("%d (%.1f%%)", sum(dat_fz$GEDI_used == "Yes", na.rm = TRUE), 100 * mean(dat_fz$GEDI_used == "Yes", na.rm = TRUE)),
                       paste0(min(dat_fz$Year, na.rm = TRUE), "–", max(dat_fz$Year, na.rm = TRUE)),
                       dplyr::n_distinct(dat_fz$Country[dat_fz$Country != "Multiple"]),
                       sprintf("%.3f (%.3f)", mean(dat_fz$Fisher_Z), sd(dat_fz$Fisher_Z)),
                       sprintf("%.3f (%.3f–%.3f)", median(dat_fz$Fisher_Z), quantile(dat_fz$Fisher_Z, .25), quantile(dat_fz$Fisher_Z, .75)),
                       if (!is.null(P1)) fmt_ci(P1$R2, P1$R2_IC_inf, P1$R2_IC_sup) else "see Table S10",
                       if (!is.null(P1)) sprintf("%.2f–%.2f", P1$R2_IPb_inf, P1$R2_IPb_sup) else "—",
                       if (!is.null(P1)) sprintf("%.2f–%.2f", P1$R2_IP_inf, P1$R2_IP_sup) else "—",
                       if (!is.null(P1)) sprintf("%.1f", P1$I2_total) else "—",
                       "—", sum(dat_fz$vi_imputed)),
  RMSE_dataset = c(nrow(dat_rmse_all) + nrem_rm, nrem_rm + nrow(dat_rmse_all) - nrow(dat_rmse), nrow(dat_rmse), dplyr::n_distinct(dat_rmse$Study_ID),
                   dplyr::n_distinct(dat_rmse$Cluster),
                   sprintf("%d (%.1f%%)", sum(dat_rmse$GEDI_used == "Yes", na.rm = TRUE), 100 * mean(dat_rmse$GEDI_used == "Yes", na.rm = TRUE)),
                   paste0(min(dat_rmse$Year, na.rm = TRUE), "–", max(dat_rmse$Year, na.rm = TRUE)),
                   dplyr::n_distinct(dat_rmse$Country[dat_rmse$Country != "Multiple"]),
                   "—", "—", "—", "—", "—", "—",
                   sprintf("%.1f (%.1f–%.1f)", median(dat_rmse$RMSE), quantile(dat_rmse$RMSE, .25), quantile(dat_rmse$RMSE, .75)), "—"))
export_table(t1, "Table1_datasets", translate = FALSE)

MODS <- c(Sensor_combination = "Sensor combination", Model_category = "Model category", GEDI_used = "GEDI used",
          Climate_zone = "Climate zone", Study_scale = "Study scale", Forest_type = "Forest type", Data_fusion = "Data fusion")
t2 <- dplyr::bind_rows(lapply(names(MODS), function(g) {
  a <- dat_fz |> dplyr::filter(!is.na(.data[[g]]), .data[[g]] != "") |> dplyr::group_by(Category = as.character(.data[[g]])) |>
    dplyr::summarise(n_obs_FisherZ = dplyr::n(), k_studies_FisherZ = dplyr::n_distinct(Study_ID), .groups = "drop")
  b <- dat_rmse |> dplyr::filter(!is.na(.data[[g]]), .data[[g]] != "") |> dplyr::group_by(Category = as.character(.data[[g]])) |>
    dplyr::summarise(n_obs_RMSE = dplyr::n(), k_studies_RMSE = dplyr::n_distinct(Study_ID), .groups = "drop")
  dplyr::full_join(a, b, by = "Category") |> dplyr::mutate(Moderator = MODS[[g]], .before = 1)
}))
export_table(t2, "TableS08_moderator_categories", translate = FALSE)

wald_z <- function(d, g) {
  d <- d[!is.na(d[[g]]) & d[[g]] != "", ]; d[[g]] <- droplevels(factor(d[[g]]))
  out <- data.frame(F = NA_real_, df_num = NA_real_, df_den = NA_real_, p = NA_real_)
  if (nlevels(d[[g]]) < 2) return(out)
  tryCatch({
    m <- metafor::rma.mv(Fisher_Z, vi_imp, mods = stats::as.formula(paste("~", g)),
                         random = ~ 1 | Cluster / Study_ID / Obs_ID, data = d, method = "REML", sparse = TRUE)
    w <- clubSandwich::Wald_test(m, constraints = clubSandwich::constrain_zero(2:m$p), vcov = "CR2", cluster = d$Cluster, test = "HTZ")
    data.frame(F = w$Fstat, df_num = w$df_num, df_den = w$df_denom, p = w$p_val)
  }, error = function(e) out)
}
wald_r <- function(d, g) {
  d <- d[!is.na(d[[g]]) & d[[g]] != "" & is.finite(d$log_RMSE), ]; d[[g]] <- droplevels(factor(d[[g]]))
  out <- data.frame(F = NA_real_, df_num = NA_real_, df_den = NA_real_, p = NA_real_)
  if (nlevels(d[[g]]) < 2) return(out)
  tryCatch({
    m <- nlme::lme(stats::as.formula(paste("log_RMSE ~", g)), random = ~ 1 | Study_ID, data = d, method = "REML")
    w <- clubSandwich::Wald_test(m, constraints = clubSandwich::constrain_zero(2:length(nlme::fixef(m))), vcov = "CR2", test = "HTZ")
    data.frame(F = w$Fstat, df_num = w$df_num, df_den = w$df_denom, p = w$p_val)
  }, error = function(e) out)
}
MOD6 <- c("Sensor_combination", "Model_category", "GEDI_used", "Climate_zone", "Study_scale", "Forest_type")
t3 <- dplyr::bind_rows(lapply(MOD6, function(g) {
  dplyr::bind_rows(
    cbind(data.frame(Moderator = MODS[[g]], Metric = "Fisher's Z", n_obs = sum(!is.na(dat_fz[[g]])),
                     k_studies = dplyr::n_distinct(dat_fz$Study_ID[!is.na(dat_fz[[g]])]),
                     levels = dplyr::n_distinct(dat_fz[[g]][!is.na(dat_fz[[g]])])), wald_z(dat_fz, g)),
    cbind(data.frame(Moderator = MODS[[g]], Metric = "log RMSE (Mg/ha)", n_obs = sum(!is.na(dat_rmse[[g]])),
                     k_studies = dplyr::n_distinct(dat_rmse$Study_ID[!is.na(dat_rmse[[g]])]),
                     levels = dplyr::n_distinct(dat_rmse[[g]][!is.na(dat_rmse[[g]])])), wald_r(dat_rmse, g)))
}))
bad <- !is.finite(t3$F) | t3$F > 1e4 | !is.finite(t3$df_den) | t3$df_den < 2
t3$Note <- ifelse(bad, "Not estimable: too few studies per level for CR2 inference", ifelse(t3$df_den < 4, "Few clusters: interpret with caution", ""))
t3[bad, c("F", "df_num", "df_den", "p")] <- NA
t3$Test <- ifelse(t3$Metric == "Fisher's Z", "Wald HTZ, 3-level rma.mv, CR2 (clusters)", "Wald HTZ, LMM (study random effect), CR2")
export_table(t3, "TableS20_moderator_tests", translate = FALSE)

DESC <- c(Obs_ID = "Unique observation identifier", Study_ID = "Study identifier (one publication or study record)",
  Site_Name = "Study site", Reference = "Short reference (author, year)", Year = "Publication year", Journal = "Journal",
  Country = "Country of study site (Multiple = multi-country)", Region = "World region", Climate_zone = "Climate zone",
  Forest_type = "Forest type", Terrain_type = "Terrain type", GEDI_used = "GEDI data used in the model (verified in full text)",
  GEDI_product = "GEDI product", Spatial_resolution_m = "Spatial resolution of predictors (m)", Resolution_class = "Resolution class (from Spatial_resolution_m)",
  Sensor_combination = "Sensor combination (7 classes)", Data_fusion = "Multi-source data fusion", Model_category = "Statistical / Machine-learning / Deep-learning",
  Model_algorithm = "Algorithm", Best_model = "Best model reported by the study", Training_samples = "Training sample size n",
  Study_scale = "Local / Regional / National / Multi-country", R2 = "Reported coefficient of determination", r = "sqrt(R²)",
  Fisher_Z = "Fisher's Z = atanh(r)", vi = "Sampling variance 1/(n − 3) as extracted", AGB_unit = "Unit of RMSE", RMSE = "Root mean square error",
  r_check = "Check: rounded r", R2_check = "Check: r_check²",
  RMSE_unit = "Unit of the analysed RMSE (Mg ha-1)",
  GEDI_coding_doubtful = "GEDI used as reference or evaluated product (not as predictor)",
  GEDI_role_unverified = "GEDI role not verified in full text", LiDAR_platform = "LiDAR platform (GEDI, ALS, UAV-LiDAR, ICESat-2…)",
  Reference_type = "Reference used to compute R²/RMSE (field plots, LiDAR-derived AGB, GEDI)",
  GEDI_role = "Role of GEDI (predictor, reference, evaluated product, not applicable)",
  Validation_design_text = "Validation scheme as described in the full text (coded into validation design by R/00d)",

  n_pseudoreplicated = "Sample size counts pixels or footprints rather than plots", R2_estimand_different = "R² not computed as 1 - SSE/SST or squared correlation",
  Map_validation = "Accuracy computed against a validation map", Simple_regression_R2 = "R² from a simple linear regression (class A)")
dict_one <- function(d, nm) data.frame(Variable = names(d), Dataset = nm,
  Type = vapply(d, function(x) { x2 <- suppressWarnings(as.numeric(x)); if (all(is.na(x2[!is.na(x) & x != ""]))) "text" else "numeric" }, ""),
  n_non_missing = vapply(d, function(x) sum(!is.na(x) & x != ""), 1L),
  n_distinct = vapply(d, function(x) dplyr::n_distinct(x[!is.na(x) & x != ""]), 1L),
  Examples = vapply(d, function(x) paste(utils::head(unique(x[!is.na(x) & x != ""]), 5), collapse = " | "), ""))
dd <- dplyr::bind_rows(dict_one(dat_fz_raw, "DATA_FisherZ.csv"), dict_one(dat_rmse_raw, "DATA_RMS.csv")) |>
  dplyr::group_by(Variable) |> dplyr::summarise(Datasets = paste(unique(Dataset), collapse = " + "), Type = dplyr::first(Type),
    n_non_missing = paste(n_non_missing, collapse = " / "), n_distinct = max(n_distinct), Examples = dplyr::first(Examples), .groups = "drop") |>
  dplyr::mutate(Description = dplyr::coalesce(unname(DESC[Variable]), "")) |>
  dplyr::arrange(match(Variable, names(dat_fz_raw)))
export_table(dd, "TableS02_data_dictionary", translate = FALSE)

st <- dplyr::bind_rows(dat_fz_all0 |> dplyr::mutate(src = "FZ"), dat_rmse_all0 |> dplyr::mutate(src = "RMSE")) |>
  dplyr::group_by(Study_ID) |>
  dplyr::summarise(Reference = dplyr::first(Reference), Year = dplyr::first(Year), Journal = dplyr::first(Journal),
                   Country = paste(unique(Country), collapse = "; "), Climate_zone = paste(unique(stats::na.omit(Climate_zone)), collapse = "; "),
                   Sensor_combinations = paste(sort(unique(stats::na.omit(Sensor_combination))), collapse = "; "),
                   Model_categories = paste(sort(unique(stats::na.omit(Model_category))), collapse = "; "),
                   GEDI_used = paste(sort(unique(stats::na.omit(GEDI_used))), collapse = "/"), .groups = "drop") |>
  dplyr::left_join(dat_fz |> dplyr::count(Study_ID, name = "n_R2"), by = "Study_ID") |>
  dplyr::left_join(dat_rmse |> dplyr::count(Study_ID, name = "n_RMSE"), by = "Study_ID") |>
  dplyr::left_join(removed_by_study |> dplyr::filter(Dataset == "Fisher_Z") |> dplyr::group_by(Study_ID) |>
                     dplyr::summarise(R2_values_removed = paste(unique(Reason), collapse = "; "), .groups = "drop"),
                   by = "Study_ID") |>
  dplyr::mutate(dplyr::across(dplyr::any_of(c("n_R2", "n_RMSE")), ~ dplyr::coalesce(.x, 0L)),
                sid = suppressWarnings(as.numeric(gsub("\\D", "", Study_ID)))) |>
  dplyr::arrange(sid, Study_ID) |> dplyr::select(-sid)
export_table(st, "TableS03_included_studies", translate = FALSE)

rem <- dplyr::select(removed_summary, Dataset, Reason, n_obs, k_studies, Studies)
export_table(rem, "TableS05_removed_observations", translate = FALSE)

files <- sort(list.files(TAB_DIR, pattern = "\\.csv$", full.names = TRUE))
ord <- order(basename(files))
wb <- openxlsx::createWorkbook()
hdr <- openxlsx::createStyle(fontColour = "white", fgFill = "#1F4E79", textDecoration = "bold", halign = "center")
for (f in files[ord]) {
  sh <- substr(sub("\\.csv$", "", basename(f)), 1, 31)
  d <- readr::read_csv(f, show_col_types = FALSE, col_types = readr::cols(.default = "c"))
  d[] <- lapply(d, function(x) { y <- suppressWarnings(as.numeric(x)); if (all(is.na(y) == is.na(x))) y else x })
  openxlsx::addWorksheet(wb, sh); openxlsx::writeData(wb, sh, d, headerStyle = hdr)
  openxlsx::setColWidths(wb, sh, cols = seq_len(ncol(d)), widths = "auto"); openxlsx::freezePane(wb, sh, firstRow = TRUE)
}
openxlsx::saveWorkbook(wb, file.path(TAB_DIR, "All_tables.xlsx"), overwrite = TRUE)
cat("✓ All_tables.xlsx :", length(files), "sheets\n=== 08_tables.R complete ===\n")
