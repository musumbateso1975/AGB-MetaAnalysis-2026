if (!exists("dat_fz")) { source("R/00_setup.R"); source("R/01_data_import.R") }
set.seed(2024)

cat("=== 02_descriptive_stats.R — Descriptive Statistics ===\n\n")

cat("--- 1. Fisher's Z Dataset Summary ---\n")

fz_summary <- dat_fz |>
  dplyr::summarise(
    n_obs          = dplyr::n(),
    n_studies      = dplyr::n_distinct(Study_ID),
    mean_Z         = mean(Fisher_Z),
    sd_Z           = sd(Fisher_Z),
    median_Z       = median(Fisher_Z),
    IQR_lo         = quantile(Fisher_Z, 0.25),
    IQR_hi         = quantile(Fisher_Z, 0.75),
    r_mean         = tanh(mean_Z),
    R2_mean        = tanh(mean_Z)^2,
    n_GEDI         = sum(GEDI_used == "Yes", na.rm = TRUE),
    pct_GEDI       = mean(GEDI_used == "Yes", na.rm = TRUE) * 100
  )

print(as.data.frame(fz_summary))

cat("\nPooled mean Z =", round(fz_summary$mean_Z, 3),
    "| SD =", round(fz_summary$sd_Z, 3),
    "| r =", round(fz_summary$r_mean, 3),
    "| R² ≈", round(fz_summary$R2_mean, 2), "\n")
cat("Median Z =", round(fz_summary$median_Z, 3),
    "| IQR:", round(fz_summary$IQR_lo, 3), "–",
    round(fz_summary$IQR_hi, 3), "\n\n")

cat("--- 2. RMSE Dataset Summary ---\n")

rmse_summary <- dat_rmse |>
  dplyr::summarise(
    n_obs       = dplyr::n(),
    n_studies   = dplyr::n_distinct(Study_ID),
    median_RMSE = median(RMSE),
    mean_RMSE   = mean(RMSE),
    sd_RMSE     = sd(RMSE),
    IQR_lo      = quantile(RMSE, 0.25),
    IQR_hi      = quantile(RMSE, 0.75),
    n_GEDI      = sum(GEDI_used == "Yes", na.rm = TRUE),
    pct_GEDI    = mean(GEDI_used == "Yes", na.rm = TRUE) * 100
  )

print(as.data.frame(rmse_summary))

cat("\nMedian RMSE =", round(rmse_summary$median_RMSE, 1),
    "Mg/ha | Mean =", round(rmse_summary$mean_RMSE, 1),
    "| SD =", round(rmse_summary$sd_RMSE, 1), "\n\n")

cat("--- 3. By Sensor Combination ---\n")

sensor_fz <- dat_fz |>
  dplyr::group_by(Sensor_combination) |>
  dplyr::summarise(
    n       = dplyr::n(),
    mean_Z  = round(mean(Fisher_Z), 3),
    sd_Z    = round(sd(Fisher_Z), 3),
    r_equiv = round(tanh(mean_Z), 3),
    R2_equiv = round(tanh(mean_Z)^2, 2),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(mean_Z))

sensor_rmse <- dat_rmse |>
  dplyr::group_by(Sensor_combination) |>
  dplyr::summarise(
    n_rmse      = dplyr::n(),
    median_RMSE = round(median(RMSE), 1),
    .groups     = "drop"
  )

table_sensor <- dplyr::left_join(sensor_fz, sensor_rmse, by = "Sensor_combination")
cat("\nTable 2 — Sensor Combination Results:\n")
print(table_sensor)

cat("\n--- 4. By Model Category ---\n")

model_fz <- dat_fz |>
  dplyr::group_by(Model_category) |>
  dplyr::summarise(
    n       = dplyr::n(),
    mean_Z  = round(mean(Fisher_Z), 3),
    sd_Z    = round(sd(Fisher_Z), 3),
    r_equiv = round(tanh(mean_Z), 3),
    R2_equiv = round(tanh(mean_Z)^2, 2),
    .groups = "drop"
  )

model_rmse <- dat_rmse |>
  dplyr::group_by(Model_category) |>
  dplyr::summarise(
    n_rmse      = dplyr::n(),
    median_RMSE = round(median(RMSE), 1),
    .groups     = "drop"
  )

table_model <- dplyr::left_join(model_fz, model_rmse, by = "Model_category")
cat("\nTable 3 — Model Category Results:\n")
print(table_model)

cat("\n--- 5. By Climate Zone ---\n")

climate_fz <- dat_fz |>
  dplyr::filter(!is.na(Climate_zone)) |>
  dplyr::group_by(Climate_zone) |>
  dplyr::summarise(
    n       = dplyr::n(),
    mean_Z  = round(mean(Fisher_Z), 3),
    R2_equiv = round(tanh(mean(Fisher_Z))^2, 2),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(mean_Z))

climate_rmse <- dat_rmse |>
  dplyr::filter(!is.na(Climate_zone)) |>
  dplyr::group_by(Climate_zone) |>
  dplyr::summarise(
    n_rmse      = dplyr::n(),
    median_RMSE = round(median(RMSE), 1),
    .groups     = "drop"
  )

table_climate <- dplyr::left_join(climate_fz, climate_rmse, by = "Climate_zone")
cat("\nTable 4 — Climate Zone Results:\n")
print(table_climate)

cat("\n--- 6. GEDI vs Non-GEDI ---\n")

gedi_fz <- dat_fz |>
  dplyr::group_by(GEDI_used) |>
  dplyr::summarise(
    n       = dplyr::n(),
    mean_Z  = round(mean(Fisher_Z), 3),
    R2_equiv = round(tanh(mean(Fisher_Z))^2, 2),
    .groups = "drop"
  )

gedi_rmse <- dat_rmse |>
  dplyr::group_by(GEDI_used) |>
  dplyr::summarise(
    n_rmse      = dplyr::n(),
    median_RMSE = round(median(RMSE), 1),
    .groups     = "drop"
  )

table_gedi <- dplyr::left_join(gedi_fz, gedi_rmse, by = "GEDI_used")
cat("\nGEDI Comparison:\n")
print(table_gedi)

cat("\n--- 7. Heatmap: Sensor × Model ---\n")

heatmap_data <- dat_fz |>
  dplyr::group_by(Sensor_combination, Model_category) |>
  dplyr::summarise(
    n       = dplyr::n(),
    mean_Z  = round(mean(Fisher_Z), 3),
    R2_equiv = round(tanh(mean(Fisher_Z))^2, 2),
    .groups = "drop"
  ) |>
  dplyr::filter(n >= 2)

cat("\nHeatmap data (n >= 2):\n")
print(heatmap_data, n = Inf)

cat("\n--- 8. Geographic Distribution ---\n")

geo_summary <- dat_fz |>
  dplyr::count(Country, sort = TRUE) |>
  dplyr::mutate(pct = round(n / sum(n) * 100, 1)) |>
  utils::head(15)
print(geo_summary)

study_geo <- dat_fz |>
  dplyr::distinct(Study_ID, Country) |>
  dplyr::count(Country, sort = TRUE)
cat("\nStudies per country (top 10):\n")
print(utils::head(study_geo, 10))

table1 <- data.frame(
  Characteristic       = c("Total observations", "Total studies",
                            "GEDI-based obs.", "Publication years",
                            "Countries", "Mean Z (SD)", "Median Z (IQR)",
                            "r (R² equiv.)", "Median RMSE (Mg/ha)"),
  Fisher_Z_Dataset     = c(
    fz_summary$n_obs,
    fz_summary$n_studies,
    paste0(fz_summary$n_GEDI, " (", round(fz_summary$pct_GEDI, 1), "%)"),
    paste0(min(dat_fz$Year), "–", max(dat_fz$Year)),
    length(unique(dat_fz$Country)),
    paste0(round(fz_summary$mean_Z, 3), " (", round(fz_summary$sd_Z, 3), ")"),
    paste0(round(fz_summary$median_Z, 3), " (",
           round(fz_summary$IQR_lo, 3), "–",
           round(fz_summary$IQR_hi, 3), ")"),
    paste0(round(fz_summary$r_mean, 3), " (R² ≈ ", round(fz_summary$R2_mean, 2), ")"),
    "—"
  ),
  RMSE_Dataset         = c(
    rmse_summary$n_obs,
    rmse_summary$n_studies,
    paste0(rmse_summary$n_GEDI, " (", round(rmse_summary$pct_GEDI, 1), "%)"),
    paste0(min(dat_rmse$Year), "–", max(dat_rmse$Year)),
    length(unique(dat_rmse$Country)),
    "—",
    "—",
    "—",
    round(rmse_summary$median_RMSE, 1)
  )
)

cat("\n=== Descriptive statistics complete ===\n")
