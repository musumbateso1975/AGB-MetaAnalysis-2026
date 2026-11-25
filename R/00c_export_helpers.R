# Working sampling variance for models fitted on untransformed R2 (equal weights).
# A small constant keeps all between-observation heterogeneity in the variance components,
# so prediction intervals on the R2 scale are not narrowed by a Fisher's Z variance (mean vi ~ 0.025).
V_R2 <- 1e-4

TAB_DIR <- "output/tables"
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)

COL_EN <- c(
  Scenario = "Scenario", k_obs = "n_obs", Etudes = "k_studies", Clusters = "k_clusters",
  IC_inf = "CI_lower", IC_sup = "CI_upper", IP_inf = "PI_lower", IP_sup = "PI_upper",
  R2_IC_inf = "R2_CI_lower", R2_IC_sup = "R2_CI_upper", R2_IP_inf = "R2_PI_lower", R2_IP_sup = "R2_PI_upper", R2_IPb_inf = "R2_PI_untransformed_lower", R2_IPb_sup = "R2_PI_untransformed_upper",
  ddl = "df_Satterthwaite", sigma2_etude = "sigma2_study", sigma2_obs = "sigma2_observation",
  I2_etude = "I2_study", I2_obs = "I2_observation", v_typique = "typical_sampling_variance",
  Inference = "Inference", Terme = "Term", SE_CR2 = "SE_CR2", Moderateur = "Moderator",
  ddl_num = "df_num", ddl_den = "df_den", tau2_nul = "tau2_null", tau2_residuel = "tau2_residual",
  Variable = "Moderator", Niveau = "Level", Fiabilite = "Reliability", Comparaison = "Comparison",
  Etudes_GEDI = "k_studies_GEDI", Etudes_sans = "k_studies_without_GEDI",
  Contraste = "Contrast", Etudes_Z = "k_studies_Z", dZ_RE = "dZ_pooled", dZ_IC_inf = "dZ_CI_lower",
  dZ_IC_sup = "dZ_CI_upper", dR2_moy = "dR2_mean", dR2_IC_inf = "dR2_CI_lower", dR2_IC_sup = "dR2_CI_upper",
  Etudes_A_meilleur_R2 = "k_studies_A_higher_R2", Etudes_RMSE = "k_studies_RMSE", RMSE_var_pct = "RMSE_change_pct",
  RMSE_IC_inf_pct = "RMSE_change_CI_lower_pct", RMSE_IC_sup_pct = "RMSE_change_CI_upper_pct",
  Etudes_A_RMSE_plus_faible = "k_studies_A_lower_RMSE", kA = "n_obs_A", kB = "n_obs_B",
  Cluster_retire = "Cluster_removed", pente = "slope", Combinaison = "Sensor_combination",
  RMSE_mediane_des_etudes = "RMSE_median_of_study_medians", Contexte = "Context", Capteurs = "Sensors",
  R2_IC = "R2_95CI", R2_IP = "R2_95PI", Statut = "Status")

VAL_EN <- c(
  "Principal (plafond n = 1000, CR2)" = "Primary (n capped at 1000; CR2)",
  "Sans n > 5 000" = "Excluding n > 5,000", "Sans plafond de n" = "No cap on n", "Plafond n = 500" = "n capped at 500",
  "Cas complets (sans n imputé)" = "Complete cases (no imputed n)",
  "Une obs./étude : meilleur modèle (borne optimiste)" = "One obs./study: best model (optimistic bound)",
  "Une obs./étude tirée au hasard (500 tirages, moyenne)" = "One obs./study at random (500 draws, mean)",
  "Sans la Chine" = "Excluding China", "Sans l'Asie de l'Est" = "Excluding East Asia", "Tropical uniquement" = "Tropical only",
  "Sans validations de cartes ni GEDI suspect" = "Excluding map validations and doubtful GEDI coding",
  "Régressions linéaires simples (R² = r² défendable)" = "Simple linear regressions only (R² = r²)",
  "Afrique centrale uniquement" = "Central Africa only", "Sans rôle GEDI non vérifié" = "Excluding unverified GEDI role", "Avec les R² d'ajustement (training fit)" = "Including training-fit R²",
  "CR2 robuste" = "CR2 robust", "modèle" = "model-based", "intervalle = 2,5–97,5 % des tirages" = "2.5–97.5% of draws",
  "Mêmes modèles sur sous-ensembles" = "Same models on subsets",
  "Combinaison de capteurs" = "Sensor combination", "Catégorie de modèle" = "Model category",
  "Échelle" = "Study scale", "Zone climatique" = "Climate zone",
  "Multi-capteurs vs mono-capteur" = "Multi-sensor vs single-sensor", "Avec LiDAR vs sans LiDAR" = "With LiDAR vs without LiDAR",
  "Apprentissage automatique vs statistique" = "Machine learning vs statistical",
  "Apprentissage profond vs autres" = "Deep learning vs other models",
  "GEDI vs sans GEDI (non suspect ; descriptif)" = "GEDI vs other LiDAR (descriptive)",
  "Sans GEDI (autre LiDAR) vs GEDI, combinaisons avec LiDAR" = "Other LiDAR (no GEDI) vs GEDI, LiDAR combinations only",
  "Egger multiniveau (Z ~ SE, 3 niveaux, CR2)" = "Multilevel Egger test (Z ~ SE, 3 levels, CR2)",
  "Avec LiDAR" = "With LiDAR", "Sans LiDAR" = "Without LiDAR", "Sans GEDI" = "No GEDI", "Régional" = "Regional",
  "National/multi-pays" = "National/multi-country", "Estimable" = "Estimable", "Aucun" = "None",
  "Preuves insuffisantes (< 5 études)" = "Insufficient evidence (< 5 studies)",
  "< 5 études : indicatif" = "< 5 studies: indicative", "≥ 5 études" = "≥ 5 studies")
VAL_RX <- c("^Climat : " = "Climate: ", "^Échelle : " = "Scale: ", "Régional" = "Regional", "multi-pays" = "multi-country",
            "^ScaleRégional$" = "ScaleRegional", "études" = "studies",
            "Terrain \\(vérifié\\)" = "Field (verified)", "Terrain \\(présumé\\)" = "Field (presumed)", "AGB dérivée de LiDAR" = "LiDAR-derived AGB",
            "^Prédicteur$" = "Predictor", "^Référence probable$" = "Probable reference", "^Référence$" = "Reference", "^Simulé$" = "Simulated",
            "^Non vérifié$" = "Not verified", "Carte globale de hauteur" = "Global canopy height map", "simulé" = "simulated",
            "Hold-out aléatoire 20 % des placettes terrain" = "Random 20% hold-out of field plots",
            "Hold-out aléatoire 20 % de pixels LiDAR" = "Random 20% hold-out of LiDAR pixels", "CV spatiale 4 plis" = "Spatial 4-fold CV",
            "Validation indépendante de la carte nationale contre les placettes de l'IFN RDC 2017–2018" = "Independent validation of the national map against DRC NFI plots (2017–2018)",
            "GEDI L4A \\(vérifié\\)" = "GEDI L4A (verified)",
            "Hold-out aléatoire" = "Random hold-out")

en_table <- function(df) {
  df <- as.data.frame(df)
  df <- df[, !duplicated(names(df)), drop = FALSE]
  for (c in names(df)) if (is.character(df[[c]]) || is.factor(df[[c]])) {
    v <- as.character(df[[c]]); hit <- v %in% names(VAL_EN); v[hit] <- VAL_EN[v[hit]]
    for (rx in names(VAL_RX)) v <- gsub(rx, VAL_RX[[rx]], v)
    df[[c]] <- v
  } else if (is.numeric(df[[c]])) df[[c]] <- round(df[[c]], 4)
  nm <- names(df); hit <- nm %in% names(COL_EN); nm[hit] <- COL_EN[nm[hit]]; names(df) <- nm
  df
}

OUT_TAB <- c(TableS02_data_dictionary = "TableS02_data_dictionary", TableS03_included_studies = "TableS03_included_studies",
  TableS05_removed_observations = "TableS04_removed_values", TableS10_pooled_model = "TableS05_pooled_model",
  TableS11_estimand_classes = "TableS06_estimand_classes", TableS12_effect_size_strategies = "TableS07_effect_size_strategies",
  TableS13_scale_and_estimand = "TableS08_scale_and_estimand", TableS15_M1_M4_two_scales = "TableS09_M1_M4_two_scales",
  TableS16_moderators_M1_M4_two_scales = "TableS10_moderators_M1_M4_two_scales", TableS17_mechanisms = "TableS11_mechanisms",
  TableS19_moderator_robustness = "TableS12_moderator_robustness", TableS20_moderator_tests = "TableS13_moderator_tests",
  TableS22_moderator_confounding = "TableS14_moderator_confounding", TableS23_meta_regression = "TableS15_meta_regression",
  TableS24_meta_regression_Wald = "TableS16_meta_regression_Wald", TableS25_variance_explained = "TableS17_variance_explained",
  TableS26_GEDI_vs_other_LiDAR = "TableS18_GEDI_vs_other_LiDAR", TableS27_GEDI_roles_by_study = "TableS19_GEDI_roles_by_study",
  TableS28_GEDI_products_verified = "TableS20_GEDI_products_verified", TableS33_contrast_precision = "TableS22_contrast_precision",
  TableS34_contrasts_robustness = "TableS23_contrasts_robustness", TableS35_LiDAR_pairing = "TableS24_LiDAR_pairing",
  TableS36_LiDAR_paired_summary = "TableS25_LiDAR_paired_summary", TableS37_LiDAR_strict_pairs = "TableS26_LiDAR_matched_pairs",
  TableS38_LiDAR_pairs_heterogeneity = "TableS27_LiDAR_pairs_heterogeneity", TableS39_LiDAR_pairs_sign_test = "TableS28_LiDAR_pairs_sign_test",
  TableS42_central_africa = "TableS29_central_africa", TableS43_sensitivity = "TableS30_sensitivity",
  TableS44_sensitivity_additional = "TableS31_sensitivity_additional", TableS45_clustering_levels = "TableS32_clustering_levels",
  TableS48_multilevel_Egger = "TableS33_multilevel_Egger", TableS50_RMSE_scale_checks = "TableS34_RMSE_scale_checks",
  TableS46_leave_one_study_out = "TableS30_leave_one_study_out", TableS47_leave_one_cluster_out = "TableS30_leave_one_cluster_out",
  TableS21_subgroups_validation_resolution_region = "FigS06_data_validation_resolution_region",
  TableS35_model_diagnostics = "TableS35_model_diagnostics", TableS35_cluster_influence = "TableS35_cluster_influence",
  TableS29_GEDI_test_verified = "FigS01_data_GEDI_test_verified", TableS41_terrain_sample_size_reference = "FigS05_data_terrain_sample_size_reference")

export_table <- function(df, name, translate = TRUE) {
  if (translate) df <- en_table(df) else { df <- as.data.frame(df)
    for (c in names(df)) if (is.numeric(df[[c]])) df[[c]] <- round(df[[c]], 4) else if (is.character(df[[c]]) && name != "TableS02_data_dictionary")
      for (rx in names(VAL_RX)) df[[c]] <- gsub(rx, VAL_RX[[rx]], df[[c]]) }
  dir_out <- TAB_DIR
  if (grepl("^TableS", name)) { if (name %in% names(OUT_TAB)) name <- OUT_TAB[[name]] else dir_out <- CACHE_DIR }
  df <- as.data.frame(df)
  for (c in names(df)) if (is.double(df[[c]])) df[[c]] <- ifelse(is.na(df[[c]]), NA_character_, as.character(df[[c]]))
  readr::write_excel_csv(df, file.path(dir_out, paste0(name, ".csv")), na = "")
  cat("  ✓ table", name, sprintf("(%d rows)\n", nrow(df)))
  invisible(df)
}

val_design <- function(txt, obs = NULL) {
  t <- tolower(ifelse(is.na(txt), "", txt))
  dplyr::case_when(
    t == "" ~ "Not reported",
    grepl("^r² d'ajustement|^r² d'entra|^ajustement", t) ~ "Training fit",
    grepl("^validation ind.pendante", t) ~ "Independent validation",
    grepl("^cv|validation crois|loo", t) ~ "Cross-validation",
    grepl("hold-out|ind.pendante|externe|validation de la carte", t) ~ "Independent validation",
    TRUE ~ "Not reported")
}
is_nested <- function(txt) grepl("sous-ensemble de l'agr", ifelse(is.na(txt), "", txt), fixed = TRUE)
is_gedi_ref <- function(ref) grepl("GEDI", ifelse(is.na(ref), "", ref))
primary_r2 <- function(d) val_design(d$Validation_design_text) != "Training fit" & !is_nested(d$Validation_design_text) & !is_gedi_ref(d$Reference_type)
primary_rmse <- function(d) val_design_rmse(d$Validation_design_text) != "Training fit" & !is_nested(d$Validation_design_text) & !is_gedi_ref(d$Reference_type)
val_design_rmse <- function(txt) {
  t <- tolower(ifelse(is.na(txt), "", txt))
  ifelse(grepl("rmse en validation crois|rmse.*leave-one-out|rmse.*k-fold", t), "Cross-validation", val_design(txt))
}
