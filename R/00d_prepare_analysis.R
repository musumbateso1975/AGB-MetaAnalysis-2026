suppressPackageStartupMessages({ library(metafor); library(clubSandwich); library(dplyr); library(tidyr); library(readr) })
AUD <- "data"; OUT <- CACHE_DIR
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
source("R/00c_export_helpers.R")
TAB_NAME <- c(R01_modele_global.csv = "TableS10_pooled_model",
  R02_meta_regression_coefficients.csv = "TableS23_meta_regression", R02b_tests_Wald_robustes.csv = "TableS24_meta_regression_Wald",
  R02c_variance_expliquee.csv = "TableS25_variance_explained",
  R03_sous_groupes_IC_IP.csv = "Table2_subgroups_CI_PI", R03b_GEDI_exploratoire.csv = "TableS26_GEDI_vs_other_LiDAR",
  R04_contrastes_intra_etude.csv = "Table4_within_study_contrasts", R04b_contrastes_par_etude_Z.csv = "TableS31_contrasts_by_study_R2",
  R04c_contrastes_par_etude_RMSE.csv = "TableS32_contrasts_by_study_RMSE", R05_sensibilites.csv = "TableS43_sensitivity",
  R05b_LOSO_cluster.csv = "TableS47_leave_one_cluster_out", R06_biais_publication.csv = "TableS48_multilevel_Egger",
  R07_RMSE_descriptif_Mgha.csv = "TableS49_RMSE_by_sensor", R08_nouvelle_Table6.csv = "TableS40_climate_scale_LiDAR")
wcsv <- function(x, f) export_table(x, TAB_NAME[[f]])
r2 <- function(z) tanh(z)^2
set.seed(20261110)

fz0 <- read_csv(file.path(AUD, "DATA_FisherZ.csv"), col_types = cols(.default = "c"), na = character(0))
rm0 <- read_csv(file.path(AUD, "DATA_RMS.csv"), col_types = cols(.default = "c"), na = character(0))
fz0$Reference <- sub("^(.*? [0-9]{4}[a-z]?)\\. .*$", "\\1", fz0$Reference)
rm0$Reference <- sub("^(.*? [0-9]{4}[a-z]?)\\. .*$", "\\1", rm0$Reference)
fz0_all <- fz0; rm0_all <- rm0
fz0 <- fz0_all[primary_r2(fz0_all), ]
rm0 <- rm0_all[primary_rmse(rm0_all), ]
fz0_tf <- fz0_all[!is_nested(fz0_all$Validation_design_text) & !is_gedi_ref(fz0_all$Reference_type), ]

clus_map <- c(S1 = "C_AfriSAR_Betul", S58 = "C_AfriSAR_Betul", S5 = "C_AfriSAR_Betul", S57 = "C_AfriSAR_Betul",
              S6 = "C_Haldwani", S12 = "C_Haldwani", S41 = "C_Sonoma", S54 = "C_Sonoma", S22 = "C_Conghua", S52 = "C_Conghua",
              S24 = "C_Maine", S32 = "C_Maine")
prep <- function(d) d %>% mutate(
  Cluster = ifelse(Study_ID %in% names(clus_map), clus_map[Study_ID], Study_ID),
  Sensor = factor(Sensor_combination, c("LiDAR+Optical", "LiDAR-only", "Optical+SAR", "Optical-only", "SAR+LiDAR", "Optical+SAR+LiDAR", "SAR-only")),
  LiDAR = ifelse(grepl("LiDAR", Sensor_combination), "Avec LiDAR", "Sans LiDAR"),
  Fusion = ifelse(grepl("\\+", Sensor_combination), "Multi-capteurs", "Mono-capteur"),
  Model = factor(Model_category, c("Statistical", "Machine-learning", "Deep-learning")),
  Scale = factor(case_when(Study_scale %in% c("Local") ~ "Local", Study_scale %in% c("Regional", "Watershed") ~ "Régional",
                           TRUE ~ "National/multi-pays"), c("Local", "Régional", "National/multi-pays")),
  Climate = factor(Climate_zone, c("Tropical", "Subtropical", "Temperate", "Arid")),
  GEDI_clean = case_when(GEDI_used == "Yes" & GEDI_coding_doubtful == "TRUE" ~ NA_character_,
                         GEDI_used == "Yes" ~ "GEDI", TRUE ~ "Sans GEDI"))
fz <- prep(fz0) %>% mutate(Z = as.numeric(Fisher_Z), R2 = as.numeric(R2), n = suppressWarnings(as.numeric(Training_samples)))
make_vi <- function(n, cap = 1000, impute = TRUE) {
  nc <- if (is.finite(cap)) pmin(n, cap) else n
  if (impute) nc[is.na(nc)] <- median(nc, na.rm = TRUE)
  ifelse(is.na(nc), NA, 1 / (nc - 3))
}
A <- fz %>% mutate(vi = make_vi(n), n_imput = is.na(n))
A_all <- prep(fz0_tf) %>% mutate(Z = as.numeric(Fisher_Z), R2 = as.numeric(R2), n = suppressWarnings(as.numeric(Training_samples))) %>% mutate(vi = make_vi(n), n_imput = is.na(n))
cat("Analysed set:", nrow(A), "obs.,", n_distinct(A$Study_ID), "studies,", n_distinct(A$Cluster), "clusters ; imputed n :", sum(A$n_imput), "\n")

fit_ml <- function(d, mods = NULL, robust = TRUE) {
  fit1 <- function(ctrl) if (is.null(mods)) rma.mv(Z, vi, random = ~ 1 | Cluster / Study_ID / Obs_ID, data = d, method = "REML", sparse = TRUE, control = ctrl)
                         else rma.mv(Z, vi, mods = mods, random = ~ 1 | Cluster / Study_ID / Obs_ID, data = d, method = "REML", sparse = TRUE, control = ctrl)
  f <- tryCatch(fit1(list()), error = function(e) {
    message("fit_ml: nlminb did not converge (", nrow(d), " obs., ", dplyr::n_distinct(d$Study_ID), " studies); retry with optim/BFGS")
    tryCatch(fit1(list(optimizer = "optim", optmethod = "BFGS")), error = function(e2) {
      message("fit_ml: BFGS failed; retry with optim/Nelder-Mead"); fit1(list(optimizer = "optim", optmethod = "Nelder-Mead")) }) })
  f
}
i2_ml <- function(m, d) {
  W <- diag(1 / d$vi); X <- model.matrix(m)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  vt <- (m$k - m$p) / sum(diag(P)); s <- m$sigma2
  c(I2_total = 100 * sum(s) / (sum(s) + vt), I2_cluster = 100 * s[1] / (sum(s) + vt),
    I2_etude = 100 * s[2] / (sum(s) + vt), I2_obs = 100 * s[3] / (sum(s) + vt), v_typique = vt)
}
pooled <- function(d, label, robust = TRUE) {
  m <- fit_ml(d); nc <- n_distinct(d$Cluster)
  if (robust && nc >= 4) {
    ct <- coef_test(m, vcov = "CR2", cluster = d$Cluster)
    se <- ct$SE; df <- ct$df_Satt; est <- ct$beta
  } else { se <- m$se; df <- Inf; est <- as.numeric(coef(m)) }
  tq <- qt(.975, df); s2 <- sum(m$sigma2); tpi <- qt(.975, max(nc - 1, 1))
  i2 <- i2_ml(m, d)
  raw <- tryCatch({ fr <- function(...) rma.mv(R2, rep(V_R2, nrow(d)), random = ~ 1 | Cluster / Study_ID / Obs_ID, data = d, method = "REML", sparse = TRUE, ...)
                    mr <- tryCatch(fr(), error = function(e) fr(control = list(optimizer = "optim", optmethod = "Nelder-Mead")))
                    if (robust && nc >= 4) { cr <- coef_test(mr, vcov = "CR2", cluster = d$Cluster); c(cr$beta, cr$SE, cr$df_Satt, sum(mr$sigma2)) } else c(coef(mr), mr$se, Inf, sum(mr$sigma2)) },
                  error = function(e) c(mean(d$R2), NA, NA, NA))
  data.frame(Scenario = label, k_obs = nrow(d), Etudes = n_distinct(d$Study_ID), Clusters = nc,
             Z = est, IC_inf = est - tq * se, IC_sup = est + tq * se, IP_inf = est - tpi * sqrt(se^2 + s2), IP_sup = est + tpi * sqrt(se^2 + s2),
             R2 = raw[1], R2_IC_inf = raw[1] - qt(.975, raw[3]) * raw[2], R2_IC_sup = raw[1] + qt(.975, raw[3]) * raw[2],
             R2_Zscale = r2(est), R2_Zscale_IC_inf = r2(est - tq * se), R2_Zscale_IC_sup = r2(est + tq * se),
             R2_IP_inf = r2(max(est - tpi * sqrt(se^2 + s2), 0)), R2_IP_sup = r2(est + tpi * sqrt(se^2 + s2)),
             R2_IPb_inf = max(raw[1] - tpi * sqrt(raw[2]^2 + raw[4]), 0), R2_IPb_sup = min(raw[1] + tpi * sqrt(raw[2]^2 + raw[4]), 1),
             ddl = df, sigma2_cluster = m$sigma2[1], sigma2_etude = m$sigma2[2], sigma2_obs = m$sigma2[3],
             tau2_total = s2, Q = m$QE, I2_total = i2["I2_total"], Inference = ifelse(robust && nc >= 4, "CR2 robuste", "modèle"))
}

strict_pairs <- function(both, colsA, colsB) {
  key <- c("Study_ID", "Site_Name", "Model_algorithm", "n", "Reference_type")
  a <- both %>% filter(L) %>% select(all_of(key), all_of(colsA), res_A = Spatial_resolution_m, sens_A = Sensor_combination)
  b <- both %>% filter(!L) %>% select(all_of(key), all_of(colsB), res_B = Spatial_resolution_m)
  inner_join(a, b, by = key, relationship = "many-to-many") %>%
    filter(!is.na(n), sens_A == "LiDAR-only" | res_A == res_B) %>% select(-res_A, -res_B, -sens_A)
}
