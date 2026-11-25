suppressPackageStartupMessages({ library(metafor); library(clubSandwich); library(dplyr); library(readr) })
source("R/00d_prepare_analysis.R", encoding = "UTF-8")
ex <- function(x, f) export_table(x, f)
rt <- function(z) tanh(z)
rmA <- prep(rm0) %>% filter(RMSE_unit == "Mg_ha") %>% mutate(lR = log(as.numeric(RMSE)))

val_ok <- val_design(A$Validation_design_text, A$Obs_ID) %in% c("Independent validation", "Cross-validation")
A <- A %>% mutate(Estimand = ifelse(Simple_regression_R2 == "TRUE", "A", "B"), B1 = Estimand == "B" & val_ok)
row_est <- function(d, lab) {
  p <- pooled(d, lab, robust = n_distinct(d$Cluster) >= 4)
  data.frame(Estimand = lab, n_obs = nrow(d), k_studies = n_distinct(d$Study_ID), R2_median = round(median(d$R2), 3),
             R2_pooled = round(p$R2, 3), R2_CI = sprintf("%.2f–%.2f", p$R2_IC_inf, p$R2_IC_sup),
             R2_PI_untransformed = sprintf("%.2f–%.2f", p$R2_IPb_inf, p$R2_IPb_sup),
             r_PI_FisherZ = sprintf("%.2f–%.2f", rt(p$IP_inf), rt(p$IP_sup)), R2_PI_FisherZ = sprintf("%.2f–%.2f", p$R2_IP_inf, p$R2_IP_sup),
             I2_total = round(p$I2_total, 1), Inference = p$Inference)
}
E <- bind_rows(row_est(A, "All comparable R² (primary synthesis)"),
               row_est(filter(A, Estimand == "A"), "A: simple linear regression (R² = r²)"),
               row_est(filter(A, Estimand == "B"), "B: multi-predictor / machine-learning / predictive R²"),
               row_est(filter(A, B1), "B1: B with documented validation design"),
               row_est(filter(A, !n_imput), "All, reported n only (no imputed variance)"))
print(E); ex(E, "TableS11_estimand_classes")

L <- A %>% mutate(L = LiDAR == "Avec LiDAR")
both <- L %>% group_by(Study_ID) %>% filter(n_distinct(L) == 2) %>% ungroup()
P <- strict_pairs(both, c(Cluster = "Cluster", Obs_A = "Obs_ID", Z_A = "Z", R2_A = "R2", v_A = "vi"),
                  c(Obs_B = "Obs_ID", Z_B = "Z", R2_B = "R2", v_B = "vi")) %>% mutate(dZ = Z_A - Z_B, dR2 = R2_A - R2_B)
lidar_ref <- function(x) grepl("LiDAR", x)
S <- both %>% group_by(Study_ID) %>% summarise(Cluster = first(Cluster), Reference = first(Reference),
  n_models_with = sum(L), n_models_without = sum(!L), Reference_data = paste(sort(unique(Reference_type)), collapse = " + "),
  dR2_unpaired = mean(R2[L]) - mean(R2[!L]), .groups = "drop") %>%
  left_join(P %>% group_by(Study_ID) %>% summarise(n_pairs = n(), dR2_paired = mean(dR2), dZ_paired = mean(dZ), .groups = "drop"), by = "Study_ID") %>%
  mutate(Pairing = case_when(!is.na(n_pairs) & !lidar_ref(Reference_data) ~ "Matched pair, field reference",
                             !is.na(n_pairs) ~ "Matched pair, LiDAR-derived reference",
                             grepl("LiDAR", Reference_data) & grepl("Terrain", Reference_data) ~ "Not paired: different reference data",
                             TRUE ~ "Not paired: different n, algorithm or unknown n"))
print(as.data.frame(S)); ex(S %>% mutate(across(where(is.numeric), ~ round(.x, 3))), "TableS35_LiDAR_pairing")
tsum <- function(x) { t <- t.test(x); c(mean = mean(x), lo = t$conf.int[1], hi = t$conf.int[2]) }
st <- S %>% filter(!is.na(n_pairs)); stf <- st %>% filter(!lidar_ref(Reference_data))
kh <- function(d, rho) { v <- with(d, 2 * vbar * (1 - rho)); r <- rma(d$dZ_paired, v, method = "REML", test = "knha", control = list(stepadj = .5, maxiter = 1000)); c(coef(r), r$ci.lb, r$ci.ub) }
vb <- P %>% group_by(Study_ID) %>% summarise(vbar = mean((v_A + v_B) / 2), .groups = "drop"); st <- st %>% left_join(vb, by = "Study_ID")
rmP <- P %>% inner_join(rmA %>% select(Obs_A = Obs_ID, lA = lR), by = "Obs_A") %>% inner_join(rmA %>% select(Obs_B = Obs_ID, lB = lR), by = "Obs_B") %>%
  mutate(lnRR = lA - lB) %>% group_by(Cluster, Study_ID) %>% summarise(lnRR = mean(lnRR), Ref = first(Reference_type), .groups = "drop") %>%
  group_by(Cluster) %>% summarise(lnRR = mean(lnRR), Ref = first(Ref), .groups = "drop")
rr <- function(x) { t <- t.test(x); 100 * (exp(c(mean(x), t$conf.int)) - 1) }
loo <- sapply(st$Study_ID, function(s) mean(st$dR2_paired[st$Study_ID != s]))
q1 <- tsum(S$dR2_unpaired); q2 <- tsum(st$dR2_paired); q3 <- tsum(stf$dR2_paired); r1 <- rr(rmP$lnRR); r2f <- rr(rmP$lnRR[!lidar_ref(rmP$Ref)])
k0 <- kh(st, 0); k8 <- kh(st, 0.8)
SUM <- data.frame(
  Analysis = c("All studies with both configurations (unpaired means)", "Matched pairs (all reference data)", "Matched pairs, field reference only",
               "Matched pairs, pooled ΔZ (Knapp–Hartung), rho = 0", "Matched pairs, pooled ΔZ (Knapp–Hartung), rho = 0.8",
               "Matched pairs, RMSE change (%), by cluster", "Matched pairs, RMSE change (%), field reference"),
  k = c(nrow(S), nrow(st), nrow(stf), nrow(st), nrow(st), nrow(rmP), sum(!lidar_ref(rmP$Ref))),
  Estimate = round(c(q1[1], q2[1], q3[1], k0[1], k8[1], r1[1], r2f[1]), 3),
  CI_lower = round(c(q1[2], q2[2], q3[2], k0[2], k8[2], r1[2], r2f[2]), 3),
  CI_upper = round(c(q1[3], q2[3], q3[3], k0[3], k8[3], r1[3], r2f[3]), 3),
  Studies_with_gain = c(sum(S$dR2_unpaired > 0), sum(st$dR2_paired > 0), sum(stf$dR2_paired > 0), NA, NA, sum(rmP$lnRR < 0), sum(rmP$lnRR[!lidar_ref(rmP$Ref)] < 0)),
  Leave_one_out = c("", sprintf("%.3f–%.3f", min(loo), max(loo)), "", "", "", "", ""),
  Metric = c("ΔR²", "ΔR²", "ΔR²", "ΔZ", "ΔZ", "ΔRMSE %", "ΔRMSE %"))
print(SUM); ex(SUM, "TableS36_LiDAR_paired_summary")

m1 <- rma(Z, vi, data = A, method = "REML"); p3 <- pooled(A, "3-level")
kw <- function(g) { d <- A[!is.na(A[[g]]), ]; kruskal.test(d$Z, factor(d[[g]]))$p.value }
qm <- function(g) { d <- A[!is.na(A[[g]]), ]; rma(Z, vi, mods = as.formula(paste("~ factor(", g, ")")), data = d, method = "REML")$QMp }
wt <- function(g) { d <- A[!is.na(A[[g]]), ]; m <- fit_ml(d, mods = as.formula(paste("~ factor(", g, ")")))
  Wald_test(m, constraints = constrain_zero(2:m$p), vcov = "CR2", cluster = d$Cluster, test = "HTZ")$p_val }
best <- A %>% group_by(Study_ID) %>% arrange(desc(Best_model == "Yes"), desc(Z)) %>% slice(1) %>% ungroup()
mb <- rma(Z, vi, data = best, method = "REML")
ml1 <- rma(Z, vi, mods = ~ I(LiDAR == "Avec LiDAR"), data = A, method = "REML"); bl <- coef(ml1)
btw <- r2(bl[1] + bl[2]) - r2(bl[1])
f3 <- function(x, d = 2) formatC(x, format = "f", digits = d)
fp <- function(p) ifelse(p < 0.001, "< 0.001", formatC(p, format = "f", digits = 3))
C5 <- data.frame(
  Question = c("Pooled R² (95% CI)", "Width of the 95% CI of pooled Z", "Sensor combination associated with R²?",
               "Model category associated with R²?", "GEDI use associated with R²?", "Advantage of configurations with LiDAR (ΔR²)",
               "Pooled R² of the best model per study vs all models"),
  Conventional = c(sprintf("%s (%s–%s); %d independent observations", f3(r2(coef(m1))), f3(r2(m1$ci.lb)), f3(r2(m1$ci.ub)), m1$k),
                   f3(m1$ci.ub - m1$ci.lb, 3),
                   sprintf("Kruskal–Wallis p %s; single-level QM p %s", fp(kw("Sensor_combination")), fp(qm("Sensor_combination"))),
                   sprintf("Kruskal–Wallis p %s; single-level QM p %s", fp(kw("Model_category")), fp(qm("Model_category"))),
                   sprintf("Kruskal–Wallis p %s; single-level QM p %s", fp(kw("GEDI_used")), fp(qm("GEDI_used"))),
                   sprintf("All observations, single-level moderator: %+.2f (p %s)", btw, fp(ml1$pval[2])),
                   sprintf("Best model: %s (%d studies)", f3(r2(coef(mb))), nrow(best))),
  Dependence_aware = c(sprintf("%s (%s–%s); %d clusters, CR2", f3(p3$R2), f3(p3$R2_IC_inf), f3(p3$R2_IC_sup), p3$Clusters),
                       f3(p3$IC_sup - p3$IC_inf, 3),
                       sprintf("Robust Wald (HTZ) p %s", fp(wt("Sensor_combination"))),
                       sprintf("Robust Wald (HTZ) p %s", fp(wt("Model_category"))),
                       sprintf("Robust Wald (HTZ) p %s", fp(wt("GEDI_used"))),
                       sprintf("Within studies, matched pairs: %+.2f (%.2f to %.2f; %d studies)", q2[1], q2[2], q2[3], nrow(st)),
                       sprintf("All models: %s; prediction interval %s–%s", f3(p3$R2), f3(p3$R2_IP_inf), f3(p3$R2_IP_sup))))
print(C5); ex(C5, "TableS18_conventional_vs_dependence_aware")

A <- A %>% mutate(CA = Region == "Central-Africa", Trop = Climate_zone == "Tropical")
tr <- function(v) { d <- A %>% filter(!is.na(.data[[v]])); m <- fit_ml(d, mods = as.formula(paste("~", v)))
  ct <- coef_test(m, vcov = "CR2", cluster = d$Cluster); c(b = ct$beta[2], lo = ct$beta[2] - qt(.975, ct$df_Satt[2]) * ct$SE[2],
  hi = ct$beta[2] + qt(.975, ct$df_Satt[2]) * ct$SE[2], df = ct$df_Satt[2], p = ct$p_Satt[2]) }
tr2 <- function(v) { d <- A %>% filter(!is.na(.data[[v]])) %>% mutate(Z = R2, vi = V_R2); m <- fit_ml(d, mods = as.formula(paste("~", v)))
  ct <- coef_test(m, vcov = "CR2", cluster = d$Cluster); c(b = ct$beta[2], lo = ct$beta[2] - qt(.975, ct$df_Satt[2]) * ct$SE[2],
  hi = ct$beta[2] + qt(.975, ct$df_Satt[2]) * ct$SE[2], df = ct$df_Satt[2], p = ct$p_Satt[2]) }
tca <- tr("CA"); ttr <- tr("Trop"); tca2 <- tr2("CA"); ttr2 <- tr2("Trop")
caS <- A %>% filter(CA) %>% group_by(Study_ID) %>% summarise(Reference = first(Reference), Country = first(Country), n_obs = n(),
  Sensors = paste(sort(unique(Sensor_combination)), collapse = "; "), Models = paste(sort(unique(Model_category)), collapse = "; "),
  Reference_data = paste(sort(unique(Reference_type)), collapse = "; "), Validation_verified = paste(sort(unique(setdiff(na.omit(Validation_design_text), ""))), collapse = "; "),
  R2_median = round(median(R2), 2), .groups = "drop")
caS <- bind_rows(caS, data.frame(Study_ID = "Test", Reference = sprintf("Central Africa vs other regions: ΔZ = %.2f (%.2f to %.2f), CR2 df = %.1f, p = %.3f%s", tca["b"], tca["lo"], tca["hi"], tca["df"], tca["p"],
                 ifelse(tca["df"] < 4, "; Satterthwaite df < 4: robust test unreliable, descriptive only", ""))),
                 data.frame(Study_ID = "Test", Reference = sprintf("Tropical vs other climates: ΔZ = %.2f (%.2f to %.2f), CR2 df = %.1f, p = %.3f", ttr["b"], ttr["lo"], ttr["hi"], ttr["df"], ttr["p"])),
                 data.frame(Study_ID = "Test", Reference = sprintf("Central Africa vs other regions, untransformed R²: ΔR² = %.2f (%.2f to %.2f), CR2 df = %.1f, p = %.3f", tca2["b"], tca2["lo"], tca2["hi"], tca2["df"], tca2["p"])),
                 data.frame(Study_ID = "Test", Reference = sprintf("Tropical vs other climates, untransformed R²: ΔR² = %.2f (%.2f to %.2f), CR2 df = %.1f, p = %.3f", ttr2["b"], ttr2["lo"], ttr2["hi"], ttr2["df"], ttr2["p"])))
print(as.data.frame(caS)); ex(caS, "TableS42_central_africa")

covs <- bind_rows(lapply(c(0, 0.3, 0.6, 0.9), function(rho) {
  V <- impute_covariance_matrix(A$vi, cluster = A$Study_ID, r = rho)
  m <- rma.mv(Z, V, random = ~ 1 | Cluster / Study_ID / Obs_ID, data = A, method = "REML", sparse = TRUE)
  ct <- coef_test(m, vcov = "CR2", cluster = A$Cluster); tq <- qt(.975, ct$df_Satt)
  data.frame(rho = rho, R2_pooled = round(r2(ct$beta), 3), R2_CI = sprintf("%.2f–%.2f", r2(ct$beta - tq * ct$SE), r2(ct$beta + tq * ct$SE)),
             tau2_total = round(sum(m$sigma2), 4)) }))
print(covs); ex(covs, "TableS14_sampling_covariance")
CTZ <- read_csv(file.path(CACHE_DIR, "TableS31_contrasts_by_study_R2.csv"), show_col_types = FALSE)
names(CTZ)[1:3] <- c("Study_ID", "dZ", "dR2"); cname <- names(CTZ)[ncol(CTZ)]
MDE <- CTZ %>% group_by(Contrast = .data[[cname]]) %>% summarise(k = n(), mean_dR2 = mean(dR2), SE_dR2 = sd(dR2) / sqrt(n()),
  CI_width_dR2 = 2 * qt(.975, n() - 1) * SE_dR2, MDE80_dR2 = (qt(.975, n() - 1) + qt(.8, n() - 1)) * SE_dR2, .groups = "drop") %>%
  mutate(across(where(is.numeric) & !k, ~ round(.x, 3)),
         Reading = ifelse(abs(mean_dR2) >= MDE80_dR2, "Difference larger than the detectable threshold",
                          "Non-significant; differences up to the MDE cannot be excluded (not evidence of equivalence)"))
print(MDE); ex(MDE, "TableS33_contrast_precision")
cat("=== 14_estimands_paired_conventional.R complete ===\n")
