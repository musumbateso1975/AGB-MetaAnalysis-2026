if (!exists("dat_fz")) { source("R/00_setup.R"); source("R/01_data_import.R") }
source("R/00b_figure_style.R"); source("R/00c_export_helpers.R")
set.seed(2024)

cat("=== 03_fisher_z_meta.R — Two-Level Meta-Analysis ===\n\n")

cat("--- 1. Two-level random-effects model (REML) ---\n")

model_reml <- metafor::rma.mv(
  yi      = Fisher_Z,
  V       = vi_imp,
  random  = ~ 1 | Cluster / Study_ID / Obs_ID,
  method  = "REML",
  data    = dat_fz,
  sparse  = TRUE,
  control = list(rel.tol = 1e-8, iter.max = 1000)
)
ct_cr2 <- clubSandwich::coef_test(model_reml, vcov = "CR2", cluster = dat_fz$Cluster)
print(ct_cr2)

cat("\nModel summary:\n")
print(summary(model_reml))

pooled_Z  <- as.numeric(model_reml$b)
pooled_r  <- tanh(pooled_Z)
pooled_R2 <- tanh(pooled_Z)^2
ci_lo     <- pooled_Z - qt(0.975, ct_cr2$df_Satt) * ct_cr2$SE
ci_hi     <- pooled_Z + qt(0.975, ct_cr2$df_Satt) * ct_cr2$SE

cat("\n=== PRIMARY RESULT ===\n")
cat("Pooled Fisher's Z =", round(pooled_Z, 3), "\n")
cat("95% CI: [", round(ci_lo, 3), ",", round(ci_hi, 3), "]\n")
cat("Equivalent r =", round(pooled_r, 3), "\n")
cat("Equivalent R² ≈", round(pooled_R2, 3), "\n")
cat("Cohen magnitude:", ifelse(pooled_r >= 0.5, "LARGE",
                          ifelse(pooled_r >= 0.3, "MEDIUM", "SMALL")), "\n")

cat("\n--- 2. Heterogeneity (I²) ---\n")

W      <- diag(1 / dat_fz$vi_imp)
X      <- model.matrix(model_reml)
P      <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2_total <- sum(model_reml$sigma2) /
            (sum(model_reml$sigma2) + (model_reml$k - model_reml$p) / sum(diag(P)))

cat("Total I² =", round(I2_total * 100, 1), "%\n")
cat("σ² (between-study)   =", round(model_reml$sigma2[1], 4), "\n")
cat("σ² (within-study)    =", round(model_reml$sigma2[2], 4), "\n")
cat("p (heterogeneity QE) =", round(model_reml$QEp, 4), "\n")

pct_between <- sum(model_reml$sigma2[1:2]) / sum(model_reml$sigma2) * 100
pct_within  <- model_reml$sigma2[3] / sum(model_reml$sigma2) * 100
cat("Between-study variance: ", round(pct_between, 1), "%\n")
cat("Within-study variance:  ", round(pct_within, 1), "%\n")

cat("\n--- 3. Robust Variance Estimation (RVE) ρ = 0.8 ---\n")

rve_08 <- robumeta::robu(
  formula  = Fisher_Z ~ 1,
  data     = dat_fz,
  studynum = Study_ID,
  var.eff.size = vi_imp,
  rho      = 0.8,
  small    = TRUE
)
print(summary(rve_08))

cat("\n--- 4. RVE Sensitivity across ρ ---\n")

rve_05 <- robumeta::robu(
  Fisher_Z ~ 1, data = dat_fz, studynum = Study_ID,
  var.eff.size = vi_imp, rho = 0.5, small = TRUE
)

rve_10 <- robumeta::robu(
  Fisher_Z ~ 1, data = dat_fz, studynum = Study_ID,
  var.eff.size = vi_imp, rho = 1.0, small = TRUE
)

rve_sensitivity <- data.frame(
  rho     = c(0.5, 0.8, 1.0),
  Z_est   = round(c(rve_05$reg_table$b.r,
                    rve_08$reg_table$b.r,
                    rve_10$reg_table$b.r), 4),
  SE      = round(c(rve_05$reg_table$SE,
                    rve_08$reg_table$SE,
                    rve_10$reg_table$SE), 4),
  CI_lo   = round(c(rve_05$reg_table$CI.L,
                    rve_08$reg_table$CI.L,
                    rve_10$reg_table$CI.L), 4),
  CI_hi   = round(c(rve_05$reg_table$CI.U,
                    rve_08$reg_table$CI.U,
                    rve_10$reg_table$CI.U), 4)
)

cat("\nRVE Sensitivity (ΔZ should be < 0.03 across ρ):\n")
print(rve_sensitivity)
cat("ΔZ (max – min):", round(diff(range(rve_sensitivity$Z_est)), 4), "\n")

cat("\n--- 5. Leave-One-Study-Out Sensitivity ---\n")

loso_results <- data.frame()
vw <- mean(dat_fz$vi_imp)

for (sid in unique(dat_fz$Study_ID)) {
  dat_tmp <- dat_fz[dat_fz$Study_ID != sid, ]
  mod_tmp <- tryCatch(
    metafor::rma.mv(
      yi     = Fisher_Z,
      V      = vi_imp,
      random = ~ 1 | Cluster / Study_ID / Obs_ID,
      method = "REML",
      data   = dat_tmp, sparse = TRUE,
      control = list(rel.tol = 1e-8)
    ),
    error = function(e) NULL
  )
  mod_raw <- tryCatch(
    metafor::rma.mv(yi = R2, V = rep(V_R2, nrow(dat_tmp)), random = ~ 1 | Cluster / Study_ID / Obs_ID,
                    method = "REML", data = dat_tmp, sparse = TRUE),
    error = function(e) NULL
  )
  if (!is.null(mod_tmp)) {
    loso_results <- rbind(loso_results, data.frame(
      study_removed = sid,
      R2_untransformed = if (!is.null(mod_raw)) as.numeric(mod_raw$b) else NA,
      Z_est         = as.numeric(mod_tmp$b),
      CI_lo         = as.numeric(mod_tmp$ci.lb),
      CI_hi         = as.numeric(mod_tmp$ci.ub)
    ))
  }
}

loso_range <- range(loso_results$Z_est)
cat("LOSO Z range:", round(loso_range[1], 3), "to", round(loso_range[2], 3),
    "| ΔZ =", round(diff(loso_range), 4), "\n")
cat("Primary estimate robustness:",
    ifelse(diff(loso_range) < 0.05, "ROBUST (ΔZ < 0.05)", "SENSITIVE"), "\n")

export_table(loso_results, "TableS46_leave_one_study_out")
cat("LOSO results saved.\n")

meta_results <- list(
  model_reml      = model_reml,
  rve_08          = rve_08,
  rve_sensitivity = rve_sensitivity,
  loso_results    = loso_results,
  pooled_Z        = pooled_Z,
  pooled_r        = pooled_r,
  pooled_R2       = pooled_R2,
  I2_total        = I2_total
)
saveRDS(meta_results, file.path(CACHE_DIR, "meta_results.rds"))
cat("\n✓ Model objects cached\n")
cat("=== Meta-analysis complete ===\n")
