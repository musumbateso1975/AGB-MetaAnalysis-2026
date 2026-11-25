suppressPackageStartupMessages({ library(metafor); library(clubSandwich); library(dplyr); library(readr); library(ggplot2); library(tidyr) })
source("R/00c_export_helpers.R")
w1 <- function(x, f) export_table(x, f, translate = FALSE)

source("R/00d_prepare_analysis.R", encoding = "UTF-8")
source("R/00b_figure_style.R")
A <- A %>% mutate(Estimand = ifelse(Simple_regression_R2 == "TRUE", "A", "B"), vc = mean(vi))
cat("Corrected analysed set:", nrow(A), "obs.,", n_distinct(A$Study_ID), "studies,", n_distinct(A$Cluster), "clusters\n")
fpv <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))
f2 <- function(x) sprintf("%.2f", x)
fit3 <- function(d, scale = "R2", mods = NULL) {
  y <- if (scale == "R2") d$R2 else d$Z; v <- if (scale == "R2") rep(V_R2, nrow(d)) else d$vi
  go <- function(ctrl) if (is.null(mods)) rma.mv(y, v, random = ~ 1 | Cluster / Study_ID / Obs_ID, data = d, method = "REML", sparse = TRUE, control = ctrl)
       else rma.mv(y, v, mods = mods, random = ~ 1 | Cluster / Study_ID / Obs_ID, data = d, method = "REML", sparse = TRUE, control = ctrl)
  tryCatch(go(list()), error = function(e) tryCatch(go(list(optimizer = "optim", optmethod = "BFGS")),
           error = function(e2) go(list(optimizer = "optim", optmethod = "Nelder-Mead", maxit = 5000)))) }
pool <- function(d, lab, scale) {
  if (n_distinct(d$Cluster) < 4) return(data.frame(Scenario = lab, Scale = scale, n_obs = nrow(d), k_studies = n_distinct(d$Study_ID), k_clusters = n_distinct(d$Cluster), R2 = NA, CI = "too few clusters", df = NA))
  m <- fit3(d, scale); ct <- coef_test(m, vcov = "CR2", cluster = d$Cluster); tq <- qt(.975, ct$df_Satt)
  bk <- if (scale == "R2") identity else r2
  data.frame(Scenario = lab, Scale = scale, n_obs = nrow(d), k_studies = n_distinct(d$Study_ID), k_clusters = n_distinct(d$Cluster),
             R2 = round(bk(ct$beta), 3), CI = sprintf("%.2f–%.2f", bk(ct$beta - tq * ct$SE), bk(ct$beta + tq * ct$SE)), df = round(ct$df_Satt, 1)) }

cl <- A %>% group_by(Cluster, Study_ID) %>% summarise(R2 = mean(R2), .groups = "drop") %>% group_by(Cluster) %>% summarise(R2 = mean(R2), .groups = "drop")
set.seed(20261119); bt <- replicate(5000, mean(sample(cl$R2, replace = TRUE)))
mR <- fit3(A, "R2"); cR <- coef_test(mR, vcov = "CR2", cluster = A$Cluster); tR <- qt(.975, cR$df_Satt)
mZ <- fit3(A, "Z");  cZ <- coef_test(mZ, vcov = "CR2", cluster = A$Cluster); tZ <- qt(.975, cZ$df_Satt)
nc <- n_distinct(A$Cluster); tpi <- qt(.975, nc - 1)
piZ <- c(cZ$beta - tpi * sqrt(sum(mZ$sigma2) + cZ$SE^2), cZ$beta + tpi * sqrt(sum(mZ$sigma2) + cZ$SE^2))
piR <- c(cR$beta - tpi * sqrt(sum(mR$sigma2) + cR$SE^2), cR$beta + tpi * sqrt(sum(mR$sigma2) + cR$SE^2))
SC <- data.frame(
  Estimate = c("Mean R², three-level CR2, untransformed (primary)", "Mean of cluster-mean R², t interval", "Mean of cluster-mean R², cluster bootstrap",
               "Back-transformed mean Fisher's Z, three-level CR2 (complementary)", "Median of cluster-mean R²", "Median of all observed R²"),
  Value = round(c(cR$beta, mean(cl$R2), mean(cl$R2), r2(cZ$beta), median(cl$R2), median(A$R2)), 3),
  CI_95 = c(sprintf("%.2f–%.2f", cR$beta - tR * cR$SE, cR$beta + tR * cR$SE), sprintf("%.2f–%.2f", t.test(cl$R2)$conf.int[1], t.test(cl$R2)$conf.int[2]),
            sprintf("%.2f–%.2f", quantile(bt, .025), quantile(bt, .975)), sprintf("%.2f–%.2f", r2(cZ$beta - tZ * cZ$SE), r2(cZ$beta + tZ * cZ$SE)), "", ""),
  Interval_for_new_study = c(sprintf("%.2f–%.2f (normal on R² scale; upper limit %s)", piR[1], min(piR[2], 1), ifelse(piR[2] > 1, "truncated at 1", "inside 0–1")),
                             sprintf("2.5–97.5%% of cluster means: %.2f–%.2f", quantile(cl$R2, .025), quantile(cl$R2, .975)), "",
                             sprintf("%.2f–%.2f (normal on Z scale, quantiles back-transformed)", r2(max(piZ[1], 0)), r2(piZ[2])), "", ""),
  Note = c(sprintf("SE = %.4f, Satterthwaite df = %.1f, sum sigma² = %.4f", cR$SE, cR$df_Satt, sum(mR$sigma2)), "Each cluster weighted equally", "5,000 resamples of clusters",
           sprintf("SE(Z) = %.4f, df = %.1f, sum sigma² = %.4f; estimates the R² at the mean of Z (close to a median), not the mean R²", cZ$SE, cZ$df_Satt, sum(mZ$sigma2)), "", ""))
print(SC[, 1:3]); w1(SC, "TableS13_scale_and_estimand")

mods_list <- c(Sensor = "Sensor_combination", Model = "Model_category", LiDAR = "LiDAR", GEDI = "GEDI_clean", Climate = "Climate_zone", Scale = "Scale")
mtest <- function(g, scale) {
  d <- A[!is.na(A[[g]]), ]; f <- as.formula(paste("~ factor(", g, ")"))
  y <- if (scale == "R2") d$R2 else d$Z; v <- if (scale == "R2") rep(V_R2, nrow(d)) else d$vi
  # M1 and M2 (conventional comparators) keep the mean Fisher's Z variance as a stand-in sampling variance on the R2 scale, as in mm()
  v12 <- if (scale == "R2") d$vc else d$vi
  a <- rma(y, v12, mods = f, data = d, method = "FE"); b <- rma(y, v12, mods = f, data = d, method = "REML")
  m <- fit3(d, scale, f); w <- Wald_test(m, constraints = constrain_zero(2:m$p), vcov = "CR2", cluster = d$Cluster, test = "HTZ")
  data.frame(Moderator = names(mods_list)[mods_list == g], Scale = scale, levels = m$p, p_M1 = a$QMp, p_M2 = b$QMp, p_M3 = m$QMp, p_M4 = w$p_val,
             F_M4 = w$Fstat, df_num = w$df_num, df_den_M4 = w$df_denom) }
MT <- bind_rows(lapply(c("R2", "Z"), function(s) bind_rows(lapply(mods_list, mtest, scale = s))))
MT <- MT %>% group_by(Scale) %>% mutate(p_M4_Holm = p.adjust(p_M4, "holm")) %>% ungroup() %>%
  mutate(Fragility = case_when(df_den_M4 < 4 ~ "HTZ denominator df < 4: unreliable", p_M4_Holm >= .05 & p_M4 < .05 ~ "Not significant after Holm", p_M4 < .05 ~ "Significant after Holm", TRUE ~ "Not significant"))
MTp <- MT %>% mutate(across(starts_with("p_"), fpv), F_M4 = round(F_M4, 2), df_den_M4 = round(df_den_M4, 1))
print(as.data.frame(MTp)); w1(MTp, "TableS16_moderators_M1_M4_two_scales")
mm <- function(scale) {
  y <- if (scale == "R2") A$R2 else A$Z; v <- if (scale == "R2") A$vc else A$vi; bk <- if (scale == "R2") identity else r2
  M1 <- rma(y, v, method = "FE"); M2 <- rma(y, v, method = "REML"); M3 <- fit3(A, scale); c4 <- coef_test(M3, vcov = "CR2", cluster = A$Cluster); q4 <- qt(.975, c4$df_Satt)
  data.frame(Scale = scale, Model = c("M1 fixed effect", "M2 random effects", "M3 three-level (model-based)", "M4 three-level CR2"),
             Estimate = round(bk(c(coef(M1), coef(M2), coef(M3), c4$beta)), 3), SE_model_scale = round(c(M1$se, M2$se, M3$se, c4$SE), 4),
             CI = c(sprintf("%.3f–%.3f", bk(M1$ci.lb), bk(M1$ci.ub)), sprintf("%.3f–%.3f", bk(M2$ci.lb), bk(M2$ci.ub)), sprintf("%.3f–%.3f", bk(M3$ci.lb), bk(M3$ci.ub)),
                    sprintf("%.3f–%.3f", bk(c4$beta - q4 * c4$SE), bk(c4$beta + q4 * c4$SE))),
             CI_width_model_scale = round(c(M1$ci.ub - M1$ci.lb, M2$ci.ub - M2$ci.lb, M3$ci.ub - M3$ci.lb, 2 * q4 * c4$SE), 4),
             df = c("z", "z", "z", sprintf("%.1f", c4$df_Satt)),
             Heterogeneity = c(sprintf("Q = %.0f", M1$QE), sprintf("tau² = %.4f, I² = %.1f%%", M2$tau2, M2$I2), sprintf("sigma² = %s", paste(round(M3$sigma2, 4), collapse = " / ")), "same fit as M3")) }
MM <- bind_rows(mm("R2"), mm("Z")); print(MM[, 1:6]); w1(MM, "TableS15_M1_M4_two_scales")

pz <- predict(mZ)
PI <- data.frame(Item = c("Formula", "Variance components (Z): cluster / study / observation", "SE of mean (CR2)", "Quantile used", "PI on Z scale",
                          "PI back-transformed (R²)", "metafor predict() (z quantile, model-based SE), R²", "PI untransformed model (R²)", "Interpretation"),
                 Value = c("mean ± t(k_clusters − 1) × sqrt(sigma²_cluster + sigma²_study + sigma²_obs + SE²)",
                           paste(round(mZ$sigma2, 4), collapse = " / "), sprintf("%.4f", cZ$SE), sprintf("t(%d) = %.3f", nc - 1, tpi),
                           sprintf("%.3f to %.3f", piZ[1], piZ[2]), sprintf("%.2f–%.2f (lower Z bound %s 0)", r2(max(piZ[1], 0)), r2(piZ[2]), ifelse(piZ[1] > 0, ">", "≤")),
                           sprintf("%.2f–%.2f", r2(max(pz$pi.lb, 0)), r2(pz$pi.ub)), sprintf("%.2f–%.2f", piR[1], piR[2]),
                           "Range of true accuracies expected for a new study drawn from the population of published studies; includes observation-level variance (a single reported model). Not a measure of geographic transferability."))
print(PI)

L <- A %>% mutate(L = LiDAR == "Avec LiDAR"); both <- L %>% group_by(Study_ID) %>% filter(n_distinct(L) == 2) %>% ungroup()
Pp <- strict_pairs(both, c(Cluster = "Cluster", ZA = "Z", RA = "R2", vA = "vi", SA_ = "Sensor_combination", Val = "Validation_design_text"),
                   c(ZB = "Z", RB = "R2", vB = "vi", SB_ = "Sensor_combination"))
ps <- Pp %>% group_by(Study_ID) %>% summarise(Ref = first(Reference_type), d = mean(RA - RB), dZ = mean(ZA - ZB), vbar = mean((vA + vB) / 2),
  dAdd = mean((RA - RB)[SA_ == "LiDAR+Optical" & SB_ == "Optical-only"]), Val = first(na.omit(Val)), .groups = "drop") %>%
  mutate(Val = val_design(ifelse(is.na(Val), "", Val), rep("", n())))
grid <- bind_rows(lapply(c(0, .2, .5, .8, .95), function(rho) { k <- rma(ps$dZ, 2 * ps$vbar * (1 - rho), method = "REML", test = "knha", control = list(stepadj = .5, maxiter = 2000))
  data.frame(rho = rho, dZ = round(coef(k), 3), CI = sprintf("%.3f–%.3f", k$ci.lb, k$ci.ub), tau2 = round(k$tau2, 4), I2 = round(k$I2, 1), max_weight = round(max(weights(k)), 1)) }))
infl <- ps %>% mutate(mean_without = sapply(seq_len(n()), function(i) mean(d[-i])), change = mean_without - mean(d)) %>% select(Study_ID, Ref, d, mean_without, change, Val) %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
st <- binom.test(sum(ps$d > 0), nrow(ps))
LP <- list(grid = grid, infl = infl)
SG <- data.frame(Item = c("Studies favouring LiDAR", "Exact sign test (two-sided)", "Studies with field reference and documented validation design", sprintf("Mean ΔR² (%d studies, t interval)", nrow(ps))),
                 Value = c(sprintf("%d of %d", sum(ps$d > 0), nrow(ps)), sprintf("p = %.4f", st$p.value),
                           sprintf("%d", sum(!grepl("LiDAR", ps$Ref) & ps$Val %in% c("Independent validation", "Cross-validation"))),
                           sprintf("%.3f (%.3f to %.3f)", mean(ps$d), t.test(ps$d)$conf.int[1], t.test(ps$d)$conf.int[2])))
print(grid); print(infl); print(SG); w1(SG, "TableS39_LiDAR_pairs_sign_test")

G <- A %>% filter(GEDI_used == "Yes") %>% mutate(Role = case_when(GEDI_role == "Non vérifié" ~ "Not verified", TRUE ~ GEDI_role))
GR <- G %>% group_by(Role) %>% summarise(n_obs = n(), k_studies = n_distinct(Study_ID), studies = paste(sort(unique(Study_ID)), collapse = ", "), mean_R2 = round(mean(R2), 3), .groups = "drop")
GP <- G %>% filter(Role != "Not verified") %>% group_by(Product = GEDI_product) %>% summarise(n_obs = n(), k_studies = n_distinct(Study_ID), mean_R2 = round(mean(R2), 3), .groups = "drop")
print(GR); print(GP); w1(GR, "TableS27_GEDI_roles_by_study"); w1(GP, "TableS28_GEDI_products_verified")
Gv <- A %>% filter(!(GEDI_used == "Yes" & GEDI_role == "Non vérifié"))
gt <- function(d, scale) { d <- d[!is.na(d$GEDI_clean), ]; m <- fit3(d, scale, ~ factor(GEDI_clean))
  w <- Wald_test(m, constraints = constrain_zero(2), vcov = "CR2", cluster = d$Cluster, test = "HTZ"); data.frame(Scale = scale, n_obs = nrow(d), k = n_distinct(d$Study_ID), beta = round(coef(m)[2], 3), p = fpv(w$p_val), df_den = round(w$df_denom, 1)) }
GT <- bind_rows(gt(Gv, "R2"), gt(Gv, "Z")) %>% mutate(Note = "GEDI used vs not, observations with unverified GEDI role excluded")
print(GT); w1(GT, "TableS29_GEDI_test_verified")

best <- A %>% group_by(Study_ID) %>% arrange(desc(Best_model == "Yes"), desc(R2)) %>% slice(1) %>% ungroup()
smean <- A %>% group_by(Cluster, Study_ID) %>% summarise(R2 = mean(R2), Z = mean(Z), vi = mean(vi), Obs_ID = first(Obs_ID), .groups = "drop")
cat("=== 16_two_scales_and_pairs.R complete ===\n")
