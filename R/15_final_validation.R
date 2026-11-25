suppressPackageStartupMessages({ library(metafor); library(clubSandwich); library(dplyr); library(readr); library(ggplot2) })
source("R/00d_prepare_analysis.R", encoding = "UTF-8")
source("R/00b_figure_style.R")
ex <- function(x, f) export_table(x, f)
f3 <- function(x, d = 3) round(x, d)
A <- A %>% mutate(Estimand = ifelse(Simple_regression_R2 == "TRUE", "A", "B"))
set.seed(20261119)

ctfit <- function(d, y = "Z", v = "vi") {
  f1 <- function(ctrl) rma.mv(d[[y]], d[[v]], random = ~ 1 | Cluster / Study_ID / Obs_ID, data = d, method = "REML", sparse = TRUE, control = ctrl)
  m <- tryCatch(f1(list()), error = function(e) { message("ctfit: nlminb did not converge (", nrow(d), " obs.); retry with optim/BFGS")
    f1(list(optimizer = "optim", optmethod = "BFGS")) })
  ct <- coef_test(m, vcov = "CR2", cluster = d$Cluster); tq <- qt(.975, ct$df_Satt)
  list(m = m, b = ct$beta, se = ct$SE, df = ct$df_Satt, lo = ct$beta - tq * ct$SE, hi = ct$beta + tq * ct$SE)
}
rowZ <- function(lab, d, note) { f <- ctfit(d)
  data.frame(Approach = lab, n_obs = nrow(d), k_studies = n_distinct(d$Study_ID), k_clusters = n_distinct(d$Cluster),
             Mean_R2 = f3(r2(f$b)), CI = sprintf("%.2f–%.2f", r2(f$lo), r2(f$hi)), df = f3(f$df, 1), Assumption = note) }
cl <- A %>% group_by(Cluster, Study_ID) %>% summarise(R2 = mean(R2), .groups = "drop") %>% group_by(Cluster) %>% summarise(R2 = mean(R2), .groups = "drop")
bt <- replicate(5000, mean(sample(cl$R2, replace = TRUE)))
tt <- t.test(cl$R2)
raw <- ctfit(A %>% mutate(vc = V_R2), y = "R2", v = "vc")
SA <- bind_rows(
  rowZ("A. Fisher's Z, n capped at 1,000, imputed n (complementary)", A, "R² treated as squared correlation; CR2 robust to variance misspecification"),
  rowZ("B. Restricted: simple linear regression only (R² = r²)", filter(A, Estimand == "A"), sprintf("Fisher's Z exact; %d studies", n_distinct(A$Study_ID[A$Estimand == "A"]))),
  data.frame(Approach = "C1. Cluster means of untransformed R², t interval", n_obs = nrow(A), k_studies = n_distinct(A$Study_ID), k_clusters = nrow(cl),
             Mean_R2 = f3(mean(cl$R2)), CI = sprintf("%.2f–%.2f", tt$conf.int[1], tt$conf.int[2]), df = nrow(cl) - 1,
             Assumption = sprintf("No sampling variance; equal cluster weights; 2.5–97.5%% of cluster means %.2f–%.2f", quantile(cl$R2, .025), quantile(cl$R2, .975))),
  data.frame(Approach = "C2. Cluster means of untransformed R², cluster bootstrap", n_obs = nrow(A), k_studies = n_distinct(A$Study_ID), k_clusters = nrow(cl),
             Mean_R2 = f3(mean(cl$R2)), CI = sprintf("%.2f–%.2f", quantile(bt, .025), quantile(bt, .975)), df = NA, Assumption = "Percentile bootstrap, 5,000 resamples of clusters"),
  data.frame(Approach = "C3. Three-level model on untransformed R², equal working weights, CR2", n_obs = nrow(A), k_studies = n_distinct(A$Study_ID), k_clusters = n_distinct(A$Cluster),
             Mean_R2 = f3(raw$b), CI = sprintf("%.2f–%.2f", raw$lo, raw$hi), df = f3(raw$df, 1), Assumption = "No transformation; dependence modelled"),
  rowZ("D1. Fisher's Z, n capped at 500", A %>% mutate(vi = make_vi(n, cap = 500)), "Weights"),
  rowZ("D2. Fisher's Z, no cap on n", A %>% mutate(vi = make_vi(n, cap = Inf)), "Weights"),
  rowZ("D3. Fisher's Z, equal working weights", A %>% mutate(vi = mean(vi)), "Weights"),
  rowZ("D4. Fisher's Z, reported n only (imputed n excluded)", filter(A, !n_imput), "Weights"),
  rowZ("D5. Fisher's Z, pseudo-replicated n excluded", filter(A, n_pseudoreplicated != "TRUE"), "Weights"))
m0 <- fit_ml(A)
tq3 <- qt(.975, n_distinct(A$Cluster) - 1); s3r <- sqrt(sum(raw$m$sigma2) + raw$se^2)
pp <- pooled(A, "p"); SA$PI_R2 <- ""; SA$PI_R2[1] <- sprintf("%.2f–%.2f", pp$R2_IP_inf, pp$R2_IP_sup); SA$PI_R2[5] <- sprintf("%.2f–%.2f (bounded at 1)", raw$b - tq3 * s3r, min(raw$b + tq3 * s3r, 1))
SA$Note <- ""; SA$Note[1] <- sprintf("Mean v_i = %.4f vs total tau² = %.3f: weights nearly equal", mean(A$vi), sum(m0$sigma2))
print(SA[, c(1, 2, 3, 5, 6)]); ex(SA, "TableS12_effect_size_strategies")

L <- A %>% mutate(L = LiDAR == "Avec LiDAR")
both <- L %>% group_by(Study_ID) %>% filter(n_distinct(L) == 2) %>% ungroup()
P <- strict_pairs(both, c(Cluster = "Cluster", ZA = "Z", RA = "R2", vA = "vi", SA_ = "Sensor_combination", Val = "Validation_design_text"),
                  c(ZB = "Z", RB = "R2", vB = "vi", SB_ = "Sensor_combination"))
ps <- P %>% group_by(Study_ID) %>% summarise(Cluster = first(Cluster), RefType = first(Reference_type), n = paste(sort(unique(n)), collapse = "/"),
  Algorithms = paste(sort(unique(Model_algorithm)), collapse = "; "), Without = paste(sort(unique(SB_)), collapse = "; "),
  With = paste(sort(unique(SA_)), collapse = "; "), Validation = paste(unique(val_design(ifelse(is.na(Val), "", Val), rep("", length(Val)))), collapse = "; "),
  n_pairs = n(), R2_without = mean(RB), R2_with = mean(RA), dR2 = mean(RA - RB), dZ = mean(ZA - ZB),
  dR2_add = mean((RA - RB)[SA_ == "LiDAR+Optical" & SB_ == "Optical-only"]),
  dR2_sub = mean((RA - RB)[SA_ == "LiDAR-only" & SB_ == "Optical-only"]),
  dR2_min = min(RA - RB), dR2_max = max(RA - RB), vbar = mean((vA + vB) / 2), .groups = "drop") %>%
  left_join(A %>% distinct(Study_ID, Ref = Reference), by = "Study_ID")
ps <- ps %>% mutate(se0 = sqrt(2 * vbar), se8 = sqrt(2 * vbar * 0.2),
  dR2_lo = tanh(atanh(sqrt(R2_without)) + dZ - 1.96 * se0)^2 - R2_without, dR2_hi = pmin(tanh(atanh(sqrt(R2_without)) + dZ + 1.96 * se0)^2, 1) - R2_without)
kh <- function(rho) rma(ps$dZ, 2 * ps$vbar * (1 - rho), method = "REML", test = "knha", control = list(stepadj = .5, maxiter = 1000))
k0 <- kh(0); k8 <- kh(0.8)
loo <- sapply(seq_len(nrow(ps)), function(i) mean(ps$dR2[-i]))
PT <- ps %>% transmute(Study_ID, Reference = Ref, Cluster, R2_without = f3(R2_without), R2_with = f3(R2_with), dR2 = f3(dR2),
  dR2_range_across_pairs = sprintf("%.2f to %.2f", dR2_min, dR2_max), dR2_95CI_sampling_rho0 = sprintf("%.2f to %.2f", dR2_lo, dR2_hi),
  n_pairs, n, Without, With, Algorithms, Reference_data = RefType, dR2_addition = f3(dR2_add), dR2_substitution = f3(dR2_sub), Validation, LOO_mean_without_study = f3(loo))
print(as.data.frame(PT[, 1:9])); ex(PT, "TableS37_LiDAR_strict_pairs")
HET <- data.frame(Statistic = c("Arithmetic mean ΔR² (t interval)", "Pooled ΔZ, Knapp–Hartung, rho = 0", "Pooled ΔZ, Knapp–Hartung, rho = 0.8",
                                "tau² (rho = 0)", "I² (rho = 0, %)", "Q test p (rho = 0)", "I² (rho = 0.8, %)", "Leave-one-out range of mean ΔR²",
                                "Studies sharing a cluster", "Largest single-study weight (rho = 0, %)",
                                "Mean ΔR², LiDAR added to optical (LiDAR+Optical vs Optical-only)", "Mean ΔR², LiDAR substituted for optical (LiDAR-only vs Optical-only)"),
  Value = c(sprintf("%.3f (%.3f to %.3f)", mean(ps$dR2), t.test(ps$dR2)$conf.int[1], t.test(ps$dR2)$conf.int[2]),
            sprintf("%.3f (%.3f to %.3f)", coef(k0), k0$ci.lb, k0$ci.ub), sprintf("%.3f (%.3f to %.3f)", coef(k8), k8$ci.lb, k8$ci.ub),
            sprintf("%.4f", k0$tau2), sprintf("%.1f", k0$I2), sprintf("%.3f", k0$QEp), sprintf("%.1f", k8$I2),
            sprintf("%.3f to %.3f", min(loo), max(loo)), ifelse(any(duplicated(ps$Cluster)), paste(ps$Study_ID[ps$Cluster %in% ps$Cluster[duplicated(ps$Cluster)]], collapse = ", "), "none"),
            sprintf("%.1f", max(weights(k0))),
            sprintf("%.3f (%d studies)", mean(ps$dR2_add, na.rm = TRUE), sum(!is.na(ps$dR2_add))),
            sprintf("%.3f (%d studies)", mean(ps$dR2_sub, na.rm = TRUE), sum(!is.na(ps$dR2_sub)))))
print(HET); ex(HET, "TableS38_LiDAR_pairs_heterogeneity")
fp <- ps %>% mutate(lab = sprintf("%s  %s", Study_ID, Ref)) %>% arrange(dR2)
fp$lab <- factor(fp$lab, fp$lab)
g <- ggplot(fp, aes(dR2, lab)) + geom_vline(xintercept = 0, colour = MUTED) +
  geom_errorbarh(aes(xmin = dR2_lo, xmax = dR2_hi), height = .2, colour = INK) +
  geom_point(aes(x = dR2_min), shape = 124, size = 3, colour = REDL) + geom_point(aes(x = dR2_max), shape = 124, size = 3, colour = REDL) +
  geom_point(shape = 15, size = 3.5, colour = INK) +
  geom_text(aes(x = 0.53, label = sprintf("%d pairs", n_pairs)), hjust = 1, size = 3.2, colour = MUTED) +
  geom_vline(xintercept = mean(ps$dR2), linetype = "dashed", colour = REDL) +
  scale_x_continuous(limits = c(-0.02, 0.53), breaks = seq(0, 0.4, 0.1)) +
  labs(x = "ΔR² (configuration with LiDAR − without LiDAR)", y = NULL,
       title = "Matched within-study LiDAR comparisons",
       subtitle = "Black: approximate sampling 95% CI (rho = 0; exploratory) · red: range across pairs · dashed: mean") +
  theme_pub() + theme(axis.text.x = element_text(face = "plain"), panel.grid.major.x = element_line(colour = GRIDC, linetype = "dashed", linewidth = .3))
save_pub(g, "FigS04_LiDAR_matched_pairs", h = 3.4)

M1 <- rma(Z, vi, data = A, method = "FE")
M2 <- rma(Z, vi, data = A, method = "REML")
M3 <- fit_ml(A); c4 <- coef_test(M3, vcov = "CR2", cluster = A$Cluster)
pi3 <- predict(M3); tpi <- qt(.975, n_distinct(A$Cluster) - 1); s3 <- sqrt(sum(M3$sigma2))
row <- function(mod, est, se, lo, hi, df, het, pil = NA, piu = NA, inf)
  data.frame(Model = mod, n_obs = nrow(A), k_studies = n_distinct(A$Study_ID), k_clusters = n_distinct(A$Cluster),
             Z = f3(est), SE = f3(se, 4), Z_CI = sprintf("%.3f–%.3f", lo, hi), CI_width_Z = f3(hi - lo), R2 = f3(r2(est)),
             R2_CI = sprintf("%.2f–%.2f", r2(lo), r2(hi)), df = df, Heterogeneity = het,
             R2_PI = ifelse(is.na(pil), "—", sprintf("%.2f–%.2f", r2(max(pil, 0)), r2(piu))), Inference = inf)
MT <- bind_rows(
  row("M1 Fixed effect, observations independent", coef(M1), M1$se, M1$ci.lb, M1$ci.ub, "∞ (z)", sprintf("Q = %.0f, p < 0.001", M1$QE), inf = "Model-based, ignores heterogeneity and dependence"),
  row("M2 Random effects, observations independent", coef(M2), M2$se, M2$ci.lb, M2$ci.ub, "∞ (z)", sprintf("tau² = %.3f; I² = %.1f%%", M2$tau2, M2$I2),
      predict(M2)$pi.lb, predict(M2)$pi.ub, "Model-based, ignores dependence"),
  row("M3 Three-level (cluster/study/observation), model-based SE", coef(M3), M3$se, M3$ci.lb, M3$ci.ub, "∞ (z)", sprintf("sum sigma² = %.3f", sum(M3$sigma2)),
      coef(M3) - tpi * sqrt(s3^2 + M3$se^2), coef(M3) + tpi * sqrt(s3^2 + M3$se^2), "Model-based, dependence modelled"),
  row("M4 = M3 with CR2 robust SE (Satterthwaite df)", c4$beta, c4$SE, c4$beta - qt(.975, c4$df_Satt) * c4$SE, c4$beta + qt(.975, c4$df_Satt) * c4$SE,
      sprintf("%.1f", c4$df_Satt), sprintf("sum sigma² = %.3f", sum(M3$sigma2)),
      coef(M3) - tpi * sqrt(s3^2 + c4$SE^2), coef(M3) + tpi * sqrt(s3^2 + c4$SE^2), "Robust to misspecified dependence, small-sample corrected"))
modt <- function(g) {
  d <- A[!is.na(A[[g]]), ]; f <- as.formula(paste("~ factor(", g, ")"))
  a <- rma(Z, vi, mods = f, data = d, method = "FE"); b <- rma(Z, vi, mods = f, data = d, method = "REML")
  m <- fit_ml(d, mods = f); w <- Wald_test(m, constraints = constrain_zero(2:m$p), vcov = "CR2", cluster = d$Cluster, test = "HTZ")
  c(M1 = a$QMp, M2 = b$QMp, M3 = m$QMp, M4 = w$p_val, df4 = w$df_denom) }
mt <- rbind(Sensor = modt("Sensor_combination"), Model = modt("Model_category"), LiDAR = modt("LiDAR"), GEDI = modt("GEDI_clean"))
fpv <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.3f", p))
MT2 <- data.frame(Model = MT$Model,
  Sensor_p = fpv(mt["Sensor", 1:4]), Model_category_p = fpv(mt["Model", 1:4]), LiDAR_p = fpv(mt["LiDAR", 1:4]), GEDI_p = fpv(mt["GEDI", 1:4]))
MT <- cbind(MT, MT2[, -1]); print(MT[, c(1, 5:11, 15:18)]); ex(MT, "Table3_models_M1_M4")

sel_best <- A %>% group_by(Study_ID) %>% arrange(desc(Best_model == "Yes"), desc(Z)) %>% slice(1) %>% ungroup()
sel_med <- A %>% group_by(Study_ID) %>% arrange(abs(Z - median(Z))) %>% slice(1) %>% ungroup()
fb <- ctfit(sel_best); fm <- ctfit(sel_med)
Ld <- A %>% mutate(Lx = as.numeric(LiDAR == "Avec LiDAR")) %>% group_by(Study_ID) %>% mutate(L_bar = mean(Lx), L_w = Lx - L_bar) %>% ungroup()
md <- fit_ml(Ld, mods = ~ L_w + L_bar); cd <- coef_test(md, vcov = "CR2", cluster = Ld$Cluster)
b0 <- coef(md)[1]
mA <- fit_ml(A, mods = ~ Estimand); cA <- coef_test(mA, vcov = "CR2", cluster = A$Cluster)
w23 <- (M3$ci.ub - M3$ci.lb) / (M2$ci.ub - M2$ci.lb)
MECH <- data.frame(
  Mechanism = c("Dependence (M2 → M3, same observations, model-based)", "Robust inference (M3 → M4, same fitted model)",
                "Heterogeneity (M1 → M2, same observations)", "Model selection (one model per study: median vs best)",
                "Comparability: LiDAR within studies (same three-level fit)", "Comparability: LiDAR between studies (same three-level fit)",
                "Metric definition (class B vs class A, same fit)", "Metric scale (Fisher's Z vs untransformed R²)"),
  What_changes = c("Random structure only", "Standard error only", "Between-observation variance only", "Which model represents each study",
                   "Within-study deviation of LiDAR use", "Study-level proportion of LiDAR configurations", "Estimand class", "Scale of synthesis"),
  Result = c(sprintf("CI width of Z ×%.2f; pooled R² %.2f → %.2f", w23, r2(coef(M2)), r2(coef(M3))),
             sprintf("SE %.4f → %.4f (×%.2f); df ∞ → %.1f", M3$se, c4$SE, c4$SE / M3$se, c4$df_Satt),
             sprintf("SE %.4f → %.4f; pooled R² %.2f → %.2f", M1$se, M2$se, r2(coef(M1)), r2(coef(M2))),
             sprintf("R² %.2f (median model) vs %.2f (best model), back-transformed Fisher's Z, %d studies", r2(fm$b), r2(fb$b), nrow(sel_best)),
             sprintf("ΔZ = %.2f (CR2 p %s); ΔR² at the intercept %+.2f", cd$beta[2], fpv(cd$p_Satt[2]), r2(b0 + cd$beta[2]) - r2(b0)),
             sprintf("ΔZ = %.2f (CR2 p %s); ΔR² %+.2f", cd$beta[3], fpv(cd$p_Satt[3]), r2(b0 + cd$beta[3]) - r2(b0)),
             sprintf("ΔZ = %.2f (CR2 p %s)", cA$beta[2], fpv(cA$p_Satt[2])),
             sprintf("R² %.2f (Z) vs %.2f (untransformed)", r2(c4$beta), raw$b)),
  Interpretation = c(ifelse(w23 > 1.2, "Ignoring dependence overstates precision", "Little effect on precision"),
                     ifelse(c4$SE / M3$se > 1.1, "Model-based SE too small", "For the pooled mean, model-based and robust SE agree; CR2 changes the moderator tests (Table 3 of the article)"),
                     "Fixed-effect precision is spurious given I² near 100%",
                     "Selecting the best model inflates the pooled accuracy",
                     ifelse(cd$p_Satt[2] < 0.05, "Within-study LiDAR difference is clear", "Within-study difference not established"),
                     ifelse(cd$p_Satt[3] < 0.05, "Between-study difference is clear", "Between-study difference masked by context"),
                     ifelse(cA$p_Satt[2] < 0.05, "Metric classes differ", "Metric classes do not differ detectably"),
                     sprintf("Back-transformed mean Z exceeds the mean of untransformed R² by %.2f (non-linear back-transformation): both are reported", r2(c4$beta) - raw$b)))
print(MECH[, c(1, 3)]); ex(MECH, "TableS17_mechanisms")
fd <- data.frame(Model = factor(c("M1", "M2", "M3", "M4"), c("M4", "M3", "M2", "M1")),
  est = r2(c(coef(M1), coef(M2), coef(M3), c4$beta)),
  lo = r2(c(M1$ci.lb, M2$ci.lb, M3$ci.lb, c4$beta - qt(.975, c4$df_Satt) * c4$SE)),
  hi = r2(c(M1$ci.ub, M2$ci.ub, M3$ci.ub, c4$beta + qt(.975, c4$df_Satt) * c4$SE)),
  pS = mt["Sensor", 1:4])
fd$lab <- ifelse(fd$pS < 0.001, "p < 0.001", sprintf("p = %.3f", fd$pS))
fd$grp <- factor(ifelse(fd$Model == "M4", "M4 (three-level + CR2)", "M1–M3 (standard models)"), c("M1–M3 (standard models)", "M4 (three-level + CR2)"))
fd$desc <- c("Fixed effect", "Random effects", "Three-level", "Three-level + CR2")
fd$ylab <- factor(sprintf("%s  %s", fd$Model, fd$desc), rev(sprintf("%s  %s", fd$Model, fd$desc)))
C5 <- c("M1–M3 (standard models)" = "#2a78d6", "M4 (three-level + CR2)" = "#eb6834")
th5 <- ggplot2::theme_minimal(base_size = 13) + ggplot2::theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
  panel.grid.major.x = element_line(colour = "#e3e6ea", linewidth = .35), axis.text = element_text(colour = INK, size = 12),
  axis.title = element_text(colour = INK, size = 12.5), plot.title = element_text(face = "bold", colour = INK, size = 13.5),
  legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_blank(), plot.title.position = "plot")
g5a <- ggplot(fd, aes(est, ylab, colour = grp)) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .22, linewidth = 1.1) + geom_point(size = 4.2, shape = 15) +
  geom_text(aes(label = sprintf("%.2f", est)), vjust = -1.05, size = 4.1, show.legend = FALSE) +
  scale_colour_manual(values = C5) + scale_x_continuous(limits = c(0.60, 0.76), breaks = seq(.60, .75, .05)) +
  labs(x = "Pooled R² (back-transformed Fisher's Z) with 95% CI", y = NULL, title = sprintf("Same %d observations, four analyses", nrow(A))) + th5
g5b <- ggplot(fd, aes(1, ylab)) +
  geom_label(aes(label = lab, fill = grp, colour = grp), size = 4.3, fontface = "bold", label.size = 0, label.padding = unit(.35, "lines"), show.legend = FALSE) +
  scale_fill_manual(values = c("M1–M3 (standard models)" = "#e8f1fc", "M4 (three-level + CR2)" = "#fdeee6")) + scale_colour_manual(values = C5) +
  scale_x_continuous(limits = c(.5, 1.5)) + labs(x = " ", y = NULL, title = "Sensor\ntest (p)") + th5 +
  theme(axis.text = element_blank(), panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank(), plot.title = element_text(face = "bold", size = 12.5, hjust = .5))
g5 <- patchwork::wrap_plots(g5a, g5b, widths = c(3, 1.1)) + patchwork::plot_layout(guides = "collect") & theme(legend.position = "bottom") & guides(colour = guide_legend(nrow = 2))
save_pub(g5, "Fig5_models_M1_M4", h = 4.2, w = 8, h1 = 4.6, w1 = 5.9)
cat("=== 15_final_validation.R complete ===\n")
