source("R/00d_prepare_analysis.R", encoding = "UTF-8")
set.seed(2024)

A <- A %>% mutate(
  Validation = droplevels(factor(val_design(Validation_design_text, Obs_ID),
                      c("Independent validation", "Cross-validation", "Training fit", "Not reported"))),
  Resolution = factor(dplyr::case_when(Resolution_class %in% c("< 1 m", "1 - 5 m") ~ "≤ 5 m", Resolution_class == "5 - 30 m" ~ "5–30 m",
                                       Resolution_class == "30 - 100 m" ~ "30–100 m", Resolution_class %in% c("100 - 500 m", "> 500 m") ~ "> 100 m",
                                       TRUE ~ NA_character_), c("≤ 5 m", "5–30 m", "30–100 m", "> 100 m")),
  Region_group = factor(dplyr::case_when(Region == "East-Asia" ~ "East Asia", Region == "South-Asia" ~ "South Asia",
                                         Region == "Southeast-Asia" ~ "Southeast Asia", grepl("Africa", Region) ~ "Africa",
                                         grepl("America", Region) ~ "Americas", TRUE ~ "Other"),
                        c("East Asia", "South Asia", "Southeast Asia", "Africa", "Americas", "Other")),
  log_n = log(n), GEDI_clean = factor(GEDI_clean, c("Sans GEDI", "GEDI"), c("No GEDI", "GEDI")),
  Country_cluster = ifelse(Country %in% c("", "Multiple"), paste0("Multi_", Study_ID), Country))
cc <- A$Country_cluster
repeat {
  old <- cc
  for (g in unique(A$Cluster)) { i <- A$Cluster == g; cc[cc %in% cc[i]] <- min(cc[i]) }
  if (identical(old, cc)) break
}
A$Country_cluster <- cc

cwb_F <- function(d, X, idx, cl, B = 9999) {
  set.seed(20261123)
  m <- tryCatch(rma.mv(Z, vi, mods = X, intercept = FALSE, random = ~ 1 | Cluster / Study_ID / Obs_ID, data = d, method = "REML", sparse = TRUE), error = function(e) tryCatch(rma.mv(Z, vi, mods = X, intercept = FALSE, random = ~ 1 | Cluster / Study_ID / Obs_ID, data = d, method = "REML", sparse = TRUE, control = list(optimizer = "optim", optmethod = "BFGS")), error = function(e2) rma.mv(Z, vi, mods = X, intercept = FALSE, random = ~ 1 | Cluster / Study_ID / Obs_ID, data = d, method = "REML", sparse = TRUE, control = list(optimizer = "optim", optmethod = "Nelder-Mead"))))
  W <- solve(as.matrix(m$M)); y <- d$Z
  X0 <- X[, -idx, drop = FALSE]
  gls <- function(Xm, yv) solve(t(Xm) %*% W %*% Xm, t(Xm) %*% W %*% yv)
  Bread <- solve(t(X) %*% W %*% X); H <- Bread %*% t(X) %*% W
  cls <- split(seq_along(y), cl); G <- length(cls); q <- length(idx)
  Fstat <- function(yv) {
    b <- H %*% yv; e <- yv - X %*% b
    meat <- Reduce(`+`, lapply(cls, function(i) { u <- t(X[i, , drop = FALSE]) %*% W[i, i, drop = FALSE] %*% e[i]; u %*% t(u) }))
    V <- G / (G - 1) * Bread %*% meat %*% Bread
    bb <- b[idx]; as.numeric(t(bb) %*% solve(V[idx, idx, drop = FALSE], bb)) / q
  }
  F0 <- Fstat(y); b0 <- gls(X0, y); f0 <- X0 %*% b0; e0 <- y - f0
  gi <- match(cl, names(cls))
  Fb <- vapply(seq_len(B), function(r) { w <- sample(c(-1, 1), G, replace = TRUE); Fstat(f0 + w[gi] * e0) }, numeric(1))
  c(F_CR1 = F0, p_CWB = (1 + sum(Fb >= F0)) / (B + 1))
}

mod_test <- function(var, label, d0 = A, cluster_var = "Cluster", B = 9999) {
  d <- d0[!is.na(d0[[var]]), ]
  if (is.factor(d[[var]])) d[[var]] <- droplevels(d[[var]])
  form <- as.formula(paste("~", var)); m <- fit_ml(d, mods = form)
  cl <- d[[cluster_var]]; idx <- 2:length(coef(m))
  w <- Wald_test(m, constraints = constrain_zero(idx), vcov = "CR2", cluster = cl, test = "HTZ")
  X <- model.matrix(form, d); cw <- cwb_F(d, X, idx, cl, B)
  d2 <- d; d2$Z <- d2$R2; d2$vi <- V_R2; m2 <- fit_ml(d2, mods = form)
  w2 <- Wald_test(m2, constraints = constrain_zero(idx), vcov = "CR2", cluster = cl, test = "HTZ"); cw2 <- cwb_F(d2, X, idx, cl, B)
  lev <- if (is.factor(d[[var]])) {
    tb <- d %>% group_by(.data[[var]]) %>% summarise(k = n_distinct(.data[[cluster_var]]), .groups = "drop")
    list(levels = nlevels(d[[var]]), min_k = min(tb$k),
         per_level = paste(sprintf("%s: %d", tb[[var]], tb$k), collapse = "; "))
  } else list(levels = NA, min_k = NA, per_level = "continuous")
  data.frame(Moderator = label, Clustering = cluster_var, n_obs = nrow(d), k_studies = n_distinct(d$Study_ID),
             k_clusters = n_distinct(cl), levels = lev$levels, min_clusters_per_level = lev$min_k,
             clusters_per_level = lev$per_level, df_num = w$df_num,
             F_HTZ_R2 = w2$Fstat, df_den_R2 = w2$df_denom, p_HTZ_R2 = ifelse(w2$df_denom > 0, w2$p_val, NA), p_CWB_R2 = cw2[["p_CWB"]],
             F_HTZ = w$Fstat, df_den = w$df_denom, p_HTZ = ifelse(w$df_denom > 0, w$p_val, NA), p_CWB = cw[["p_CWB"]],
             Reliability = ifelse(!is.finite(min(w$df_denom, w2$df_denom)) | min(w$df_denom, w2$df_denom) < 4, "df < 4: HTZ unreliable, rely on CWB", "HTZ adequate"))
}

MT <- bind_rows(
  mod_test("Sensor", "Sensor combination"), mod_test("Model", "Model category"),
  mod_test("GEDI_clean", "GEDI used"), mod_test("Climate", "Climate zone"),
  mod_test("Scale", "Study scale"), mod_test("Region_group", "Region"),
  mod_test("Resolution", "Spatial resolution"), mod_test("Validation", "Validation design"),
  mod_test("log_n", "log(training sample size)"),
  mod_test("Sensor", "Sensor combination", cluster_var = "Country_cluster"))
export_table(MT, "TableS19_moderator_robustness")
print(MT[, c("Moderator", "Clustering", "k_clusters", "min_clusters_per_level", "df_den", "p_HTZ", "p_CWB")])

lev_est <- function(var) {
  d <- A[!is.na(A[[var]]), ]; d[[var]] <- droplevels(d[[var]])
  bind_rows(lapply(levels(d[[var]]), function(l) {
    s <- d[d[[var]] == l, ]; if (n_distinct(s$Study_ID) < 2) return(NULL)
    p <- pooled(s, l, robust = n_distinct(s$Cluster) >= 4)
    data.frame(Moderator = var, Level = l, n_obs = p$k_obs, k_studies = p$Etudes, k_clusters = p$Clusters,
               R2 = p$R2, R2_CI_lower = p$R2_IC_inf, R2_CI_upper = p$R2_IC_sup, df_Z = p$ddl,
               R2_PI_untransformed_lower = p$R2_IPb_inf, R2_PI_untransformed_upper = p$R2_IPb_sup,
               R2_PI_lower = p$R2_IP_inf, R2_PI_upper = p$R2_IP_sup,
               Reliability = ifelse(p$Clusters < 5 | p$ddl < 4, "Indicative (< 5 clusters or df < 4)", "Adequate"))
  }))
}
SGa <- bind_rows(lev_est("Validation"), lev_est("Resolution"), lev_est("Region_group"))
SGa$Moderator <- dplyr::recode(SGa$Moderator, Validation = "Validation design", Resolution = "Spatial resolution", Region_group = "Region")
export_table(SGa, "TableS21_subgroups_validation_resolution_region")

source("R/00b_figure_style.R")
P0 <- pooled(A, "all")
fg <- SGa %>% mutate(Moderator = factor(Moderator, c("Validation design", "Spatial resolution", "Region")),
                     lab_y = sprintf("%s  (k = %d)", Level, k_studies),
                     txt = sprintf("%.2f (%.2f–%.2f)   n = %d", R2, R2_CI_lower, R2_CI_upper, n_obs),
                     open = grepl("Indicative", Reliability)) %>%
  group_by(Moderator) %>% mutate(ord = rev(seq_len(n()))) %>% ungroup()
pS6 <- ggplot2::ggplot(fg, ggplot2::aes(y = reorder(lab_y, ord))) +
  ggplot2::geom_vline(xintercept = P0$R2, linetype = "dashed", colour = "#444444") +
  ggplot2::geom_segment(ggplot2::aes(x = R2_PI_untransformed_lower, xend = R2_PI_untransformed_upper, yend = reorder(lab_y, ord)), colour = "grey60", linewidth = .7) +
  ggplot2::geom_segment(ggplot2::aes(x = pmax(R2_PI_lower, 0), xend = pmin(R2_PI_upper, 1), yend = reorder(lab_y, ord)), colour = "grey35", linewidth = .35, linetype = "dotted", position = ggplot2::position_nudge(y = -.25)) +
  ggplot2::geom_segment(ggplot2::aes(x = R2_CI_lower, xend = R2_CI_upper, yend = reorder(lab_y, ord)), colour = "#1f5fae", linewidth = 2) +
  ggplot2::geom_point(ggplot2::aes(x = R2, shape = open), size = 3.2, colour = "#1f5fae", fill = "white", stroke = 1.1) +
  ggplot2::geom_text(ggplot2::aes(x = 1.03, label = txt), hjust = 0, size = 3.1) +
  ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21), guide = "none") +
  ggplot2::facet_wrap(~ Moderator, ncol = 1, scales = "free_y") +
  ggplot2::scale_x_continuous(breaks = seq(0, 1, .2), expand = c(0, 0)) + ggplot2::scale_y_discrete(expand = ggplot2::expansion(add = .6)) +
  ggplot2::coord_cartesian(xlim = c(0, 1), clip = "off") +
  ggplot2::labs(x = "Pooled R\u00b2 (thick: 95% CI; thin: 95% PI, untransformed; dotted: 95% PI, Fisher's Z)", y = NULL,
                title = "Pooled R\u00b2 by validation design, spatial resolution and region",
                subtitle = sprintf("Three-level model, CR2 inference; dashed line: overall pooled R\u00b2 (%.2f);\nopen symbol: < 5 clusters or df < 4", P0$R2)) +
  theme_pub() + ggplot2::theme(strip.text = ggplot2::element_text(face = "bold", hjust = 0, size = 11), strip.background = ggplot2::element_blank(),
                               plot.margin = ggplot2::margin(5, 125, 5, 5), panel.spacing = grid::unit(8, "pt"))
save_pub(pS6, "FigS06_validation_resolution_region", h = 7.2, w = 8.2)

pc <- function(d, lab, clv, rnd = ~ 1 | Cluster / Study_ID / Obs_ID) {
  m <- tryCatch(rma.mv(Z, vi, random = rnd, data = d, method = "REML", sparse = TRUE), error = function(e) tryCatch(rma.mv(Z, vi, random = rnd, data = d, method = "REML", sparse = TRUE, control = list(optimizer = "optim", optmethod = "BFGS")), error = function(e2) rma.mv(Z, vi, random = rnd, data = d, method = "REML", sparse = TRUE, control = list(optimizer = "optim", optmethod = "Nelder-Mead")))); ct <- coef_test(m, vcov = "CR2", cluster = d[[clv]])
  mr <- tryCatch(rma.mv(R2, rep(V_R2, nrow(d)), random = rnd, data = d, method = "REML", sparse = TRUE), error = function(e) tryCatch(rma.mv(R2, rep(V_R2, nrow(d)), random = rnd, data = d, method = "REML", sparse = TRUE, control = list(optimizer = "optim", optmethod = "BFGS")), error = function(e2) rma.mv(R2, rep(V_R2, nrow(d)), random = rnd, data = d, method = "REML", sparse = TRUE, control = list(optimizer = "optim", optmethod = "Nelder-Mead"))))
  cr <- coef_test(mr, vcov = "CR2", cluster = d[[clv]])
  data.frame(Clustering = lab, n_clusters = n_distinct(d[[clv]]),
             R2 = cr$beta, R2_CI_lower = cr$beta - qt(.975, cr$df_Satt) * cr$SE, R2_CI_upper = cr$beta + qt(.975, cr$df_Satt) * cr$SE,
             df_R2 = cr$df_Satt, Z = ct$beta, SE_Z = ct$SE, df_Z = ct$df_Satt)
}
CC <- bind_rows(pc(A, "Clusters of studies sharing field data (primary)", "Cluster"),
                pc(A, "Studies (shared field data ignored)", "Study_ID", ~ 1 | Study_ID / Obs_ID), pc(A, "Countries (countries linked by shared field data merged)", "Country_cluster"))
export_table(CC, "TableS45_clustering_levels")
print(CC)

rmA <- prep(rm0) %>% filter(RMSE_unit == "Mg_ha") %>% mutate(lR = log(as.numeric(RMSE)))
contrast_rob <- function(var, Aset, Bset, nom) {
  z <- A %>% filter(.data[[var]] %in% c(Aset, Bset)) %>% mutate(g = ifelse(.data[[var]] %in% Aset, "A", "B")) %>%
    group_by(Cluster, Study_ID) %>% filter(n_distinct(g) == 2) %>%
    summarise(dZ = mean(Z[g == "A"]) - mean(Z[g == "B"]), dR2 = mean(R2[g == "A"]) - mean(R2[g == "B"]),
              vA = mean(vi[g == "A"]) / sum(g == "A"), vB = mean(vi[g == "B"]) / sum(g == "B"), .groups = "drop")
  zc <- z %>% group_by(Cluster) %>% summarise(dR2 = mean(dR2), .groups = "drop")
  rho_fit <- function(rho) {
    v <- z$vA + z$vB - 2 * rho * sqrt(z$vA * z$vB)
    f <- tryCatch(rma(z$dZ, v, method = "REML", test = "knha", control = list(stepadj = 0.5, maxiter = 1000)), error = function(e) rma(z$dZ, v, method = "DL", test = "knha"))
    c(dZ = as.numeric(coef(f)), lo = f$ci.lb, hi = f$ci.ub, p = f$pval)
  }
  r0 <- rho_fit(0); r5 <- rho_fit(.5); r9 <- rho_fit(.9)
  tc <- if (nrow(zc) >= 3) t.test(zc$dR2) else NULL
  st <- binom.test(sum(z$dR2 > 0), sum(z$dR2 != 0))
  data.frame(Contrast = nom, k_studies = nrow(z), k_clusters = nrow(zc),
             dZ_rho0 = r0[["dZ"]], p_rho0 = r0[["p"]], dZ_rho05 = r5[["dZ"]], CI_rho05 = sprintf("%.2f to %.2f", r5[["lo"]], r5[["hi"]]), p_rho05 = r5[["p"]],
             dZ_rho09 = r9[["dZ"]], p_rho09 = r9[["p"]],
             dR2_cluster_mean = mean(zc$dR2), dR2_cluster_CI_lower = if (!is.null(tc)) tc$conf.int[1] else NA,
             dR2_cluster_CI_upper = if (!is.null(tc)) tc$conf.int[2] else NA,
             k_A_higher = sum(z$dR2 > 0), p_sign_test = st$p.value)
}
CR <- bind_rows(contrast_rob("Fusion", "Multi-capteurs", "Mono-capteur", "Multi-sensor vs single-sensor"),
                contrast_rob("LiDAR", "Avec LiDAR", "Sans LiDAR", "With LiDAR vs without LiDAR"),
                contrast_rob("Model_category", "Machine-learning", "Statistical", "Machine learning vs statistical"),
                contrast_rob("Model_category", "Deep-learning", c("Machine-learning", "Statistical"), "Deep learning vs other models"))
export_table(CR, "TableS34_contrasts_robustness")
print(CR)

P0 <- pooled(A, "Primary")
sc <- function(d, lab) { p <- pooled(d, lab)
  data.frame(Scenario = lab, n_obs = p$k_obs, k_studies = p$Etudes, k_clusters = p$Clusters, R2 = p$R2,
             R2_CI_lower = p$R2_IC_inf, R2_CI_upper = p$R2_IC_sup, df_Z = p$ddl, R2_PI_untransformed_lower = p$R2_IPb_inf, R2_PI_untransformed_upper = p$R2_IPb_sup,
             R2_PI_lower = p$R2_IP_inf, R2_PI_upper = p$R2_IP_sup, Change_vs_primary = p$R2 - P0$R2) }
mkA <- function(d) prep(d) %>% mutate(Z = as.numeric(Fisher_Z), R2 = as.numeric(R2), n = suppressWarnings(as.numeric(Training_samples))) %>% mutate(vi = make_vi(n), n_imput = is.na(n))
sc_w <- function(d, lab) { A_keep <- A; A <<- d; on.exit(A <<- A_keep); sc(d, lab) }
nest_rows <- fz0_all[is_nested(fz0_all$Validation_design_text), ]
agg_ids <- unique(sub(".*agrégat (M[0-9]+).*", "\\1", nest_rows$Validation_design_text))
gedi_rows <- fz0_all[is_gedi_ref(fz0_all$Reference_type) & val_design(fz0_all$Validation_design_text) != "Training fit", ]
n_train <- c(S46 = 130, S8 = 62, S20 = 77, S37 = 78, S38 = 30)
A_a <- A; i_a <- A_a$Study_ID %in% names(n_train); A_a$n[i_a] <- n_train[A_a$Study_ID[i_a]]; A_a$vi <- make_vi(A_a$n)
A_c <- A; A_c$vi <- make_vi(A_c$n, cap = Inf)
SA <- bind_rows(
  sc(A, "Primary analysis"),
  sc(A_all, "Including training-fit R²"),
  sc_w(A_a, "Training n instead of test-set n"),
  sc(mkA(bind_rows(fz0[!fz0$Obs_ID %in% agg_ids, ], nest_rows)), "Nested subsets instead of aggregates"),
  sc(mkA(bind_rows(fz0, nest_rows)), "Aggregates and nested subsets"),
  sc_w(A_c, "No cap on n"),
  sc(mkA(bind_rows(fz0, gedi_rows)), "Including R² computed against GEDI"),
  sc(filter(A, Map_validation != "TRUE"), "Excluding map validations"),
  sc(filter(A, Reference_type != "Terrain (présumé)"), "Excluding presumed field references"),
  sc(filter(A, Validation %in% c("Independent validation", "Cross-validation")), "Cross-validation or independent validation only"),
  sc(filter(A, Validation == "Independent validation"), "Independent validation only"),
  sc(filter(A, R2_estimand_different != "TRUE"), "Excluding R² flagged as a different estimand"),
  sc(filter(A, Region_group == "Africa"), "Africa only"),
  sc(filter(A, Region_group == "Americas"), "Americas only"),
  sc(filter(A, Region_group == "South Asia"), "South Asia only"),
  sc(filter(A, Region_group == "East Asia"), "East Asia only"))
export_table(SA, "TableS44_sensitivity_additional")
print(SA[, c("Scenario", "n_obs", "k_clusters", "R2", "R2_CI_lower", "R2_CI_upper", "df_Z")])

fp <- function(p) ifelse(p < 0.001, "< 0.001", sprintf("%.2f", p))
sk <- function(x) mean((x - mean(x))^3) / sd(x)^3
rs <- rmA %>% mutate(R = as.numeric(RMSE))
st_med <- rs %>% group_by(Study_ID) %>% summarise(R = median(R), .groups = "drop")
pair <- inner_join(select(A, Obs_No, R2), select(rs, Obs_No, R), by = "Obs_No")
RM <- data.frame(
  Check = c("Unit of all analysed RMSE values", "Observations / studies", "Skewness, observation level (Mg/ha)",
            "Skewness, observation level (log scale)", "Skewness of study medians (Mg/ha)", "Skewness of study medians (log scale)",
            "Shapiro-Wilk p, study medians (Mg/ha)", "Shapiro-Wilk p, study medians (log scale)",
            "Spearman rho between R² and RMSE for the same model", "Mean AGB of the validation plots available"),
  Value = c(paste(unique(rmA$AGB_unit), collapse = ", "), sprintf("%d / %d", nrow(rs), n_distinct(rs$Study_ID)),
            sprintf("%.2f", sk(rs$R)), sprintf("%.2f", sk(log(rs$R))), sprintf("%.2f", sk(st_med$R)), sprintf("%.2f", sk(log(st_med$R))),
            fp(shapiro.test(st_med$R)$p.value), fp(shapiro.test(log(st_med$R))$p.value),
            sprintf("%.2f (%d models)", cor(pair$R2, pair$R, method = "spearman"), nrow(pair)), "No (relative RMSE cannot be computed)"))
export_table(RM, "TableS50_RMSE_scale_checks")
print(RM)

cv <- function(a, b) { t <- table(a, b); t <- t[rowSums(t) > 0, colSums(t) > 0, drop = FALSE]
  x <- suppressWarnings(chisq.test(t, correct = FALSE)$statistic); sqrt(as.numeric(x) / (sum(t) * (min(dim(t)) - 1))) }
prs <- list(c("Sensor", "Region_group"), c("Sensor", "Model"), c("Sensor", "Climate"), c("Sensor", "Resolution"), c("Sensor", "Validation"),
            c("Region_group", "Climate"), c("Region_group", "Model"), c("Region_group", "Validation"), c("Model", "Validation"))
lab <- c(Sensor = "Sensor combination", Region_group = "Region", Model = "Model category", Climate = "Climate zone",
         Resolution = "Spatial resolution", Validation = "Validation design")
CF <- bind_rows(lapply(prs, function(p) { d <- A[!is.na(A[[p[1]]]) & !is.na(A[[p[2]]]), ]
  data.frame(Pair = paste(lab[[p[1]]], "×", lab[[p[2]]]), n_obs = nrow(d), Cramers_V = round(cv(d[[p[1]]], d[[p[2]]]), 2),
             Reading = ifelse(cv(d[[p[1]]], d[[p[2]]]) >= .5, "Strong association: effects not separable",
                       ifelse(cv(d[[p[1]]], d[[p[2]]]) >= .3, "Moderate association: partial confounding", "Weak association"))) }))
export_table(CF, "TableS22_moderator_confounding")
print(CF)
