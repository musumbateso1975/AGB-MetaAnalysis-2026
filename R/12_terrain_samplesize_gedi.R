suppressPackageStartupMessages({ library(metafor); library(clubSandwich); library(dplyr); library(readr); library(ggplot2) })
source("R/00d_prepare_analysis.R", encoding = "UTF-8")
source("R/00b_figure_style.R")

rmA <- prep(rm0) %>% filter(RMSE_unit == "Mg_ha") %>%
  mutate(lR = log(as.numeric(RMSE)), n = suppressWarnings(as.numeric(Training_samples)))
terr <- function(x) factor(case_when(x %in% c("Flat / Lowland", "Coastal", "Plateau") ~ "Flat / lowland / plateau",
                                     x == "Hilly" ~ "Hilly",
                                     x %in% c("Mountainous", "Complex terrain") ~ "Mountainous / complex",
                                     TRUE ~ NA_character_),
                           c("Flat / lowland / plateau", "Hilly", "Mountainous / complex"))
ncls <- function(n) factor(cut(n, c(0, 50, 100, 200, Inf), labels = c("≤ 50", "51–100", "101–200", "> 200")),
                           c("≤ 50", "51–100", "101–200", "> 200"))
reft <- function(x) factor(case_when(grepl("^Terrain", x) ~ "Field plots",
                                     grepl("LiDAR", x) ~ "LiDAR-derived AGB", TRUE ~ NA_character_),
                           c("Field plots", "LiDAR-derived AGB"))
A   <- A   %>% mutate(Terrain = terr(Terrain_type), Ncls = ncls(n), RefT = reft(Reference_type))
rmA <- rmA %>% mutate(Terrain = terr(Terrain_type), Ncls = ncls(n), RefT = reft(Reference_type))

wz <- function(d, g, raw = FALSE) {
  d <- d[!is.na(d[[g]]), ]; d[[g]] <- droplevels(d[[g]])
  if (raw) { d$Z <- d$R2; d$vi <- V_R2 }
  m <- rma.mv(Z, vi, mods = as.formula(paste("~", g)), random = ~ 1 | Cluster / Study_ID / Obs_ID, data = d, method = "REML", sparse = TRUE)
  w <- Wald_test(m, constraints = constrain_zero(2:m$p), vcov = "CR2", cluster = d$Cluster, test = "HTZ")
  c(F = w$Fstat, df_num = w$df_num, df_den = w$df_denom, p = w$p_val)
}
wr <- function(d, g) {
  d <- d[!is.na(d[[g]]) & is.finite(d$lR), ]; d[[g]] <- droplevels(d[[g]])
  out <- tryCatch({
    m <- nlme::lme(as.formula(paste("lR ~", g)), random = ~ 1 | Study_ID, data = d, method = "REML")
    w <- Wald_test(m, constraints = constrain_zero(2:length(nlme::fixef(m))), vcov = "CR2", test = "HTZ")
    c(F = w$Fstat, df_num = w$df_num, df_den = w$df_denom, p = w$p_val) }, error = function(e) c(F = NA, df_num = NA, df_den = NA, p = NA))
  if (!is.finite(out["F"]) || out["F"] > 1e4 || !is.finite(out["df_den"]) || out["df_den"] < 2) out[] <- NA
  out
}
dn <- A %>% filter(!n_imput, n > 3)
mn <- rma.mv(Z, vi, mods = ~ log(n), random = ~ 1 | Cluster / Study_ID / Obs_ID, data = dn, method = "REML", sparse = TRUE)
cn <- coef_test(mn, vcov = "CR2", cluster = dn$Cluster)
mn2 <- tryCatch(rma.mv(R2, rep(V_R2, nrow(dn)), mods = ~ log(n), random = ~ 1 | Cluster / Study_ID / Obs_ID, data = dn, method = "REML", sparse = TRUE), error = function(e) tryCatch(rma.mv(R2, rep(V_R2, nrow(dn)), mods = ~ log(n), random = ~ 1 | Cluster / Study_ID / Obs_ID, data = dn, method = "REML", sparse = TRUE, control = list(optimizer = "optim", optmethod = "BFGS")), error = function(e2) rma.mv(R2, rep(V_R2, nrow(dn)), mods = ~ log(n), random = ~ 1 | Cluster / Study_ID / Obs_ID, data = dn, method = "REML", sparse = TRUE, control = list(optimizer = "optim", optmethod = "Nelder-Mead"))))
cn2 <- coef_test(mn2, vcov = "CR2", cluster = dn$Cluster)
dnr <- rmA %>% filter(is.finite(n), n > 3, is.finite(lR))
mnr <- nlme::lme(lR ~ log(n), random = ~ 1 | Study_ID, data = dnr, method = "REML")
cnr <- coef_test(mnr, vcov = "CR2")

rows <- list()
for (g in c("Terrain", "Ncls", "RefT")) {
  lab <- c(Terrain = "Terrain", Ncls = "Training sample size (n)", RefT = "Reference data")[[g]]
  tz <- wz(A, g); t2 <- wz(A, g, raw = TRUE); tr <- wr(rmA, g)
  for (l in levels(A[[g]])) {
    d <- A[!is.na(A[[g]]) & A[[g]] == l, ]; r <- rmA[!is.na(rmA[[g]]) & rmA[[g]] == l, ]
    k <- n_distinct(d$Study_ID)
    pz <- if (k >= 2) pooled(d, l, robust = n_distinct(d$Cluster) >= 4) else NULL
    rs <- r %>% group_by(Study_ID) %>% summarise(m = median(as.numeric(RMSE)), .groups = "drop")
    rows[[paste(g, l)]] <- data.frame(Moderator = lab, Level = l, n_obs = nrow(d), k_studies = k,
      R2 = if (!is.null(pz)) pz$R2 else NA, R2_CI_lower = if (!is.null(pz)) pz$R2_IC_inf else NA, R2_CI_upper = if (!is.null(pz)) pz$R2_IC_sup else NA,
      R2_PI_untransformed_lower = if (!is.null(pz)) pz$R2_IPb_inf else NA, R2_PI_untransformed_upper = if (!is.null(pz)) pz$R2_IPb_sup else NA,
      R2_PI_lower = if (!is.null(pz)) pz$R2_IP_inf else NA, R2_PI_upper = if (!is.null(pz)) pz$R2_IP_sup else NA,
      n_obs_RMSE = nrow(r), k_studies_RMSE = nrow(rs), RMSE_median = if (nrow(rs)) median(rs$m) else NA,
      RMSE_Q1 = if (nrow(rs)) quantile(rs$m, .25) else NA, RMSE_Q3 = if (nrow(rs)) quantile(rs$m, .75) else NA,
      Reliability = ifelse(k >= 5, "≥ 5 studies", "< 5 studies: indicative"),
      Test_R2 = sprintf("F = %.2f (%.0f, %.1f), p = %.3f", t2["F"], t2["df_num"], t2["df_den"], t2["p"]),
      Test_Z = sprintf("F = %.2f (%.0f, %.1f), p = %.3f", tz["F"], tz["df_num"], tz["df_den"], tz["p"]),
      Test_logRMSE = if (is.na(tr["F"])) "not estimable" else sprintf("F = %.2f (%.0f, %.1f), p = %.3f", tr["F"], tr["df_num"], tr["df_den"], tr["p"]))
  }
}
T8 <- bind_rows(rows)
T8 <- bind_rows(T8, data.frame(Moderator = "Training sample size (n)", Level = "log(n), continuous (observed n only)",
  n_obs = nrow(dn), k_studies = n_distinct(dn$Study_ID), n_obs_RMSE = nrow(dnr), k_studies_RMSE = n_distinct(dnr$Study_ID),
  Test_R2 = sprintf("slope = %.3f (SE %.3f, df %.1f), p = %.3f", cn2$beta[2], cn2$SE[2], cn2$df_Satt[2], cn2$p_Satt[2]),
  Test_Z = sprintf("slope = %.3f (SE %.3f, df %.1f), p = %.3f", cn$beta[2], cn$SE[2], cn$df_Satt[2], cn$p_Satt[2]),
  Test_logRMSE = sprintf("slope = %.3f (SE %.3f, df %.1f), p = %.3f", cnr$beta[2], cnr$SE[2], cnr$df_Satt[2], cnr$p_Satt[2])))
export_table(T8, "TableS41_terrain_sample_size_reference", translate = FALSE)
print(T8[, c("Moderator", "Level", "n_obs", "k_studies", "R2", "R2_CI_lower", "R2_CI_upper", "R2_PI_lower", "R2_PI_upper", "RMSE_median", "Test_R2", "Test_Z", "Test_logRMSE")])

gp <- function(x) case_when(x %in% c("GEDI L2A", "GEDI L2B", "GEDI L4A", "GEDI L4B") ~ x,
                            x == "GEDI Multi-level" ~ "Several GEDI levels", TRUE ~ "Other / not reported")
G1 <- A %>% filter(GEDI_used == "Yes") %>% mutate(Product = gp(GEDI_product))
R1 <- rmA %>% filter(GEDI_used == "Yes") %>% mutate(Product = gp(GEDI_product))
TS14 <- G1 %>% group_by(Product) %>% summarise(n_obs = n(), k_studies = n_distinct(Study_ID),
          R2_mean = mean(R2), R2_median = median(R2), R2_min = min(R2), R2_max = max(R2), .groups = "drop") %>%
  left_join(R1 %>% group_by(Product) %>% summarise(n_obs_RMSE = n(), RMSE_median = median(as.numeric(RMSE)), .groups = "drop"), by = "Product") %>%
  mutate(Note = ifelse(k_studies < 5, "< 5 studies: descriptive only", ""))
print(TS14)

role_en <- c("Prédicteur" = "Predictor", "Référence" = "Reference", "Référence probable" = "Probable reference", "Produit évalué" = "Evaluated product",
             "Non vérifié" = "Not verified", "Sans objet" = "GEDI not used", "Simulé" = "Simulated GEDI")
TS15 <- fz %>% mutate(Role = dplyr::recode(GEDI_role, !!!role_en), Analysed = TRUE) %>%
  group_by(Role) %>% summarise(n_obs_database = n(), n_obs_analysed = sum(Analysed), k_studies = n_distinct(Study_ID),
                               R2_mean_analysed = ifelse(any(Analysed), mean(as.numeric(R2)[Analysed]), NA), .groups = "drop")
print(TS15)

pf <- T8 %>% filter(!is.na(R2)) %>% mutate(Moderator = factor(sub("Training sample size \\(n\\)", "Training\nsize (n)", sub("Reference data", "Reference\ndata", Moderator)), c("Terrain", "Training\nsize (n)", "Reference\ndata")),
                                           Level = factor(paste0(Level, "  (k = ", k_studies, ")"), rev(paste0(Level, "  (k = ", k_studies, ")"))))
p <- ggplot(pf, aes(y = Level)) +
  geom_vline(xintercept = pooled(A, "Primary")$R2, linetype = 2, colour = "grey40") +
  geom_errorbarh(aes(xmin = R2_PI_untransformed_lower, xmax = R2_PI_untransformed_upper), height = 0, linewidth = 0.6, colour = "grey60") +
  geom_errorbarh(aes(xmin = R2_PI_lower, xmax = R2_PI_upper), height = 0, linewidth = 0.35, colour = "grey35", linetype = "dotted", position = position_nudge(y = -.25)) +
  geom_errorbarh(aes(xmin = R2_CI_lower, xmax = R2_CI_upper), height = 0, linewidth = 1.6, colour = "#1f5fa8") +
  geom_point(aes(x = R2, shape = Reliability), size = 3.4, colour = "#1f5fa8", fill = "white") +
  scale_shape_manual(values = c("≥ 5 studies" = 16, "< 5 studies: indicative" = 21), name = NULL) +
  facet_grid(Moderator ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  labs(x = "Pooled R\u00b2 (thick: 95% CI; thin: 95% PI, untransformed;\ndotted: 95% PI, Fisher's Z)", y = NULL,
       title = "Pooled R² by terrain, training size and reference data",
       subtitle = "Three-level model, CR2 inference; dashed line: overall pooled R²") +
  theme_bw(base_size = 15) + theme(legend.position = if (all(pf$k_studies >= 5)) "none" else "bottom", plot.title.position = "plot", strip.text.y = element_text(angle = 0, face = "bold"),
                                   plot.title = element_text(face = "bold"))
ggsave("output/figures/FigS05_terrain_sample_size_reference.png", p, width = 9, height = 7.2, dpi = 1000, bg = "white")
ggsave("output/figures/FigS05_terrain_sample_size_reference.pdf", p, width = 9, height = 7.2, device = grDevices::cairo_pdf, bg = "white")
cat("12_terrain_samplesize_gedi.R: done\n")
