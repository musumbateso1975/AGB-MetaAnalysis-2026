source("R/00d_prepare_analysis.R", encoding = "UTF-8")

m0 <- fit_ml(A); P1 <- pooled(A, "Principal (plafond n = 1000, CR2)")
i2 <- i2_ml(m0, A)
wcsv(cbind(P1, t(round(i2, 2)), mean_vi = round(mean(A$vi), 4)), "R01_modele_global.csv"); print(P1[, c("k_obs", "Etudes", "Z", "R2", "R2_IC_inf", "R2_IC_sup", "R2_IP_inf", "R2_IP_sup", "I2_total")])

mm <- fit_ml(A, mods = ~ Sensor + Model + Scale + Climate)
ct <- coef_test(mm, vcov = "CR2", cluster = A$Cluster)
ci <- conf_int(mm, vcov = "CR2", cluster = A$Cluster)
reg <- data.frame(Terme = rownames(mm$beta), beta = ct$beta, SE_CR2 = ct$SE, t = ct$tstat, ddl = ct$df_Satt, p = ct$p_Satt,
                  IC_inf = ci$CI_L, IC_sup = ci$CI_U)
wald <- function(pat, nom) { idx <- grep(pat, rownames(mm$beta))
  w <- Wald_test(mm, constraints = constrain_zero(idx), vcov = "CR2", cluster = A$Cluster, test = "HTZ")
  data.frame(Moderateur = nom, q = length(idx), F = w$Fstat, ddl_num = w$df_num, ddl_den = w$df_denom, p = w$p_val) }
wt <- bind_rows(wald("^Sensor", "Combinaison de capteurs"), wald("^Model", "Catégorie de modèle"),
                wald("^Scale", "Échelle"), wald("^Climate", "Zone climatique"))
R2meta <- 1 - sum(mm$sigma2) / sum(m0$sigma2)
wcsv(reg, "R02_meta_regression_coefficients.csv"); wcsv(wt, "R02b_tests_Wald_robustes.csv")
wcsv(data.frame(R2_meta = R2meta, tau2_nul = sum(m0$sigma2), tau2_residuel = sum(mm$sigma2), k = mm$k, Etudes = n_distinct(A$Study_ID[mm$not.na])), "R02c_variance_expliquee.csv")
print(wt); cat("R2_meta =", round(R2meta, 3), "\n")

subg <- function(var) {
  lv <- levels(factor(A[[var]])); out <- list()
  for (l in lv) { d <- A[A[[var]] == l & !is.na(A[[var]]), ]
    if (n_distinct(d$Study_ID) < 2) next
    out[[l]] <- cbind(Variable = var, Niveau = l, pooled(d, l, robust = n_distinct(d$Cluster) >= 4)) }
  bind_rows(out)
}
SG <- bind_rows(subg("Sensor"), subg("Model"), subg("Scale"), subg("Climate"), subg("LiDAR"), subg("GEDI_clean"))
SG$Fiabilite <- ifelse(SG$Etudes >= 5, "≥ 5 études", "< 5 études : indicatif")
wcsv(SG, "R03_sous_groupes_IC_IP.csv")

G <- A %>% filter(LiDAR == "Avec LiDAR", !is.na(GEDI_clean))
mg <- fit_ml(G, mods = ~ GEDI_clean); ctg <- coef_test(mg, vcov = "CR2", cluster = G$Cluster)
wcsv(data.frame(Comparaison = "Sans GEDI (autre LiDAR) vs GEDI, combinaisons avec LiDAR", k_obs = nrow(G),
                Etudes_GEDI = n_distinct(G$Study_ID[G$GEDI_clean == "GEDI"]), Etudes_sans = n_distinct(G$Study_ID[G$GEDI_clean == "Sans GEDI"]),
                beta = ctg$beta[2], SE_CR2 = ctg$SE[2], ddl = ctg$df_Satt[2], p = ctg$p_Satt[2]), "R03b_GEDI_exploratoire.csv")

rmA <- prep(rm0) %>% filter(RMSE_unit == "Mg_ha") %>% mutate(lR = log(as.numeric(RMSE)))
contrast <- function(var, Aset, Bset, nom, dz = A, dr = rmA) {
  z <- dz %>% filter(.data[[var]] %in% c(Aset, Bset), !is.na(.data[[var]])) %>% mutate(g = ifelse(.data[[var]] %in% Aset, "A", "B")) %>%
    group_by(Study_ID) %>% filter(n_distinct(g) == 2) %>%
    summarise(dZ = mean(Z[g == "A"]) - mean(Z[g == "B"]), dR2 = mean(R2[g == "A"]) - mean(R2[g == "B"]),
              vi = mean(vi[g == "A"]) / sum(g == "A") + mean(vi[g == "B"]) / sum(g == "B"), kA = sum(g == "A"), kB = sum(g == "B"), .groups = "drop")
  r <- dr %>% filter(.data[[var]] %in% c(Aset, Bset), !is.na(.data[[var]])) %>% mutate(g = ifelse(.data[[var]] %in% Aset, "A", "B")) %>%
    group_by(Study_ID) %>% filter(n_distinct(g) == 2) %>%
    summarise(lnRR = mean(lR[g == "A"]) - mean(lR[g == "B"]), .groups = "drop")
  tt <- function(x) if (length(x) >= 3) { t <- t.test(x); c(mean(x), t$conf.int) } else c(if (length(x)) mean(x) else NA, NA, NA)
  rz <- if (nrow(z) >= 3) tryCatch(rma(dZ, vi, data = z, method = "REML", test = "knha", control = list(stepadj = 0.5, maxiter = 1000)), error = function(e) rma(dZ, vi, data = z, method = "DL", test = "knha")) else NULL
  q <- tt(z$dR2); l <- tt(r$lnRR)
  list(etudes_z = z %>% mutate(Contraste = nom), etudes_r = r %>% mutate(Contraste = nom),
       resume = data.frame(Contraste = nom, Etudes_Z = nrow(z), dZ_RE = if (!is.null(rz)) as.numeric(coef(rz)) else NA,
         dZ_IC_inf = if (!is.null(rz)) rz$ci.lb else NA, dZ_IC_sup = if (!is.null(rz)) rz$ci.ub else NA, dZ_p = if (!is.null(rz)) rz$pval else NA,
         dR2_moy = q[1], dR2_IC_inf = q[2], dR2_IC_sup = q[3], Etudes_A_meilleur_R2 = sum(z$dR2 > 0),
         Etudes_RMSE = nrow(r), RMSE_var_pct = 100 * (exp(l[1]) - 1), RMSE_IC_inf_pct = 100 * (exp(l[2]) - 1), RMSE_IC_sup_pct = 100 * (exp(l[3]) - 1),
         Etudes_A_RMSE_plus_faible = sum(r$lnRR < 0)))
}
CT <- list(contrast("Fusion", "Multi-capteurs", "Mono-capteur", "Multi-capteurs vs mono-capteur"),
           contrast("LiDAR", "Avec LiDAR", "Sans LiDAR", "Avec LiDAR vs sans LiDAR"),
           contrast("Model_category", "Machine-learning", "Statistical", "Apprentissage automatique vs statistique"),
           contrast("Model_category", "Deep-learning", c("Machine-learning", "Statistical"), "Apprentissage profond vs autres"))
gl <- A %>% filter(LiDAR == "Avec LiDAR", !is.na(GEDI_clean)) %>% group_by(Study_ID) %>% filter(n_distinct(GEDI_clean) == 2) %>% pull(Study_ID) %>% unique()
cat("GEDI vs other LiDAR (LiDAR configurations):", length(gl), "study(ies):", paste(gl, collapse = ", "), "\n")
CTR <- bind_rows(lapply(CT, `[[`, "resume")); wcsv(CTR, "R04_contrastes_intra_etude.csv")
wcsv(bind_rows(lapply(CT, `[[`, "etudes_z")), "R04b_contrastes_par_etude_Z.csv"); wcsv(bind_rows(lapply(CT, `[[`, "etudes_r")), "R04c_contrastes_par_etude_RMSE.csv")
print(CTR[, c("Contraste", "Etudes_Z", "dZ_RE", "dZ_p", "dR2_moy", "Etudes_RMSE", "RMSE_var_pct")])

best1 <- A %>% group_by(Study_ID) %>% arrange(desc(Best_model == "Yes"), desc(Z)) %>% slice(1) %>% ungroup()
rand_once <- function() {
  d <- A %>% group_by(Study_ID) %>% slice_sample(n = 1) %>% ungroup()
  c(as.numeric(coef(rma(Z, vi, data = d, method = "REML", test = "knha"))),
    as.numeric(coef(rma(R2, rep(V_R2, nrow(d)), data = d, method = "REML", test = "knha"))))
}
rr <- replicate(500, rand_once()); rz <- rr[1, ]; rq <- rr[2, ]
SENS <- bind_rows(
  P1,
  pooled(A %>% filter(is.na(n) | n <= 5000), "Sans n > 5 000"),
  pooled(A %>% mutate(vi = make_vi(n, cap = Inf)), "Sans plafond de n"),
  pooled(A %>% mutate(vi = make_vi(n, cap = 500)), "Plafond n = 500"),
  pooled(A %>% filter(!n_imput), "Cas complets (sans n imputé)"),
  pooled(best1, "Une obs./étude : meilleur modèle (borne optimiste)", robust = TRUE),
  pooled(A %>% filter(Country != "China"), "Sans la Chine"),
  pooled(A %>% filter(Region != "East-Asia"), "Sans l'Asie de l'Est"),
  pooled(A %>% filter(Climate_zone == "Tropical"), "Tropical uniquement"),
  pooled(A %>% filter(Map_validation != "TRUE", GEDI_coding_doubtful != "TRUE"), "Sans validations de cartes ni GEDI suspect"),
  pooled(A %>% filter(GEDI_role != "Non vérifié"), "Sans rôle GEDI non vérifié"),
  pooled(A_all, "Avec les R² d'ajustement (training fit)"),
  pooled(A %>% filter(Simple_regression_R2 == "TRUE"), "Régressions linéaires simples (R² = r² défendable)"),
  pooled(A %>% filter(Region == "Central-Africa"), "Afrique centrale uniquement"))
SENS <- bind_rows(SENS, data.frame(Scenario = "Une obs./étude tirée au hasard (500 tirages, moyenne)", k_obs = n_distinct(A$Study_ID),
  Etudes = n_distinct(A$Study_ID), Z = mean(rz), IC_inf = quantile(rz, .025), IC_sup = quantile(rz, .975), R2 = mean(rq),
  R2_IC_inf = quantile(rq, .025), R2_IC_sup = quantile(rq, .975),
  R2_Zscale = r2(mean(rz)), R2_Zscale_IC_inf = r2(quantile(rz, .025)), R2_Zscale_IC_sup = r2(quantile(rz, .975)), Inference = "intervalle = 2,5–97,5 % des tirages"))
wcsv(SENS, "R05_sensibilites.csv"); print(SENS[, c("Scenario", "k_obs", "Etudes", "R2", "R2_IC_inf", "R2_IC_sup", "R2_IP_inf", "R2_IP_sup")])

loso <- bind_rows(lapply(unique(A$Cluster), function(cl) { d <- A[A$Cluster != cl, ]; m <- fit_ml(d)
  mr <- tryCatch(rma.mv(R2, rep(V_R2, nrow(d)), random = ~ 1 | Cluster / Study_ID / Obs_ID, data = d, method = "REML", sparse = TRUE), error = function(e) tryCatch(rma.mv(R2, rep(V_R2, nrow(d)), random = ~ 1 | Cluster / Study_ID / Obs_ID, data = d, method = "REML", sparse = TRUE, control = list(optimizer = "optim", optmethod = "BFGS")), error = function(e2) rma.mv(R2, rep(V_R2, nrow(d)), random = ~ 1 | Cluster / Study_ID / Obs_ID, data = d, method = "REML", sparse = TRUE, control = list(optimizer = "optim", optmethod = "Nelder-Mead"))))
  data.frame(Cluster_retire = cl, R2_untransformed = as.numeric(coef(mr)), Z = as.numeric(coef(m)), R2_Zscale = r2(coef(m))) }))
wcsv(loso, "R05b_LOSO_cluster.csv")

A$sei <- sqrt(A$vi); me <- fit_ml(A, mods = ~ sei); cte <- coef_test(me, vcov = "CR2", cluster = A$Cluster)
eg <- data.frame(Test = "Egger multiniveau (Z ~ SE, 3 niveaux, CR2)", pente = cte$beta[2], SE = cte$SE[2], ddl = cte$df_Satt[2], p = cte$p_Satt[2])
wcsv(eg, "R06_biais_publication.csv"); print(eg)

rd <- rmA %>% group_by(Sensor_combination, Study_ID) %>% summarise(med = median(as.numeric(RMSE)), .groups = "drop") %>%
  group_by(Combinaison = Sensor_combination) %>% summarise(Etudes = n(), RMSE_mediane_des_etudes = median(med),
  Q1 = quantile(med, .25), Q3 = quantile(med, .75), .groups = "drop") %>% arrange(RMSE_mediane_des_etudes)
wcsv(rd, "R07_RMSE_descriptif_Mgha.csv")

cell <- function(d, lab1, lab2) { ke <- n_distinct(d$Study_ID)
  if (ke < 5) return(data.frame(Contexte = lab1, Capteurs = lab2, k_obs = nrow(d), Etudes = ke, R2 = NA, R2_IC = "", R2_IP = "", Statut = "Preuves insuffisantes (< 5 études)"))
  p <- pooled(d, "", robust = n_distinct(d$Cluster) >= 4)
  data.frame(Contexte = lab1, Capteurs = lab2, k_obs = nrow(d), Etudes = ke, R2 = round(p$R2, 2),
             R2_IC = sprintf("%.2f–%.2f", p$R2_IC_inf, p$R2_IC_sup), R2_IP = sprintf("%.2f–%.2f", p$R2_IP_inf, p$R2_IP_sup), Statut = "Estimable") }
T6 <- list()
for (cz in levels(A$Climate)) for (li in c("Avec LiDAR", "Sans LiDAR")) T6[[length(T6) + 1]] <- cell(A %>% filter(Climate == cz, LiDAR == li), paste("Climat :", cz), li)
for (sc in levels(A$Scale)) for (li in c("Avec LiDAR", "Sans LiDAR")) T6[[length(T6) + 1]] <- cell(A %>% filter(Scale == sc, LiDAR == li), paste("Échelle :", sc), li)
T6 <- bind_rows(T6); wcsv(T6, "R08_nouvelle_Table6.csv"); print(T6)

saveRDS(list(A = A, SG = SG, CT = CT, CTR = CTR, SENS = SENS, P1 = P1, m0 = m0, rmA = rmA, fz = fz, rm = prep(rm0), fz_all = prep(fz0_all), rm_all = prep(rm0_all)), file.path(CACHE_DIR, "reanalyse.rds"))

cat("OK\n")
