# Model diagnostics of the primary three-level model (untransformed R2 and Fisher's Z):
# standardized residuals, cluster-level externally studentized residuals and Cook's distances,
# normality of residuals and of cluster random effects, profile-likelihood CIs of the variance components.
# Output: output/tables/TableS35_model_diagnostics.csv, TableS35_cluster_influence.csv, output/figures/FigS08_model_diagnostics.*
source("R/00d_prepare_analysis.R", encoding = "UTF-8")
source("R/00b_figure_style.R", encoding = "UTF-8")
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

d <- A
fit_r2 <- function(ctrl = list()) rma.mv(R2, rep(V_R2, nrow(d)), random = ~ 1 | Cluster / Study_ID / Obs_ID,
                                         data = d, method = "REML", sparse = TRUE, control = ctrl)
mods <- list(
  "Untransformed R²" = tryCatch(fit_r2(), error = function(e) fit_r2(list(optimizer = "optim", optmethod = "Nelder-Mead"))),
  "Fisher's Z" = fit_ml(d))

cook_crit <- qchisq(0.5, df = 1)   # usual reference value for one fixed coefficient
diag_rows <- list(); infl_rows <- list(); qq_rows <- list(); prof_rows <- list()
for (sc in names(mods)) {
  m <- mods[[sc]]
  rs <- rstandard(m)$z
  sw <- shapiro.test(rs)
  rst <- rstudent(m, cluster = d$Cluster, reestimate = TRUE)
  cd  <- cooks.distance(m, cluster = d$Cluster, reestimate = TRUE)
  rc  <- data.frame(Cluster = names(cd), cooks_d = as.numeric(cd))
  ix <- match(rc$Cluster, rst$cluster$slab)
  rc$X2 <- as.numeric(rst$cluster$X2[ix]); rc$k <- as.numeric(rst$cluster$k[ix])
  rc$p_outlier <- pchisq(rc$X2, df = rc$k, lower.tail = FALSE)
  re  <- ranef(m)$Cluster
  ci  <- confint(m)
  ci_s <- function(i) { x <- ci[[i]]$random[1, ]; sprintf("%.4f (%.4f–%.4f)", x[1], x[2], x[3]) }
  top <- rc[which.max(rc$cooks_d), ]
  diag_rows[[sc]] <- data.frame(
    Scale = sc, k_obs = m$k, Clusters = nrow(rc),
    Obs_abs_std_resid_gt_1.96 = sum(abs(rs) > 1.96), Obs_abs_std_resid_gt_3 = sum(abs(rs) > 3),
    Shapiro_W_std_resid = round(unname(sw$statistic), 3), Shapiro_p_std_resid = signif(sw$p.value, 3),
    Shapiro_p_cluster_effects = signif(shapiro.test(re$intrcpt)$p.value, 3),
    Clusters_outlier_p_lt_0.05 = sum(rc$p_outlier < 0.05, na.rm = TRUE),
    Clusters_outlier_p_lt_Bonferroni = sum(rc$p_outlier < 0.05 / nrow(rc), na.rm = TRUE),
    Max_cooks_d = round(top$cooks_d, 3), Max_cooks_cluster = top$Cluster,
    Clusters_cooks_gt_chisq50 = sum(rc$cooks_d > cook_crit, na.rm = TRUE),
    sigma2_cluster_profileCI = ci_s(1), sigma2_study_profileCI = ci_s(2), sigma2_obs_profileCI = ci_s(3),
    check.names = FALSE)
  infl_rows[[sc]] <- data.frame(Scale = sc, rc[order(-rc$cooks_d), ])
  qq_rows[[sc]] <- rbind(
    data.frame(Scale = sc, Type = "Standardized residuals (observations)", z = sort(rs), q = qnorm(ppoints(length(rs)))),
    data.frame(Scale = sc, Type = "Cluster random effects", z = sort(re$intrcpt / sd(re$intrcpt)), q = qnorm(ppoints(nrow(re)))))
  if (sc == "Untransformed R²") {
    for (i in 1:3) {
      ub <- max(ci[[i]]$random[1, 3], m$sigma2[i]) * 1.5
      pr <- profile(m, sigma2 = i, xlim = c(0, ub), plot = FALSE, progbar = FALSE, steps = 25)
      prof_rows[[i]] <- data.frame(Component = c("σ² cluster", "σ² study", "σ² observation")[i], s2 = pr$sigma2, ll = pr$ll, est = m$sigma2[i])
    }
  }
}
tab <- bind_rows(diag_rows); infl <- bind_rows(infl_rows)
export_table(tab, "TableS35_model_diagnostics", translate = FALSE)
export_table(infl, "TableS35_cluster_influence", translate = FALSE)

qq <- bind_rows(qq_rows) %>% mutate(Scale = factor(Scale, names(mods)),
  Type = factor(Type, c("Standardized residuals (observations)", "Cluster random effects"))); infl_r2 <- infl %>% filter(Scale == "Untransformed R²") %>% mutate(lab = ifelse(cooks_d > cook_crit, Cluster, ""))
pa <- ggplot(qq, aes(q, z)) + geom_abline(colour = MUTED, linetype = "dashed") + geom_point(size = .9, colour = INK, alpha = .7) +
  facet_grid(Scale ~ Type) + labs(title = "(a) Normal Q–Q plots", x = "Theoretical quantile", y = "Standardized value") + theme_pub(9.5)
pb <- ggplot(infl_r2, aes(reorder(Cluster, -cooks_d), cooks_d)) + geom_col(fill = "#4c72b0", width = .7) +
  geom_hline(yintercept = cook_crit, colour = REDL, linetype = "dashed") +
  labs(title = "(b) Cook's distance by cluster (untransformed R²)", x = "Cluster (ordered)", y = "Cook's distance") +
  theme_pub(9.5) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
pc <- ggplot(bind_rows(prof_rows), aes(s2, ll)) + geom_line(colour = INK) + geom_vline(aes(xintercept = est), colour = REDL, linetype = "dashed") +
  facet_wrap(~ factor(Component, c("σ² cluster", "σ² study", "σ² observation")), scales = "free_x") +
  scale_x_continuous(n.breaks = 3) + labs(title = "(c) Profile likelihood of the variance components (untransformed R²)",
  x = "Variance component", y = "Restricted log-likelihood") + theme_pub(9.5)
g <- pa / pb / pc + plot_layout(heights = c(2, 1, 1))
save_pub(g, "FigS08_model_diagnostics", h = 9, w = FIG_W)
print(tab)
