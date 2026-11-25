suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(maps); library(patchwork); library(readr) })
OUT <- CACHE_DIR
FIG <- "output/figures"; dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
R <- readRDS(file.path(CACHE_DIR, "reanalyse.rds")); A <- R$A
BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; INK <- "#1f2328"; MUTED <- "#6b7178"; GRID <- "#e6e8eb"
th <- theme_minimal(base_size = 11, base_family = "sans") + theme(
  panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = GRID, linewidth = .3),
  axis.text = element_text(colour = "#2b3035", size = 10.5), axis.title = element_text(colour = INK, size = 11), plot.title = element_text(face = "bold", colour = INK, size = 12.5),
  plot.subtitle = element_text(colour = "#4b5157", size = 10), legend.text = element_text(size = 10.5), legend.position = "top", strip.text = element_text(face = "bold", colour = INK, hjust = 0, size = 11.5), plot.title.position = "plot")
save <- function(p, f, w, h, h1 = NULL, w1 = NULL) { ggsave(file.path(FIG, paste0(f, ".png")), p, width = w, height = h, dpi = 1000, bg = "white"); ggsave(file.path(FIG, paste0(f, ".pdf")), p, width = w, height = h, device = grDevices::cairo_pdf, bg = "white"); export_columns(p, f, w, h, h1, FIG, w1) }

if (!exists("removed_by_study")) source("R/01_data_import.R")
auF <- dplyr::filter(removed_by_study, Dataset == "Fisher_Z"); auR <- dplyr::filter(removed_by_study, Dataset == "RMSE")
nAn <- nrow(A); kAn <- n_distinct(A$Study_ID)
nFZ <- nrow(R$fz) + sum(removed_summary$n_obs[removed_summary$Dataset == "Fisher_Z"]); kFZ <- n_distinct(c(R$fz$Study_ID, auF$Study_ID))
nRM <- nrow(R$rm) + sum(removed_summary$n_obs[removed_summary$Dataset == "RMSE"]); kRM <- n_distinct(c(R$rm$Study_ID, auR$Study_ID)); nRMc <- sum(R$rm$RMSE_unit == "Mg_ha")
nEx <- sum(auF$n_obs[!duplicated(auF[c("Dataset", "Reason")])]); kEx <- n_distinct(auF$Study_ID)
kInc <- n_distinct(c(R$fz_all$Study_ID, R$rm_all$Study_ID))
X_EXTR <- c(gedi_ref = 5L, metric = 5L, height = 1L, gedi_product = 1L)
kAll <- kInc + sum(X_EXTR)
C_MAIN <- "#e8f1fc"; C_MAIN_B <- "#2a78d6"; C_EXC <- "#f2f3f5"; C_EXC_B <- "#8a9097"; C_RES <- "#d6e7fb"; C_RES_B <- "#123e73"
N_ID <- 5412L; N_DUP <- 2193L; N_SCR <- N_ID - N_DUP; N_EXS <- 2997L; N_SOUGHT <- N_SCR - N_EXS; N_NR <- 2L
N_FT <- N_SOUGHT - N_NR; N_EXF <- 128L; N_ELIG <- N_FT - N_EXF; N_ST <- kAll
stopifnot(N_SOUGHT == 222L, N_ELIG == 92L, N_ST - sum(X_EXTR) == kInc)
box <- function(x, y, w, h, txt, fill = C_MAIN, col = C_MAIN_B, size = 3.6, face = "plain", hjust = .5) list(
  annotation_custom(grid::roundrectGrob(r = unit(2, "mm"), gp = grid::gpar(fill = fill, col = col, lwd = 1.2)),
                    xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2),
  annotate("text", x = if (hjust == 0) x - w / 2 + .18 else x, y = y, label = txt, size = size, fontface = face, colour = INK, lineheight = 1.05, hjust = hjust))
arr <- function(x1, y1, x2, y2) annotate("segment", x = x1, y = y1, xend = x2, yend = y2, linewidth = .6, colour = "#3b4046",
  arrow = arrow(length = unit(2, "mm"), type = "closed"))
band <- function(y1, y2, lab) list(
  annotation_custom(grid::roundrectGrob(r = unit(2, "mm"), gp = grid::gpar(fill = "#1f5fae", col = NA)), xmin = -1.05, xmax = -.35, ymin = y1, ymax = y2),
  annotate("text", x = -.7, y = (y1 + y2) / 2, label = lab, angle = 90, fontface = "bold", size = 4.1, colour = "white"))
fmt <- function(x) format(x, big.mark = ",")
XL <- 2.6; WL <- 5.0; XR <- 8.15; WR <- 3.9
p1 <- ggplot() + coord_cartesian(xlim = c(-1.1, 10.2), ylim = c(1.9, 13.7), expand = FALSE) + theme_void() +
  band(11.6, 13.6, "Identification") + band(7.85, 11.45, "Screening") + band(5.95, 7.7, "Eligibility") + band(1.95, 5.8, "Included") +
  arr(XL + WL / 2, 12.6, XR - WR / 2, 12.6) + arr(XL, 11.75, XL, 11.25) + arr(XL + WL / 2, 10.75, XR - WR / 2, 10.75) +
  arr(XL, 10.25, XL, 9.4) + arr(XL + WL / 2, 8.9, XR - WR / 2, 8.9) + arr(XL, 8.4, XL, 7.4) +
  arr(XL + WL / 2, 6.9, XR - WR / 2, 6.9) + arr(XL, 6.4, XL, 5.675) +
  arr(XL - 1.0, 4.725, 1.25, 3.875) + arr(XL + 1.0, 4.725, 4.15, 3.875) +
  box(XL, 12.6, WL, 1.7, size = 3.4, sprintf("Records identified\n(Web of Science, Scopus, Google Scholar\nand citation searching)\nn = %s", fmt(N_ID))) +
  box(XR, 12.6, WR, 1.0, sprintf("Duplicates removed\nn = %s", fmt(N_DUP)), fill = C_EXC, col = C_EXC_B) +
  box(XL, 10.75, WL, 1.0, sprintf("Records screened (title and abstract)\nn = %s", fmt(N_SCR))) +
  box(XR, 10.75, WR, 1.25, sprintf("Records excluded, n = %s\n• Title and abstract (%s)\n• Full text not retrieved (%d)", fmt(N_EXS + N_NR), fmt(N_EXS), N_NR), fill = C_EXC, col = C_EXC_B, size = 3.3, hjust = 0) +
  box(XL, 8.9, WL, 1.0, sprintf("Full-text reports assessed for eligibility\nn = %d", N_FT)) +
  box(XR, 8.85, WR, 1.8, sprintf("Reports excluded, n = %d\n• No forest AGB estimate\n• No usable R² or RMSE\n• No validation against\n   reference data\n• Review, simulation or\n   non-forest study", N_EXF), fill = C_EXC, col = C_EXC_B, size = 3.2, hjust = 0) +
  box(XL, 6.9, WL, 1.0, sprintf("Eligible reports, n = %d\n(%d studies)", N_ELIG, N_ST)) +
  box(XR, 6.85, WR, 1.55, sprintf("Excluded at data extraction, n = %d\n• Accuracy vs GEDI L4A/L4B (%d)\n• No comparable R² or RMSE (%d)\n• Canopy height only (%d)\n• GEDI product validation only (%d)",
    sum(X_EXTR), X_EXTR[["gedi_ref"]], X_EXTR[["metric"]], X_EXTR[["height"]], X_EXTR[["gedi_product"]]), fill = C_EXC, col = C_EXC_B, size = 3.2, hjust = 0) +
  box(XL, 5.2, WL, .95, sprintf("Studies included in the review\nn = %d", kInc), fill = C_RES, col = C_RES_B, face = "bold") +
  box(1.25, 2.95, 2.75, 1.85, sprintf("R² dataset\n%d values extracted\n(%d studies)\n%d analysed\n(%d studies)", nFZ, kFZ, nAn, kAn), fill = C_RES, col = C_RES_B, size = 3.3) +
  box(4.15, 2.95, 2.75, 1.85, sprintf("RMSE dataset\n%d values extracted\n(%d studies)\n%d analysed in Mg/ha\n(%d studies)", nRM, kRM, nrow(R$rm), n_distinct(R$rm$Study_ID)), fill = C_RES, col = C_RES_B, size = 3.3)
save(p1, "Fig1_PRISMA_flow_diagram", 6.5, 7.4)

source("R/11b_map.R")
carte_etudes(R, FIG)

sg <- R$SG %>% filter(Variable == "Sensor") %>% mutate(Niveau = factor(Niveau, rev(levels(A$Sensor))),
  lab = sprintf("%d obs.\n%d studies", k_obs, Etudes), faible = Etudes < 5)
p3 <- ggplot(sg, aes(y = Niveau)) +
  geom_vline(xintercept = R$P1$R2, linetype = "dashed", colour = MUTED, linewidth = .3) +
  geom_errorbarh(aes(xmin = R2_IPb_inf, xmax = R2_IPb_sup), height = 0, colour = "#9fb9db", linewidth = .6) +
  geom_errorbarh(aes(xmin = R2_IP_inf, xmax = R2_IP_sup), height = 0, colour = MUTED, linewidth = .35, linetype = "dotted", position = position_nudge(y = -.22)) +
  geom_errorbarh(aes(xmin = R2_IC_inf, xmax = R2_IC_sup), height = 0, colour = BLUE, linewidth = 2) +
  geom_point(aes(x = R2, shape = faible), size = 3, colour = BLUE, fill = "white", stroke = .9) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21), labels = c("≥ 5 studies", "< 5 studies (indicative)"), name = NULL) +
  geom_text(aes(x = 1.03, label = lab), hjust = 0, size = 3.5, colour = "#2b3035", lineheight = .9) +
  scale_x_continuous(limits = c(0, 1.25), breaks = seq(0, 1, .2), name = expression("Pooled R"^2)) +
  labs(y = NULL, title = "Reported accuracy by sensor configuration",
       subtitle = "Thick bar: 95% CI (CR2 robust). Thin bar: 95% prediction interval, untransformed R².\nDotted bar: prediction interval on Fisher's Z (complementary). Dashed line: overall estimate.") + th
save(p3, "Fig4_sensor_configuration_CI_PI", 6.5, 5, h1 = 5, w1 = 5.3)

EN <- c("Multi-capteurs vs mono-capteur" = "Multi-sensor vs single-sensor", "Avec LiDAR vs sans LiDAR" = "With LiDAR vs without LiDAR",
  "Apprentissage automatique vs statistique" = "Machine learning vs statistical", "Apprentissage profond vs autres" = "Deep learning vs other models",
  "GEDI vs sans GEDI (non suspect ; descriptif)" = "GEDI vs other LiDAR (descriptive)")
R$CTR$Contraste <- EN[R$CTR$Contraste]; for (i in seq_along(R$CT)) { R$CT[[i]]$etudes_z$Contraste <- EN[R$CT[[i]]$etudes_z$Contraste]; R$CT[[i]]$etudes_r$Contraste <- EN[R$CT[[i]]$etudes_r$Contraste] }
ctr <- R$CTR %>% mutate(Contraste = factor(Contraste, rev(Contraste)), desc = grepl("descriptif", Contraste))
ez <- bind_rows(lapply(R$CT, `[[`, "etudes_z")) %>% mutate(Contraste = factor(Contraste, levels(ctr$Contraste)))
er <- bind_rows(lapply(R$CT, `[[`, "etudes_r")) %>% mutate(Contraste = factor(Contraste, levels(ctr$Contraste)), pct = 100 * (exp(lnRR) - 1))
thL <- theme_minimal(base_size = 13, base_family = "sans") + theme(
  panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), panel.grid.major.x = element_line(colour = GRID, linewidth = .35),
  axis.text = element_text(colour = INK, size = 12), axis.title = element_text(colour = INK, size = 12.5), plot.title = element_text(face = "bold", colour = INK, size = 13.5),
  plot.subtitle = element_text(colour = "#4b5157", size = 11.5), legend.position = "bottom", legend.text = element_text(size = 12), legend.title = element_blank(),
  plot.title.position = "plot", plot.margin = margin(6, 8, 6, 6))
ctr$kz <- sprintf("k = %d", ctr$Etudes_Z); ctr$kr <- sprintf("k = %d", ctr$Etudes_RMSE)
pb <- ggplot(ctr, aes(y = as.numeric(Contraste))) + geom_vline(xintercept = 0, colour = INK, linewidth = .4) +
  geom_point(data = er, aes(x = pct, y = as.numeric(Contraste)), colour = "#f5b59b", size = 2.2, position = position_jitter(height = .12, seed = 1)) +
  geom_errorbarh(aes(xmin = RMSE_IC_inf_pct, xmax = RMSE_IC_sup_pct), height = 0, colour = ORANGE, linewidth = 2.2) +
  geom_point(aes(x = RMSE_var_pct), colour = ORANGE, size = 4.2) +
  scale_y_continuous(breaks = seq_len(nlevels(ctr$Contraste)), labels = levels(ctr$Contraste), expand = expansion(add = .6),
                     sec.axis = sec_axis(~ ., breaks = as.numeric(ctr$Contraste), labels = ctr$kr)) +
  scale_x_continuous(limits = c(-70, 60), name = "RMSE change within study (%, A vs B)") +
  labs(y = NULL, title = "b  RMSE (Mg/ha)") + thL
pa <- ggplot(ctr, aes(y = as.numeric(Contraste))) + geom_vline(xintercept = 0, colour = INK, linewidth = .4) +
  geom_point(data = ez, aes(x = dR2, y = as.numeric(Contraste)), colour = "#9fb9db", size = 2.2, position = position_jitter(height = .12, seed = 1)) +
  geom_errorbarh(aes(xmin = dR2_IC_inf, xmax = dR2_IC_sup), height = 0, colour = BLUE, linewidth = 2.2) +
  geom_point(aes(x = dR2_moy), colour = BLUE, size = 4.2) +
  scale_y_continuous(breaks = seq_len(nlevels(ctr$Contraste)), labels = levels(ctr$Contraste), expand = expansion(add = .6),
                     sec.axis = sec_axis(~ ., breaks = as.numeric(ctr$Contraste), labels = ctr$kz)) +
  scale_x_continuous(limits = c(-.35, .55), breaks = seq(-.2, .4, .2), name = expression(Delta*"R"^2*" within study (A − B)")) +
  labs(y = NULL, title = "a  Coefficient of determination (R²)") + thL
p4 <- (pa / pb) + plot_annotation(title = "Within-study contrasts (same plots, same study)",
  subtitle = "Large point: mean across studies with 95% CI; small points: individual studies; k: number of studies.",
  theme = theme(plot.title = element_text(face = "bold", size = 14), plot.subtitle = element_text(colour = "#4b5157", size = 11.5)))
save(p4, "Fig6_within_study_contrasts", 8, 9)

SEN <- c("Manuscrit soumis (314 lignes avec copies ; 2 niveaux ; sans RVE)" = "Submitted analysis (314 rows incl. duplicates; 2 levels; no RVE)",
  "Principal (plafond n = 1000, CR2)" = "Primary (n capped at 1000; CR2)", "Sans n > 5 000" = "Excluding n > 5,000", "Sans plafond de n" = "No cap on n",
  "Plafond n = 500" = "n capped at 500", "Cas complets (sans n imputé)" = "Complete cases (no imputed n)",
  "Une obs./étude : meilleur modèle (borne optimiste)" = "One obs./study: best model (optimistic bound)", "Sans la Chine" = "Excluding China",
  "Sans l'Asie de l'Est" = "Excluding East Asia", "Tropical uniquement" = "Tropical only",
  "Sans validations de cartes ni GEDI suspect" = "Excluding map validations and doubtful GEDI coding",
  "Régressions linéaires simples (R² = r² défendable)" = "Simple linear regressions only", "Afrique centrale uniquement" = "Central Africa only", "Sans rôle GEDI non vérifié" = "Excluding unverified GEDI role", "Avec les R² d'ajustement (training fit)" = "Including training-fit R²",
  "Une obs./étude tirée au hasard (500 tirages, moyenne)" = "One obs./study at random (500 draws)")
R$SENS$Scenario <- SEN[R$SENS$Scenario]
se <- R$SENS %>% mutate(prim = grepl("^Primary", Scenario), lab = sprintf("%d obs. · %d st.", k_obs, Etudes)) %>%
  arrange(desc(prim)) %>% mutate(pos = rev(seq_len(n())))
p5 <- ggplot(se, aes(y = pos)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = max(se$pos) - .5, ymax = max(se$pos) + .5, fill = "#fdeee6") +
  geom_vline(xintercept = R$P1$R2, linetype = "dashed", colour = MUTED, linewidth = .4) +
  geom_errorbarh(aes(xmin = R2_IC_inf, xmax = R2_IC_sup, colour = prim), height = 0, linewidth = 1.5) +
  geom_point(aes(x = R2, colour = prim), size = 3.4) + scale_colour_manual(values = c(`FALSE` = BLUE, `TRUE` = ORANGE), guide = "none") +
  scale_y_continuous(breaks = se$pos, labels = se$Scenario, expand = expansion(add = .6), sec.axis = sec_axis(~ ., breaks = se$pos, labels = se$lab)) +
  scale_x_continuous(limits = c(.35, .9), breaks = seq(.4, .9, .1), name = "Pooled R² with 95% CI") +
  labs(y = NULL, title = "Sensitivity of the pooled estimate",
       subtitle = sprintf("Orange: primary analysis; dashed line: primary estimate (R² = %.2f). st. = studies", R$P1$R2)) +
  thL + theme(axis.text.y.left = element_text(size = 11.5), axis.text.y.right = element_text(size = 11, colour = "#4b5157"))
save(p5, "Fig7_sensitivity_analyses", 9.5, 8)

mu <- R$P1$Z; seg <- seq(0, max(sqrt(A$vi)) * 1.05, length.out = 50)
p6 <- ggplot(A, aes(Z, sqrt(vi))) +
  geom_line(data = data.frame(s = seg), aes(mu - 1.96 * s, s), inherit.aes = FALSE, colour = MUTED, linewidth = .3, linetype = "dotted") +
  geom_line(data = data.frame(s = seg), aes(mu + 1.96 * s, s), inherit.aes = FALSE, colour = MUTED, linewidth = .3, linetype = "dotted") +
  geom_vline(xintercept = mu, colour = ORANGE, linewidth = .4) + geom_point(colour = BLUE, alpha = .6, size = 1.7) +
  scale_y_reverse(name = "Standard error of Fisher's Z") + xlab("Fisher's Z") +
  labs(title = "Funnel plot", subtitle = sprintf("Multilevel Egger test (CR2): p = %.2f.\nSE from capped training n (≤ 1000); %d imputed.",
                                                   read_csv("output/tables/TableS33_multilevel_Egger.csv", show_col_types = FALSE)$p, sum(A$n_imput))) +
  th + theme(panel.grid.major.y = element_line(colour = GRID, linewidth = .3))
save(p6, "FigS07_funnel_plot", 6.5, 5.2)

cat("Figures OK\n")
