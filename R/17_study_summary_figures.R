suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(ggplot2); library(readr) })
TAB <- "output/tables"; FIG <- "output/figures"
dir.create(TAB, showWarnings = FALSE, recursive = TRUE); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
INK <- "#1f2328"; GRID <- "#e6e8eb"; BLUE <- "#2a78d6"; ORANGE <- "#eb6834"; GREY <- "#8a9097"
th <- theme_minimal(base_size = 10.5) + theme(panel.grid.minor = element_blank(), plot.title.position = "plot",
  plot.title = element_text(face = "bold", colour = INK), legend.position = "top")
sv <- function(p, f, w, h, h1 = NULL) { ggsave(file.path(FIG, paste0(f, ".png")), p, width = w, height = h, dpi = 1000, bg = "white"); ggsave(file.path(FIG, paste0(f, ".pdf")), p, width = w, height = h, device = grDevices::cairo_pdf, bg = "white"); export_columns(p, f, w, h, h1, FIG) }

m  <- read_csv("data/AGB_master_database.csv", show_col_types = FALSE)
fz <- read_csv("data/DATA_FisherZ.csv", show_col_types = FALSE)
rm <- read_csv("data/DATA_RMS.csv", show_col_types = FALSE)
source("R/00c_export_helpers.R")
A  <- fz[primary_r2(fz), ]
RM <- rm %>% filter(RMSE_unit == "Mg_ha", primary_rmse(rm))
stopifnot(nrow(A) == 258, nrow(RM) == 294)

EN_from <- c("Non v\u00e9rifi\u00e9", "Produit \u00e9valu\u00e9", "Pr\u00e9dicteur", "R\u00e9f\u00e9rence", "R\u00e9f\u00e9rence probable", "Simul\u00e9", "Auxiliary (sample expansion via ALS)", "Indirect predictor (canopy height trained on GEDI L2A RH100)", "Training reference (GEDI L4A); F3 simulated GEDI", "ALS (terrain seulement)", "Carte globale de hauteur", "GEDI simul\u00e9 (ALS)", "GLAS (produit)", "AGB d\u00e9riv\u00e9e de LiDAR", "GEDI probable (\u00e0 confirmer)", "Terrain (pr\u00e9sum\u00e9)", "Terrain (v\u00e9rifi\u00e9)", "Multiple", "Democratic Republic of Congo", "GEDI L4A (empreintes)", "Training reference (GEDI L4A)", "GEDI L4B (grille 1 km)", "Training reference (GEDI L4B)")
EN_to <- c("not verified", "evaluated product", "predictor", "reference", "probable reference", "simulated", "auxiliary (sample expansion via ALS)", "indirect predictor (GEDI L2A height model)", "training reference (L4A); simulated", "ALS (terrain only)", "global canopy-height map", "simulated GEDI (ALS)", "GLAS (product)", "LiDAR-derived AGB", "GEDI (probable)", "field (presumed)", "field (verified)", "Multi-country", "DR Congo", "GEDI L4A footprints", "training reference (L4A)", "GEDI L4B (1-km grid)", "training reference (L4B)")
tr <- function(x) { i <- match(enc2utf8(x), enc2utf8(EN_from)); ifelse(is.na(i), x, EN_to[i]) }
uq <- function(x) { x <- tr(as.character(x)); x <- unique(na.omit(as.character(x))); x <- x[!x %in% c("", "Not applicable", "Sans objet", "Sans LiDAR")]
                    if (length(x)) paste(sort(x), collapse = "; ") else "\u2014" }
rng <- function(x, d = 2) { x <- na.omit(x); if (!length(x)) "\u2014" else if (length(x) == 1 || min(x) == max(x))
                    formatC(x[1], format = "f", digits = d) else paste0(formatC(min(x), format = "f", digits = d), "\u2013", formatC(max(x), format = "f", digits = d)) }

ids <- sort(unique(c(fz$Study_ID, rm$Study_ID)))
keep <- c("Study_ID","Reference","Year","Country","Climate_zone","Forest_type","Study_scale","Sensor_combination","LiDAR_platform","GEDI_role","Model_algorithm","Reference_type","Training_samples")
sel <- function(d) d %>% select(all_of(keep)) %>% mutate(across(-c(Year, Training_samples), as.character), Year = as.integer(Year), Training_samples = suppressWarnings(as.numeric(Training_samples)))
info <- bind_rows(sel(fz), sel(rm)) %>% group_by(Study_ID) %>% summarise(
  Reference = first(Reference), Year = first(Year), Country = uq(Country), Climate = uq(Climate_zone),
  Forest_type = uq(Forest_type), Scale = uq(Study_scale), Sensors = uq(Sensor_combination),
  LiDAR_platform = uq(LiDAR_platform), GEDI_role = uq(GEDI_role), Algorithms = uq(Model_algorithm),
  Reference_data = uq(Reference_type), n_train = ifelse(all(is.na(Training_samples)), NA_real_, median(Training_samples, na.rm = TRUE)), .groups = "drop")
r2 <- A %>% group_by(Study_ID) %>% summarise(k_R2 = n(), R2_range = rng(R2), .groups = "drop")
best <- A %>% group_by(Study_ID) %>% slice_max(R2, n = 1, with_ties = FALSE) %>%
  transmute(Study_ID, Best_model = paste0(Model_algorithm, " (", Sensor_combination, "; R\u00b2 = ", formatC(R2, format = "f", digits = 2), ")"))
rmse <- RM %>% group_by(Study_ID) %>% summarise(k_RMSE = n(), RMSE_range = rng(RMSE, 1), .groups = "drop")
S <- info %>% left_join(r2, by = "Study_ID") %>% left_join(rmse, by = "Study_ID") %>% left_join(best, by = "Study_ID") %>%
  mutate(across(c(k_R2, k_RMSE), ~ replace_na(.x, 0L)), across(c(R2_range, RMSE_range, Best_model), ~ replace_na(.x, "\u2014")),
         n_train = ifelse(is.na(n_train), "\u2014", formatC(n_train, format = "f", digits = 0))) %>% arrange(Year, Reference)
stopifnot(nrow(S) == 72)
source("R/00c_export_helpers.R"); export_table(S, "TableS06_study_characteristics", translate = FALSE)

desc <- function(x, lab) { x <- na.omit(x); tibble(Variable = lab, N = length(x), Mean = mean(x), SD = sd(x), Min = min(x),
  Q1 = unname(quantile(x, .25)), Median = median(x), Q3 = unname(quantile(x, .75)), Max = max(x)) }
D <- bind_rows(desc(A$Training_samples, "Training sample size (reported n)"), desc(A$R2, "R\u00b2"), desc(sqrt(A$R2), "r = \u221aR\u00b2"),
               desc(A$Fisher_Z, "Fisher's Z"), desc(A$vi, "Sampling variance of Z (vi)"), desc(RM$RMSE, "RMSE (Mg ha\u207b\u00b9)"))
export_table(D, "TableS07_descriptive_statistics", translate = FALSE)

P <- read_csv(file.path(TAB, "TableS05_pooled_model.csv"), show_col_types = FALSE)
mu <- P$R2[1]; pi_lo <- P$R2_PI_untransformed_lower[1]; pi_hi <- P$R2_PI_untransformed_upper[1]; pz_lo <- P$R2_PI_lower[1]; pz_hi <- P$R2_PI_upper[1]
shortref <- function(x) sub("^(.*?(19|20)[0-9]{2}[a-z]?)\\..*$", "\\1", x)
st <- A %>% group_by(Study_ID) %>% summarise(Reference = shortref(first(Reference)), Region = first(Region), mean = mean(R2), lo = min(R2), hi = max(R2), k = n(), .groups = "drop") %>%
  mutate(grp = ifelse(Region == "Central-Africa", "Central Africa", "Other regions"),
         lab = paste0(Reference, " [", Study_ID, "] (k = ", k, ")")) %>% arrange(mean) %>% mutate(lab = factor(lab, levels = unique(lab)))
st <- st %>% arrange(desc(mean)) %>% mutate(rank = row_number(), lab = paste0(Reference, "  [", Study_ID, "]  (k = ", k, ")"))
half <- ceiling(nrow(st) / 2)
part_plot <- function(d, title) {
  d <- d %>% mutate(y = rev(seq_len(n())))
  ggplot(d) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = d$y[d$rank %% 2 == 0] - .5, ymax = d$y[d$rank %% 2 == 0] + .5, fill = "#f3f5f8") +
    geom_rect(aes(xmin = pi_lo, xmax = pi_hi, ymin = -Inf, ymax = Inf, fill = "pi"), data = d[1, ], alpha = .5) +
    geom_vline(aes(xintercept = mu, linetype = "mu"), colour = INK, linewidth = .55) +
    geom_vline(aes(xintercept = c(pz_lo, pz_hi), linetype = "pz"), data = data.frame(i = 1:2), colour = GREY, linewidth = .45) +
    geom_segment(aes(x = lo, xend = hi, y = y, yend = y, colour = grp), linewidth = .8) +
    geom_point(aes(x = mean, y = y, colour = grp), size = 2.6) +
    scale_y_continuous(breaks = d$y, labels = d$lab, expand = expansion(add = .6)) +
    scale_colour_manual(values = c("Central Africa" = ORANGE, "Other regions" = BLUE), name = NULL, limits = c("Central Africa", "Other regions")) +
    scale_fill_manual(values = c(pi = "#c9daf0"), name = NULL, labels = sprintf("95%% PI, untransformed R\u00b2 (%.2f\u2013%.2f)", pi_lo, pi_hi)) +
    scale_linetype_manual(values = c(mu = "dashed", pz = "dotted"), breaks = c("mu", "pz"), name = NULL, labels = c(sprintf("Pooled R\u00b2 = %.2f", mu), sprintf("95%% PI, Fisher's Z (%.2f\u2013%.2f)", pz_lo, pz_hi))) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2), expand = expansion(mult = c(0, .01))) +
    guides(colour = guide_legend(order = 1, override.aes = list(size = 3)), fill = guide_legend(order = 2), linetype = guide_legend(order = 3)) +
    labs(x = "Reported R\u00b2 (point: study mean; line: range of models)", y = NULL, title = title,
         subtitle = sprintf("%d model observations from %d studies; ranges are descriptive, not CIs", nrow(A), nrow(st))) +
    theme_minimal(base_size = 12) +
    theme(plot.title.position = "plot", plot.title = element_text(face = "bold", size = 14, colour = INK),
          plot.subtitle = element_text(size = 10.5, colour = "#4a5059"), axis.text.y = element_text(size = 10.5, colour = INK),
          axis.text.x = element_text(size = 11, colour = INK), axis.title.x = element_text(size = 11.5, margin = margin(t = 8)),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), panel.grid.major.x = element_line(colour = GRID),
          axis.line.x = element_line(colour = "#9aa0a6"), legend.position = "bottom", legend.direction = "horizontal",
          legend.box = "vertical", legend.text = element_text(size = 10.5), legend.spacing.y = unit(0, "pt"),
          legend.margin = margin(0, 0, 0, 0), plot.margin = margin(8, 14, 8, 8))
}
sv(part_plot(st %>% filter(rank <= half), sprintf("Reported R\u00b2 by study (studies 1\u2013%d of %d)", half, nrow(st))), "Fig3a_accuracy_by_study", 7.5, 9.6)
sv(part_plot(st %>% filter(rank > half), sprintf("Reported R\u00b2 by study (studies %d\u2013%d of %d)", half + 1, nrow(st), nrow(st))), "Fig3b_accuracy_by_study", 7.5, 9.6)

cat("R/17: study summary (", nrow(S), " studies), descriptive statistics, 2 figures written\n")
