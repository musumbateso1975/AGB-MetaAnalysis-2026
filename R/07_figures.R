if (!exists("dat_fz")) { source("R/00_setup.R"); source("R/01_data_import.R") }
source("R/00b_figure_style.R")
if (!file.exists(file.path(CACHE_DIR, "meta_results.rds"))) source("R/03_fisher_z_meta.R")
meta <- readRDS(file.path(CACHE_DIR, "meta_results.rds"))
POOLED_Z <- meta$pooled_Z
HLAB <- sprintf("Pooled Z = %.2f", POOLED_Z)
FIGD <- "output/figures"
set.seed(2024)
cat("=== 07_figures.R — Descriptive figures (v3) ===\n\n")

ann_z <- function(d, g, ytop) d |> dplyr::filter(!is.na(.data[[g]])) |>
  dplyr::group_by(.data[[g]]) |>
  dplyr::summarise(lab = sprintf("n = %d\nk = %d\nZ = %.2f\nR² = %.2f", dplyr::n(), dplyr::n_distinct(Study_ID),
                                 mean(Fisher_Z), tanh(mean(Fisher_Z))^2), .groups = "drop") |>
  dplyr::mutate(ypos = ytop)
ann_r <- function(d, g, ytop) d |> dplyr::filter(!is.na(.data[[g]])) |>
  dplyr::group_by(.data[[g]]) |>
  dplyr::summarise(lab = sprintf("n = %d\nk = %d\nMd = %.1f", dplyr::n(), dplyr::n_distinct(Study_ID),
                                 stats::median(RMSE)), .groups = "drop") |>
  dplyr::mutate(ypos = ytop)

two_panel <- function(g, cols, label, name, min_n = 1, stack = TRUE, h = NULL, order = TRUE) {
  dz <- dat_fz |> dplyr::filter(!is.na(.data[[g]]), .data[[g]] != "")
  dr <- dat_rmse |> dplyr::filter(!is.na(.data[[g]]), .data[[g]] != "")
  keep_z <- names(which(table(as.character(dz[[g]])) >= min_n)); keep_r <- names(which(table(as.character(dr[[g]])) >= min_n))
  dz <- dz |> dplyr::filter(as.character(.data[[g]]) %in% keep_z)
  dr <- dr |> dplyr::filter(as.character(.data[[g]]) %in% keep_r)
  if (order) {
    oz <- dz |> dplyr::group_by(.data[[g]]) |> dplyr::summarise(m = stats::median(Fisher_Z)) |> dplyr::arrange(dplyr::desc(m))
    dz[[g]] <- factor(as.character(dz[[g]]), levels = as.character(oz[[g]]))
    orr <- dr |> dplyr::group_by(.data[[g]]) |> dplyr::summarise(m = stats::median(RMSE)) |> dplyr::arrange(m)
    dr[[g]] <- factor(as.character(dr[[g]]), levels = as.character(orr[[g]]))
  } else { dz[[g]] <- droplevels(factor(dz[[g]])); dr[[g]] <- droplevels(factor(dr[[g]])) }
  pa <- violin_pub(dz, g, "Fisher_Z", cols, sprintf("(a) Fisher's Z by %s", label), ml_test_z(dz, g), "Fisher's Z",
                   ann_z(dz, g, max(dz$Fisher_Z) + .5), hline = POOLED_Z, hlab = NULL)
  pb <- violin_pub(dr, g, "RMSE", cols, sprintf("(b) RMSE by %s", label), ml_test_rmse(dr, g),
                   expression(bold("RMSE (Mg ha"^-1*", log scale)")), ann_r(dr, g, max(dr$RMSE) * 2.2), log_y = TRUE)
  if (!stack) { pa$labels$subtitle <- sub(": ", ":\n", pa$labels$subtitle); pb$labels$subtitle <- sub(": ", ":\n", sub(" on log RMSE", "\non log RMSE", pb$labels$subtitle)) }
  if (max(nlevels(droplevels(factor(dz[[g]]))), nlevels(droplevels(factor(dr[[g]])))) > 7) {
    sm <- ggplot2::theme(axis.text.x = ggplot2::element_text(size = 9, face = "bold", lineheight = .85))
    pa <- pa + sm; pb <- pb + sm }
  p <- if (stack) pa / pb else (pa | pb)
  save_pub(p, name, h = if (is.null(h)) (if (stack) 8.6 else 4.6) else h, dir = FIGD)
  invisible(p)
}

dat_fz$GEDI_lab   <- factor(ifelse(dat_fz$GEDI_used == "Yes", "GEDI used", "No GEDI"), levels = c("GEDI used", "No GEDI"))
dat_rmse$GEDI_lab <- factor(ifelse(dat_rmse$GEDI_used == "Yes", "GEDI used", "No GEDI"), levels = c("GEDI used", "No GEDI"))
two_panel("GEDI_lab", c("GEDI used" = unname(GEDI_COLORS["Yes"]), "No GEDI" = unname(GEDI_COLORS["No"])),
          "GEDI Use", "FigS01_GEDI_use", stack = FALSE, order = FALSE, h = 5.4)

cat("--- algorithms ---\n")
MCOL <- c("Statistical" = "#1f77b4", "Machine-learning" = "#2ca02c", "Deep-learning" = "#d62728")
algo_plot <- function(d, val, xlab, title, name, ref = NULL, log_x = FALSE) {
  a <- d |> dplyr::filter(!is.na(Model_algorithm), Model_algorithm != "") |> dplyr::group_by(Model_algorithm) |>
    dplyr::summarise(n = dplyr::n(), k = dplyr::n_distinct(Study_ID), m = mean(.data[[val]]),
                     se = stats::sd(.data[[val]]) / sqrt(dplyr::n()), cat = dplyr::first(as.character(Model_category)), .groups = "drop") |>
    dplyr::filter(n >= 5, k >= 2) |>
    dplyr::mutate(lo = m - stats::qt(.975, n - 1) * se, hi = m + stats::qt(.975, n - 1) * se, lab = sprintf("n = %d · k = %d", n, k))
  xr <- range(c(if (val == "RMSE") pmax(a$lo, 0) else a$lo, a$hi)); pad <- diff(xr) * .42
  p <- ggplot2::ggplot(a, ggplot2::aes(m, stats::reorder(Model_algorithm, m), colour = cat)) +
    { if (!is.null(ref)) ggplot2::geom_vline(xintercept = ref, linetype = "dashed", colour = INK, linewidth = .6) } +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = if (val == "RMSE") pmax(lo, 0) else lo, xmax = hi), height = .3, linewidth = .7) +
    ggplot2::geom_point(ggplot2::aes(size = n)) +
    ggplot2::geom_text(ggplot2::aes(x = xr[2] + pad * .08, label = lab), hjust = 0, size = 2.9, colour = INK, fontface = "bold") +
    ggplot2::scale_colour_manual(values = MCOL, name = NULL) + ggplot2::scale_size_area(max_size = 6, name = "Obs.") +
    ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(6)) +
    ggplot2::guides(colour = ggplot2::guide_legend(order = 1, nrow = 3), size = ggplot2::guide_legend(order = 2)) +
    ggplot2::coord_cartesian(xlim = c(xr[1], xr[2] + pad), clip = "off") +
    ggplot2::labs(title = title, subtitle = "Algorithms with ≥ 5 observations from ≥ 2 studies\nDescriptive 95% CI (observations within a study are not independent)",
                  x = xlab, y = NULL) +
    theme_pub() + ggplot2::theme(panel.grid.major.y = ggplot2::element_line(colour = GRIDC, linewidth = .3, linetype = "solid"),
                                 legend.position = "bottom", legend.box = "horizontal", axis.text.y = ggplot2::element_text(size = 9.5))
  save_pub(p, name, h = 2.6 + .32 * nrow(a), dir = FIGD)
}
algo_plot(dat_fz, "Fisher_Z", "Mean Fisher's Z (95% CI)", "Mean Fisher's Z by algorithm", "FigS02_R2_by_algorithm", ref = POOLED_Z)
algo_plot(dat_rmse, "RMSE", expression(bold("Mean RMSE (Mg ha"^-1*", 95% CI)")), "Mean RMSE by algorithm", "FigS03_RMSE_by_algorithm")

cat("\n=== 07_figures.R complete — output/figures/descriptive/ ===\n")
