FIG_W   <- 6.5
FIG_DPI <- 1000
INK <- "#1f2328"; MUTED <- "#4b5157"; GRIDC <- "#d9dcdf"; REDL <- "#d62728"

theme_pub <- function(base = 10.5) {
  ggplot2::theme_bw(base_size = base) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = base + 1.5, colour = INK),
      plot.subtitle    = ggplot2::element_text(face = "bold", size = base - 0.5, colour = MUTED),
      plot.caption     = ggplot2::element_text(size = base - 1, colour = MUTED, hjust = 0),
      axis.title       = ggplot2::element_text(face = "bold", size = base, colour = INK),
      axis.text        = ggplot2::element_text(size = base - 1, colour = INK),
      axis.text.x      = ggplot2::element_text(face = "bold", lineheight = .9),
      panel.border     = ggplot2::element_rect(colour = INK, linewidth = .6),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = GRIDC, linetype = "dashed", linewidth = .3),
      panel.grid.minor = ggplot2::element_blank(),
      legend.title     = ggplot2::element_text(face = "bold", size = base - 1),
      legend.text      = ggplot2::element_text(size = base - 1.5),
      plot.title.position = "plot"
    )
}

COL_IN <- c(one = 90, two = 190) / 25.4
COL_DPI <- 1000

write_scaled <- function(p, base, w, h, target_w) {
  s <- target_w / w
  ragg::agg_png(paste0(base, ".png"), width = round(target_w * COL_DPI), height = round(h * s * COL_DPI),
                units = "px", res = COL_DPI * s, background = "white")
  print(p); grDevices::dev.off()
  magick::image_write(magick::image_read(paste0(base, ".png")), paste0(base, ".png"), density = COL_DPI)
  svglite::svglite(paste0(base, ".svg"), width = w, height = h, bg = "white"); print(p); grDevices::dev.off()
  grDevices::cairo_pdf(paste0(base, ".pdf"), width = w, height = h, bg = "white"); print(p); grDevices::dev.off()
  sv <- readLines(paste0(base, ".svg"), warn = FALSE)
  i <- grep("<svg ", sv)[1]
  sv[i] <- sub("width='[^']*'", sprintf("width='%.1fmm'", target_w * 25.4), sv[i])
  sv[i] <- sub("height='[^']*'", sprintf("height='%.1fmm'", h * s * 25.4), sv[i])
  writeLines(sv, paste0(base, ".svg"))
}

export_columns <- function(p, name, w, h, h1 = NULL, dir = "output/figures", w1 = NULL) {
  if (!grepl("^Fig[0-9]", name) || !isTRUE(getOption("agb.article_figures", FALSE))) return(invisible())
  dir.create(file.path(dir, "two_column"), showWarnings = FALSE, recursive = TRUE)
  tw <- min(COL_IN[["two"]], (240 / 25.4) * w / h)
  write_scaled(p, file.path(dir, "two_column", name), w, h, tw)
  if (!is.null(w1)) {
    dir.create(file.path(dir, "one_column"), showWarnings = FALSE, recursive = TRUE)
    write_scaled(p, file.path(dir, "one_column", name), w1, if (is.null(h1)) h else h1, COL_IN[["one"]])
  }
}

save_pub <- function(p, name, h, w = FIG_W, dir = "output/figures", h1 = NULL, w1 = NULL) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(file.path(dir, paste0(name, ".png")), p, width = w, height = h, dpi = FIG_DPI, bg = "white")
  ggplot2::ggsave(file.path(dir, paste0(name, ".pdf")), p, width = w, height = h, device = grDevices::cairo_pdf, bg = "white")
  export_columns(p, name, w, h, h1, dir, w1)
  cat("✓", name, sprintf("(%.1f × %.1f in, %d dpi)\n", w, h, FIG_DPI))
}

fmt_p <- function(p) ifelse(is.na(p), "n.e.", ifelse(p < 0.001, "p < 0.001", sprintf("p = %.3f", p)))

ml_test_z <- function(d, g) {
  d <- d[!is.na(d[[g]]), ]; d[[g]] <- droplevels(factor(d[[g]]))
  if (nlevels(d[[g]]) < 2) return(NA_character_)
  out <- tryCatch({
    m <- metafor::rma.mv(Fisher_Z, vi_imp, mods = stats::as.formula(paste("~", g)),
                         random = ~ 1 | Cluster / Study_ID / Obs_ID, data = d, method = "REML", sparse = TRUE)
    w <- clubSandwich::Wald_test(m, constraints = clubSandwich::constrain_zero(2:m$p),
                                 vcov = "CR2", cluster = d$Cluster, test = "HTZ")
    if (is.na(w$p_val) || !is.finite(w$Fstat) || w$Fstat > 1e4 || w$df_denom < 2) "Multilevel robust test (CR2): not estimable (too few studies per level)"
    else sprintf("Multilevel robust test (CR2): F = %.2f, %s", w$Fstat, fmt_p(w$p_val))
  }, error = function(e) "Multilevel robust test: not estimable")
  out
}

ml_test_rmse <- function(d, g) {
  d <- d[!is.na(d[[g]]) & is.finite(d$log_RMSE), ]; d[[g]] <- droplevels(factor(d[[g]]))
  if (nlevels(d[[g]]) < 2) return(NA_character_)
  tryCatch({
    m <- nlme::lme(stats::as.formula(paste("log_RMSE ~", g)), random = ~ 1 | Study_ID, data = d, method = "REML")
    p <- length(nlme::fixef(m))
    w <- clubSandwich::Wald_test(m, constraints = clubSandwich::constrain_zero(2:p), vcov = "CR2", test = "HTZ")
    if (is.na(w$p_val) || !is.finite(w$Fstat) || w$Fstat > 1e4 || w$df_denom < 2) "Multilevel robust test on log RMSE (CR2): not estimable (too few studies per level)"
    else sprintf("Multilevel robust test on log RMSE (CR2): F = %.2f, %s", w$Fstat, fmt_p(w$p_val))
  }, error = function(e) "Multilevel robust test: not estimable")
}

violin_pub <- function(d, x, y, cols, title, subtitle, ylab, ann, log_y = FALSE, hline = NULL, hlab = NULL,
                       ann_size = 2.75) {
  p <- ggplot2::ggplot(d, ggplot2::aes(.data[[x]], .data[[y]], fill = .data[[x]])) +
    ggplot2::geom_violin(alpha = .38, trim = FALSE, scale = "width", colour = INK, linewidth = .35, width = .85) +
    ggplot2::geom_jitter(ggplot2::aes(colour = .data[[x]]), width = .09, size = .7, alpha = .55, show.legend = FALSE) +
    ggplot2::geom_boxplot(width = .17, alpha = .9, colour = "black", linewidth = .45, outlier.shape = NA) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::scale_colour_manual(values = cols, guide = "none") +
    ggplot2::scale_x_discrete(labels = function(l) { l <- gsub("\\+", "\n+", gsub("-", "-\n", l)); vapply(l, function(z) if (grepl("\n", z)) z else stringr::str_wrap(z, 9), "") }) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = ylab) + theme_pub()
  if (log_y) p <- p + ggplot2::scale_y_log10(labels = scales::label_number(accuracy = 1, big.mark = " "),
                                             expand = ggplot2::expansion(mult = c(.03, .30)))
  else p <- p + ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.03, .34)))
  if (!is.null(hline)) p <- p +
    ggplot2::geom_hline(yintercept = hline, colour = REDL, linetype = "dashed", linewidth = .7)
  if (!is.null(hlab)) p <- p +
    ggplot2::annotate("label", x = Inf, y = hline, label = hlab, hjust = 1.02, vjust = -.25, size = 2.7,
                      colour = REDL, fontface = "bold", label.size = 0, fill = "white", alpha = .85)
  p + ggplot2::geom_text(data = ann, ggplot2::aes(x = .data[[x]], y = ypos, label = lab), inherit.aes = FALSE,
                         vjust = 0, size = ann_size, fontface = "bold", colour = INK, lineheight = .92)
}
