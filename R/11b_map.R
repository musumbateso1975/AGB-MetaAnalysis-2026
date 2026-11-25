suppressPackageStartupMessages({ library(sf); library(ggrepel); library(patchwork) })
shadow_lab <- function(d) geom_text_repel(data = d, aes(lx, ly, label = nom), size = 2.7, fontface = "bold", colour = "#1f2328",
  bg.color = "white", bg.r = .12, force = 0, force_pull = 10, seed = 7, max.overlaps = Inf)
carte_etudes <- function(R, FIG) {
  sf_use_s2(FALSE)
  d <- bind_rows(R$fz %>% dplyr::select(Study_ID, Reference, Journal, Country), R$rm %>% dplyr::select(Study_ID, Reference, Journal, Country)) %>%
    group_by(Study_ID) %>% mutate(Journal = dplyr::first(na.omit(Journal))) %>% ungroup() %>%
    mutate(Reference = paste(Reference, Journal, sep = " | ")) %>% distinct(Study_ID, Reference, Country)
  nmult <- n_distinct(setdiff(d$Study_ID[d$Country == "Multiple"], d$Study_ID[d$Country != "Multiple"]))
  cnt <- d %>% filter(Country != "Multiple") %>% group_by(Country) %>%
    summarise(Etudes = n_distinct(Study_ID), Art = n_distinct(Reference), .groups = "drop") %>%
    mutate(region = dplyr::recode(Country, "Democratic Republic of Congo" = "Democratic Republic of the Congo", "Cabo Verde" = "Cape Verde",
                           "United States" = "USA", "South Korea" = "South Korea"))
  w <- st_as_sf(maps::map("world", fill = TRUE, plot = FALSE)) %>% rename(region = ID) %>% filter(region != "Antarctica") %>%
    st_make_valid() %>% left_join(cnt, by = "region")
  miss <- setdiff(cnt$region, w$region); if (length(miss)) warning("Unmatched countries: ", paste(miss, collapse = ", "))
  brk <- c(0, 5, 10, 15, 20, 30); labs_cl <- c("1–5", "6–10", "11–15", "16–20", "21–30")
  PAL <- c("#ffff3d", "#f2b705", "#f08a24", "#f04e1c", "#e00000")
  w$cl <- cut(w$Etudes, brk, labels = labs_cl)
  lab <- w %>% filter(!is.na(Etudes)) %>% st_cast("POLYGON", warn = FALSE)
  lab$a <- as.numeric(st_area(st_set_crs(st_geometry(lab), NA)))
  lab <- lab %>% group_by(region) %>% slice_max(a, n = 1, with_ties = FALSE) %>% ungroup()
  xy <- st_coordinates(st_point_on_surface(st_geometry(lab))); lab <- lab %>% st_drop_geometry() %>% mutate(x = xy[, 1], y = xy[, 2],
    nom = toupper(dplyr::recode(region, "Democratic Republic of the Congo" = "DR Congo", "South Korea" = "Korea")))
  man <- tribble(~region, ~mx, ~my, "USA", -99, 39, "Canada", -105, 58, "China", 101, 34, "India", 78.5, 18.5, "Brazil", -52, -9,
                 "Australia", 134, -25, "Argentina", -66, -36, "Democratic Republic of the Congo", 22, 0, "Mexico", -103, 24,
                 "Finland", 26, 64, "Spain", -5, 39.5, "France", 2.5, 47, "Zambia", 26, -13.5, "Botswana", 23, -22.3)
  rap <- tribble(~region, ~lx, ~ly, "South Korea", 143, 46, "Japan", 156, 31, "Taiwan", 140, 18, "Vietnam", 122, 5,
                 "Myanmar", 82, 2, "Thailand", 102, -6, "Nepal", 74, 45, "Pakistan", 57, 14, "Iran", 42, 44,
                 "Cape Verde", -32, 24, "Ghana", -14, 0, "Gabon", -2, -11, "Cameroon", 2, 13, "Tanzania", 55, -10,
                 "Denmark", -4, 61, "Italy", 18, 34)
  lab <- lab %>% left_join(man, by = "region") %>% left_join(rap, by = "region") %>%
    mutate(fixe = !is.na(mx), x = ifelse(fixe, mx, x), y = ifelse(fixe, my, y), rappel = !is.na(lx),
           x = ifelse(region == "Japan", 138.5, x), y = ifelse(region == "Japan", 36.5, y))
  OCEAN <- "#b9daf1"
  oc <- tribble(~x, ~y, ~t, -42, 41, "North\nAtlantic\nOcean", -18, -28, "South\nAtlantic\nOcean", 78, -28, "Indian\nOcean",
                -150, 12, "North\nPacific\nOcean", -120, -32, "South\nPacific\nOcean")
  m <- ggplot() + geom_sf(data = w, fill = "white", colour = "#7a858f", linewidth = .08) +
    geom_sf(data = w %>% filter(!is.na(cl)), aes(fill = cl), colour = "#4b5157", linewidth = .12) +
    geom_text(data = oc, aes(x, y, label = t), fontface = "italic", size = 2.3, colour = "#4f7fa6", lineheight = .85) +
    scale_fill_manual(values = setNames(PAL, labs_cl), drop = FALSE,
                      name = "Studies per country") +
    geom_text_repel(data = lab %>% filter(fixe), aes(x, y, label = nom), size = 2.9, fontface = "bold", colour = "#1f2328",
                    bg.color = "white", bg.r = .12, force = 0, force_pull = 10, seed = 7, max.overlaps = Inf) +
    geom_segment(data = lab %>% filter(rappel), aes(x, y, xend = lx, yend = ly), colour = "#4b5157", linewidth = .25) +
    geom_point(data = lab %>% filter(rappel), aes(x, y), size = .4, colour = "#4b5157") +
    shadow_lab(lab %>% filter(rappel)) +
    geom_text_repel(data = lab %>% filter(!fixe, !rappel), aes(x, y, label = nom), size = 2.6, fontface = "bold", colour = "#1f2328",
                    bg.color = "white", bg.r = .12, seed = 7, min.segment.length = .3, segment.colour = "#5f6b75",
                    segment.size = .3, box.padding = .35, force = 2, max.overlaps = Inf, xlim = c(-165, 178), ylim = c(-52, 75)) +
    coord_sf(xlim = c(-170, 180), ylim = c(-56, 78), expand = FALSE, datum = NA) +
    theme_void(base_size = 8.5) +
    theme(panel.background = element_rect(fill = OCEAN, colour = "#1f2328", linewidth = .35),
          legend.position = c(.005, .02), legend.justification = c(0, 0),
          legend.background = element_rect(fill = "white", colour = "#5f6b75", linewidth = .3),
          legend.title = element_text(face = "bold", size = 8.5), legend.text = element_text(size = 8.5),
          legend.key.size = unit(.33, "cm"), legend.margin = margin(3, 4, 3, 4))
  top <- cnt %>% arrange(desc(Etudes), desc(Art), Country) %>% slice_head(n = 6) %>%
    mutate(rang = row_number(), col = PAL[as.integer(cut(Etudes, brk))], nom = dplyr::recode(Country, "United States" = "USA"), cx = 0:5, cy = 0)
  tp <- ggplot(top) + geom_point(aes(cx, cy, fill = col), shape = 21, size = 7, colour = "#5f6b75", stroke = .3) +
    scale_fill_identity() + geom_text(aes(cx, cy, label = rang), fontface = "bold", size = 3.6, colour = "#1f2328") +
    geom_text(aes(cx + .17, cy + .3, label = nom), hjust = 0, vjust = .5, fontface = "bold", size = 3.6, colour = "#1f2328") +
    geom_text(aes(cx + .17, cy - .2, label = sprintf("%d studies\n%d articles", Etudes, Art)), hjust = 0, vjust = .5, size = 3.1,
              colour = "#30353b", lineheight = .9) +
    scale_x_continuous(limits = c(-.2, 5.95)) + scale_y_continuous(limits = c(-.75, .6)) + theme_void() +
    labs(title = "Top 6 countries by number of studies") +
    theme(plot.title = element_text(face = "bold", size = 10, colour = "#30353b", margin = margin(b = 2)))
  p <- (m / tp + plot_layout(heights = c(3.3, 1))) & theme(plot.margin = margin(2, 3, 2, 3))
  p <- p + plot_annotation(title = "Geographic distribution of studies included in the meta-analysis",
    theme = theme(plot.title = element_text(face = "bold", size = 11.5, colour = "#1f2328"),
                  plot.background = element_rect(fill = "white", colour = "#4a4f55", linewidth = 0.8),
                  plot.margin = margin(8, 9, 8, 9)))
  ggsave(file.path(FIG, "Fig2_study_map.png"), p, width = 6.5, height = 4.05, dpi = 1000, bg = "white")
  ggsave(file.path(FIG, "Fig2_study_map.pdf"), p, width = 6.5, height = 4.05, device = grDevices::cairo_pdf, bg = "white")
  export_columns(p, "Fig2_study_map", 6.5, 4.05, dir = FIG)
}
