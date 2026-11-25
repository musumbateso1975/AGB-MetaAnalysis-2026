suppressPackageStartupMessages({ library(dplyr); library(readr) })
source("R/00d_prepare_analysis.R", encoding = "UTF-8")

rmA <- prep(rm0) %>% filter(RMSE_unit == "Mg_ha") %>% mutate(RMSE = as.numeric(RMSE))

CANON <- c("Random Forest" = "Breiman (2001)", "Quantile Random Forest" = "Breiman (2001)",
           "SVM" = "Cortes and Vapnik (1995)", "SVR" = "Drucker et al. (1997)",
           "XGBoost" = "Chen and Guestrin (2016)",
           "Artificial Neural Network (ANN)" = "Rumelhart et al. (1986)", "Multilayer Perceptron" = "Rumelhart et al. (1986)",
           "Backpropagation Neural Network" = "Rumelhart et al. (1986)",
           "CNN" = "LeCun et al. (1998)", "CNN LSTM" = "LeCun et al. (1998); Hochreiter and Schmidhuber (1997)",
           "GWR" = "Brunsdon et al. (1996)")

short_ref <- function(x) sub("^(.*?\\b(19|20)\\d{2}[a-c]?)\\b.*$", "\\1", x)
refs <- fz %>% distinct(Study_ID, Reference) %>% group_by(Study_ID) %>% slice(1) %>% ungroup() %>%
  mutate(Ref = short_ref(Reference))

alg <- A %>% filter(!is.na(Model_algorithm), Model_algorithm != "") %>%
  group_by(Model_category, Model_algorithm) %>%
  summarise(n_obs = n(), k_studies = n_distinct(Study_ID),
            R2_min = min(R2), R2_max = max(R2), studies = list({u <- unique(Study_ID); u[order(as.numeric(gsub("\\D", "", u)), u)]}), .groups = "drop")
rm_alg <- rmA %>% filter(!is.na(Model_algorithm), Model_algorithm != "") %>% group_by(Model_algorithm) %>%
  summarise(n_obs_RMSE = n(), RMSE_min = min(RMSE), RMSE_max = max(RMSE), .groups = "drop")

TS8 <- alg %>% left_join(rm_alg, by = "Model_algorithm") %>%
  mutate(Method_reference = unname(ifelse(Model_algorithm %in% names(CANON), CANON[Model_algorithm], "")),
         Study_IDs = vapply(studies, function(s) paste(s, collapse = "; "), ""),
         Included_studies = vapply(studies, function(s) paste(refs$Ref[match(s, refs$Study_ID)], collapse = "; "), "")) %>%
  select(-studies) %>%
  mutate(Model_category = factor(Model_category, c("Statistical", "Machine-learning", "Deep-learning"))) %>%
  arrange(Model_category, desc(n_obs), Model_algorithm) %>%
  mutate(across(c(R2_min, R2_max), ~ round(.x, 2)), across(c(RMSE_min, RMSE_max), ~ round(.x, 1)))

export_table(TS8, "TableS09_algorithms_references", translate = FALSE)
cat(sprintf("13_algorithm_table.R: %d algorithms, %d observations, %d studies\n", nrow(TS8), sum(TS8$n_obs), n_distinct(A$Study_ID)))
