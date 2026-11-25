suppressPackageStartupMessages({ library(dplyr); library(readr) })
if (!file.exists("AGB_MetaAnalysis.Rproj")) stop("Run from the project folder.")
source("R/00c_export_helpers.R")
res <- list()
chk <- function(group, test, ok, detail = "") res[[length(res) + 1]] <<- data.frame(group = group, test = test, result = ifelse(isTRUE(ok), "PASS", "FAIL"), detail = detail)

fz <- read_csv("data/DATA_FisherZ.csv", col_types = cols(.default = "c"), na = character(0))
rm <- read_csv("data/DATA_RMS.csv", col_types = cols(.default = "c"), na = character(0))
ms <- read_csv("data/AGB_master_database.csv", col_types = cols(.default = "c"), na = character(0))
num <- function(x) suppressWarnings(as.numeric(x))

chk("Identifiers", "Obs_ID unique in master database", !anyDuplicated(ms$Obs_ID), sprintf("%d rows", nrow(ms)))
chk("Identifiers", "Obs_ID unique in R2 database", !anyDuplicated(fz$Obs_ID), sprintf("%d rows", nrow(fz)))
chk("Identifiers", "Obs_ID unique in RMSE database", !anyDuplicated(rm$Obs_ID), sprintf("%d rows", nrow(rm)))
chk("Identifiers", "Obs_No unique in master database", !anyDuplicated(ms$Obs_No))
chk("Identifiers", "Obs_ID format M<number>", all(grepl("^M[0-9]+$", c(ms$Obs_ID, fz$Obs_ID, rm$Obs_ID))))
chk("Identifiers", "Study_ID format S1-S72", all(grepl("^S([1-9]|[1-6][0-9]|7[0-2])$", c(ms$Study_ID, fz$Study_ID, rm$Study_ID))))
chk("Identifiers", "72 studies in master database", n_distinct(ms$Study_ID) == 72, sprintf("%d", n_distinct(ms$Study_ID)))
chk("Identifiers", "R2 observations present in master database", all(fz$Obs_ID %in% ms$Obs_ID))
chk("Identifiers", "RMSE observations present in master database", all(rm$Obs_ID %in% ms$Obs_ID))
m1 <- inner_join(select(fz, Obs_ID, S1 = Study_ID, N1 = Obs_No), select(ms, Obs_ID, S2 = Study_ID, N2 = Obs_No), by = "Obs_ID")
chk("Identifiers", "Study_ID and Obs_No identical between R2 and master databases", all(m1$S1 == m1$S2 & m1$N1 == m1$N2))
m2 <- inner_join(select(rm, Obs_ID, S1 = Study_ID, N1 = Obs_No), select(ms, Obs_ID, S2 = Study_ID, N2 = Obs_No), by = "Obs_ID")
chk("Identifiers", "Study_ID and Obs_No identical between RMSE and master databases", all(m2$S1 == m2$S2 & m2$N1 == m2$N2))
chk("Identifiers", "One reference per study", all(tapply(sub("^(.*? [0-9]{4}[a-z]?)\\. .*$", "\\1", ms$Reference), ms$Study_ID, function(x) length(unique(x))) == 1))
dupv <- fz %>% count(Study_ID, Site_Name, Sensor_combination, Model_algorithm, Training_samples, R2) %>% filter(n > 1)
dupo <- fz %>% semi_join(dupv, by = c("Study_ID", "Site_Name", "Sensor_combination", "Model_algorithm", "Training_samples", "R2")) %>%
  group_by(Study_ID, R2) %>% summarise(ids = paste(Obs_ID, collapse = "/"), .groups = "drop")
res[[length(res) + 1]] <- data.frame(group = "Duplicates", test = "Identical R2 for the same study, site, sensors, algorithm and n (distinct models with equal rounded R2, or duplicates)",
  result = ifelse(nrow(dupv) == 0, "PASS", "REVIEW"), detail = paste(sprintf("%s R2 %s: %s", dupo$Study_ID, dupo$R2, dupo$ids), collapse = "; "))

r2 <- num(fz$R2); z <- num(fz$Fisher_Z)
chk("Ranges", "0 < R2 <= 1", all(r2 > 0 & r2 <= 1, na.rm = TRUE) && !anyNA(r2))
chk("Ranges", "Fisher_Z = atanh(sqrt(R2)) (tolerance 1e-6)", max(abs(z - atanh(sqrt(r2)))) < 1e-6, sprintf("max |diff| = %.2e", max(abs(z - atanh(sqrt(r2))))))
n <- num(fz$Training_samples)
chk("Ranges", "Training sample size > 3 where reported", all(n[!is.na(n)] > 3), sprintf("%d missing (median imputed in the models)", sum(is.na(n))))
chk("Missing values", "22 R2 values without training n", sum(is.na(n)) == 22, sprintf("%d", sum(is.na(n))))
rmv <- num(rm$RMSE)
chk("Units", "RMSE units are Mg_ha or listed non-comparable units", all(rm$RMSE_unit %in% c("Mg_ha")), paste(unique(rm$RMSE_unit), collapse = ", "))
chk("Ranges", "RMSE > 0", all(rmv > 0, na.rm = TRUE) && !anyNA(rmv))
yr <- num(ms$Year); chk("Ranges", "Publication year 2008-2026", all(yr >= 2008 & yr <= 2026, na.rm = TRUE))

vf <- val_design(fz$Validation_design_text); vr <- val_design_rmse(rm$Validation_design_text)
chk("Validation design", "10 training-fit R2 values (S2, S6, S17, S31)", sum(vf == "Training fit") == 10 && setequal(unique(fz$Study_ID[vf == "Training fit"]), c("S2", "S6", "S17", "S31")))
chk("Validation design", "2 training-fit RMSE values (S2)", sum(vr == "Training fit") == 2 && all(rm$Study_ID[vr == "Training fit"] == "S2"))
pf <- primary_r2(fz); pr <- primary_rmse(rm)
chk("Exclusion rules", "7 nested subsets (rule b) and 1 GEDI-referenced R2 (rule e) excluded", sum(is_nested(fz$Validation_design_text)) == 7 && sum(is_gedi_ref(fz$Reference_type)) == 1)
chk("Validation design", "Primary R2 set: 258 values, 62 studies", sum(pf) == 258 && n_distinct(fz$Study_ID[pf]) == 62)
chk("Validation design", "Primary RMSE set: 294 values, 66 studies", sum(pr) == 294 && n_distinct(rm$Study_ID[pr]) == 66)
chk("Validation design", "Primary R2 set by design: 87 independent, 71 cross-validation, 100 not reported",
    identical(as.integer(table(factor(vf[pf], c("Independent validation", "Cross-validation", "Not reported")))), c(87L, 71L, 100L)))

rem <- read_csv("output/tables/TableS04_removed_values.csv", show_col_types = FALSE)
chk("Counts", "328 R2 values extracted = 258 analysed + 70 not analysed", 258 + sum(rem$n_obs[rem$Dataset == "Fisher_Z"]) == 328)
chk("Counts", "343 RMSE values extracted = 294 analysed + 49 not analysed", 294 + sum(rem$n_obs[rem$Dataset == "RMSE"]) == 343)

S5 <- read_csv("output/tables/TableS05_pooled_model.csv", show_col_types = FALSE)
chk("Data-model correspondence", "Pooled model fitted on 258 values, 62 studies, 56 clusters", S5$n_obs == 258 && S5$k_studies == 62 && S5$k_clusters == 56)
S30 <- read_csv("output/tables/TableS30_sensitivity.csv", show_col_types = FALSE)
S31 <- read_csv("output/tables/TableS31_sensitivity_additional.csv", show_col_types = FALSE)
chk("Sensitivity concordance", "Primary R2 identical in Tables S5, S30 and S31 (tolerance 1e-4)",
    abs(S5$R2 - S30$R2[1]) < 1e-4 && abs(S5$R2 - S31$R2[1]) < 1e-4)
chk("Sensitivity concordance", "Including training-fit values: 268 values, 65 studies, 59 clusters",
    with(S31[S31$Scenario == "Including training-fit R²", ], n_obs == 268 && k_studies == 65 && k_clusters == 59))
chk("Sensitivity concordance", "Training-fit scenario identical in Tables S30 and S31 (tolerance 1e-4)",
    abs(S30$R2[S30$Scenario == "Including training-fit R²"] - S31$R2[S31$Scenario == "Including training-fit R²"]) < 1e-4)

expected <- c(sprintf("output/tables/%s.csv", OUT_TAB), "output/tables/Table1_datasets.csv", "output/tables/Table2_subgroups_CI_PI.csv",
              "output/tables/Table3_models_M1_M4.csv", "output/tables/Table4_within_study_contrasts.csv",
              sprintf("output/figures/Fig%s.png", c("1_PRISMA_flow_diagram", "2_study_map", "3a_accuracy_by_study", "3b_accuracy_by_study",
                      "4_sensor_configuration_CI_PI", "5_models_M1_M4", "6_within_study_contrasts", "7_sensitivity_analyses")),
              sprintf("output/figures/FigS0%d_%s.png", 1:7, c("GEDI_use", "R2_by_algorithm", "RMSE_by_algorithm", "LiDAR_matched_pairs",
                      "terrain_sample_size_reference", "validation_resolution_region", "funnel_plot")))
miss <- expected[!file.exists(expected)]
chk("Outputs", "Every table and figure cited in the article and S2 was produced", length(miss) == 0, paste(miss, collapse = "; "))

man <- read_csv("tests/results_manifest.csv", col_types = cols(.default = "c"))
cmp <- lapply(seq_len(nrow(man)), function(i) {
  m <- man[i, ]; f <- file.path("output/tables", m$output_file)
  d <- read_csv(f, col_types = cols(.default = "c"), na = character(0))
  if (!is.na(m$key_column) && nzchar(m$key_column)) {
    kc <- strsplit(m$key_column, "\\|")[[1]]; kv <- strsplit(m$key_value, "\\|")[[1]]
    sel <- Reduce(`&`, lapply(seq_along(kc), function(j) d[[kc[j]]] == kv[j]))
    d <- d[sel, , drop = FALSE]
  }
  raw <- if (nrow(d) == 1) d[[m$value_column]] else NA_character_
  nums <- regmatches(raw, gregexpr("-?[0-9]+\\.?[0-9]*", raw))[[1]]
  out <- if (length(nums) >= as.integer(m$number_index)) as.numeric(nums[as.integer(m$number_index)]) else NA_real_
  rep <- as.numeric(m$reported_value)
  dec <- ifelse(grepl("\\.", m$reported_value), nchar(sub("^-?[0-9]+\\.", "", m$reported_value)), 0)
  tol <- 0.5 * 10^(-dec) + 1e-9
  data.frame(id = m$id, manuscript_location = m$manuscript_location, result = m$result, reported = m$reported_value,
             output_value = out, tolerance = tol, status = ifelse(!is.na(out) && abs(out - rep) <= tol, "PASS", "FAIL"),
             s2_item = m$s2_item, script = m$script, output_file = m$output_file)
})
cmp <- bind_rows(cmp)
chk("Manuscript vs outputs", sprintf("%d reported results match the outputs within half a unit of the last reported digit", nrow(cmp)),
    all(cmp$status == "PASS"), paste(cmp$id[cmp$status != "PASS"], collapse = ", "))

out <- bind_rows(res)
cat(sprintf("Consistency checks: %d PASS, %d REVIEW, %d FAIL\n", sum(out$result == "PASS"), sum(out$result == "REVIEW"), sum(out$result == "FAIL")))
cat(sprintf("Manuscript vs outputs: %d PASS, %d FAIL\n", sum(cmp$status == "PASS"), sum(cmp$status == "FAIL")))
print(out[out$result != "PASS", ])
print(cmp[cmp$status == "FAIL", ])
