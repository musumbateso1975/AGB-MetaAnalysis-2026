# Beyond pooled accuracy: A dependence-aware meta-analysis of LiDAR, GEDI and multi-sensor forest biomass estimation for carbon monitoring

Data and R scripts of the article by Musumba Teso, Okole, Birindwa, Muhindo Sahani, Mwanjalolo and Karume (Remote Sensing Applications: Society and Environment). Running `run_all.R` regenerates every table and figure of the article and of Supplementary Material S2.

## How to run

1. Open `AGB_MetaAnalysis.Rproj` in RStudio, or set the working directory to this folder.
2. Run `source("run_all.R", encoding = "UTF-8")` (R ≥ 4.2; about 15–20 min). Missing CRAN packages are installed by `R/00_setup.R`.
3. Results are written to `output/figures/` (PNG 1000 dpi and vector PDF) and `output/tables/` (CSV, and `All_tables.xlsx` with one sheet per table). Publication-size versions of the article figures (PNG 1000 dpi, SVG and PDF) are written to `output/figures/one_column/` and `output/figures/two_column/` only after `options(agb.article_figures = TRUE)`. Random steps use fixed seeds (`set.seed()` in the scripts).
4. Run `Rscript tests/consistency_checks.R` from this folder after `run_all.R`. It checks identifiers, counts, ranges, units and the data–model correspondence, compares the results reported in the article (`tests/results_manifest.csv`) with the outputs, and prints the results to the console.
5. `renv.lock` records the package versions used for the submitted results (R 4.3.3). To restore them: `install.packages("renv"); renv::restore()`.

## Data (`data/`)

| File | Content |
|---|---|
| `AGB_master_database.csv` | All observations of the 72 included studies (318 rows: one row per reported model); version v16 of 23 November 2026. |
| `DATA_FisherZ.csv` | R² values, converted to Fisher's Z (276 rows, 66 studies). |
| `DATA_RMS.csv` | RMSE values in Mg ha⁻¹ (304 rows, 68 studies). |

Validation design is coded for each value from `Validation_design_text` (`val_design()` and `val_design_rmse()` in `R/00c_export_helpers.R`). Three kinds of values stay in the files but are excluded from the primary analysis by `primary_r2()` and `primary_rmse()` (`R/00c_export_helpers.R`): training-fit values (10 R² values from S2, S6, S17 and S31; 2 RMSE values from S2; rule d), subsets nested in an aggregate retained from the same study (7 R² and 7 RMSE values from S40, S50 and S60, marked « sous-ensemble de l'agrégat » in `Validation_design_text`; rule b) and values computed against GEDI (1 R² and 1 RMSE value from S4, `Reference_type` containing « GEDI »; rule e). The primary analysis uses 258 R² values (62 studies, 56 clusters) and 294 RMSE values (66 studies). Each excluded group enters a sensitivity analysis (Table S31). The v16 coding corrections are summarised in Supplementary Material S2 (Methods S6) and listed, with their sources, in `Journal_corrections_v16.csv`, delivered with the v16 databases.

`Obs_ID` and `Obs_No` identify each observation uniquely and are the same in the three databases, so an R² and an RMSE reported for the same model share the same identifier. `Study_ID` (S1–S72) identifies the study. Variables are defined in Table S2 of the supplement (`output/tables/TableS02_data_dictionary.csv`).

## Scripts (`R/`)

`run_all.R` runs the scripts in the order below; `R/08_tables.R` runs last because it gathers all tables in the workbook. `00b`–`00d` are helpers sourced by the other scripts.

| Script | Content |
|---|---|
| `00_setup.R` | Packages, global options, output folders and helper functions. |
| `00b_figure_style.R` | Shared figure theme, PNG (1000 dpi) and PDF export and optional publication-size exports of the article figures (PNG 1000 dpi and SVG). |
| `00c_export_helpers.R` | Export of tables to output/tables (English labels, rounded values, file names numbered as in Supplementary Material S2). |
| `00d_prepare_analysis.R` | Analysed dataset, clusters, sampling variances and model helpers shared by R/10 and R/12-R/16. |
| `01_data_import.R` | Reads the databases, harmonises variables, defines clusters and sampling variances. |
| `02_descriptive_stats.R` | Descriptive statistics of the analysed observations. |
| `03_fisher_z_meta.R` | Two-level random-effects model of Fisher's Z with robust variance; leave-one-study-out. |
| `10_multilevel_reanalysis.R` | Primary analysis: three-level model with CR2 inference, subgroups, contrasts, meta-regression, sensitivity. |
| `07_figures.R` | GEDI use and accuracy by algorithm (Figs S1–S3). |
| `11_main_figures.R` | PRISMA diagram, map, sensor configurations, within-study contrasts, sensitivity and funnel plot. |
| `11b_map.R` | Map of the included studies by country (sourced by 11_main_figures.R). |
| `12_terrain_samplesize_gedi.R` | Accuracy by terrain, training sample size and reference data. |
| `13_algorithm_table.R` | Algorithms, method references and the included studies using them. |
| `14_estimands_paired_conventional.R` | Estimand classes, matched LiDAR pairs, conventional vs dependence-aware analyses, Central Africa, precision of contrasts. |
| `15_final_validation.R` | Effect-size strategies, matched LiDAR pairs by study, models M1-M4 and mechanisms. |
| `16_two_scales_and_pairs.R` | Untransformed vs Fisher's Z scale, M1-M4 on both scales, matched LiDAR pairs (sign test) and GEDI roles. |
| `17_study_summary_figures.R` | Accuracy by study (Fig. 3). |
| `18_robustness_analyses.R` | Robustness analyses: validation design, spatial resolution, region and sample size as moderators; cluster wild bootstrap of moderator tests; country-level clustering; within-study contrasts under assumed sampling correlations; added sensitivity scenarios; RMSE scale; associations between moderators. |
| `19_model_diagnostics.R` | Diagnostics of the primary three-level model on both scales: standardized residuals, cluster-level outlier tests and Cook's distances, normality of residuals and cluster random effects, profile-likelihood confidence intervals of the variance components. |
| `08_tables.R` | Dataset table, moderator categories and tests, data dictionary, included studies, removed observations, workbook. |

## Outputs

File names follow the numbering of the article (Fig1–Fig7, with Fig. 3 in two parts; Table1–Table4) and of Supplementary Material S2 (FigS01–FigS08; TableS02–TableS35). Tables S1 and S21 of the supplement are documentary tables and are not produced by the scripts.

| Figure | Content |
|---|---|
| `Fig1_PRISMA_flow_diagram.png` | PRISMA 2020 flow diagram |
| `Fig2_study_map.png` | Included studies by country |
| `Fig3a_accuracy_by_study.png` | Reported R² by study, studies 1–33 |
| `Fig3b_accuracy_by_study.png` | Reported R² by study, studies 34–66 |
| `Fig4_sensor_configuration_CI_PI.png` | Mean R² by sensor configuration with CI and PI |
| `Fig5_models_M1_M4.png` | Models M1–M4 and sensor-combination test |
| `Fig6_within_study_contrasts.png` | Within-study contrasts |
| `Fig7_sensitivity_analyses.png` | Sensitivity of the pooled R² |
| `FigS01_GEDI_use.png` | Fig. S1: Fisher's Z and RMSE with and without GEDI |
| `FigS02_R2_by_algorithm.png` | Fig. S2: mean Fisher's Z by algorithm |
| `FigS03_RMSE_by_algorithm.png` | Fig. S3: mean RMSE by algorithm |
| `FigS04_LiDAR_matched_pairs.png` | Fig. S4: matched LiDAR comparisons |
| `FigS05_terrain_sample_size_reference.png` | Fig. S5: pooled R² by terrain, training sample size and reference data |
| `FigS06_validation_resolution_region.png` | Fig. S6: pooled R² by validation design, spatial resolution and region |
| `FigS07_funnel_plot.png` | Fig. S7: funnel plot and multilevel Egger test |
| `FigS08_model_diagnostics.png` | Fig. S8: diagnostics of the primary model (Q–Q plots, Cook's distances, profile likelihoods) |

| Table | Content |
|---|---|
| `Table1_datasets.csv` | Table 1: analytical datasets |
| `Table2_subgroups_CI_PI.csv` | Table 2: mean R² by subgroup |
| `Table3_models_M1_M4.csv` | Table 3: models M1–M4 |
| `Table4_within_study_contrasts.csv` | Table 4: within-study contrasts |
| `TableS02_data_dictionary.csv` | Table S2: data dictionary |
| `TableS03_included_studies.csv` | Table S3: included studies |
| `TableS04_removed_values.csv` | Table S4: values not analysed, by reason |
| `TableS05_pooled_model.csv` | Table S5: pooled three-level model |
| `TableS06_estimand_classes.csv` | Table S6: pooled R² by estimand class |
| `TableS07_effect_size_strategies.csv` | Table S7: effect-size and weighting strategies |
| `TableS08_scale_and_estimand.csv` | Table S8: untransformed R² and Fisher's Z |
| `TableS09_M1_M4_two_scales.csv` | Table S9: models M1–M4 on both scales |
| `TableS10_moderators_M1_M4_two_scales.csv` | Table S10: moderator tests under M1–M4 on both scales |
| `TableS11_mechanisms.csv` | Table S11: mechanisms, one element changed at a time |
| `TableS12_moderator_robustness.csv` | Table S12: HTZ and cluster wild bootstrap moderator tests |
| `TableS13_moderator_tests.csv` | Table S13: robust tests of each moderator, Fisher's Z and log RMSE |
| `TableS14_moderator_confounding.csv` | Table S14: associations between moderators |
| `TableS15_meta_regression.csv` | Table S15: joint meta-regression |
| `TableS16_meta_regression_Wald.csv` | Table S16: Wald tests of the joint meta-regression |
| `TableS17_variance_explained.csv` | Table S17: heterogeneity explained |
| `TableS18_GEDI_vs_other_LiDAR.csv` | Table S18: GEDI vs other LiDAR sources |
| `TableS19_GEDI_roles_by_study.csv` | Table S19: GEDI roles by study |
| `TableS20_GEDI_products_verified.csv` | Table S20: GEDI products |
| `TableS22_contrast_precision.csv` | Table S22: precision of the within-study contrasts |
| `TableS23_contrasts_robustness.csv` | Table S23: robustness of the within-study contrasts |
| `TableS24_LiDAR_pairing.csv` | Table S24: LiDAR pairing by study |
| `TableS25_LiDAR_paired_summary.csv` | Table S25: unpaired and matched LiDAR estimates |
| `TableS26_LiDAR_matched_pairs.csv` | Table S26: matched LiDAR pairs by study |
| `TableS27_LiDAR_pairs_heterogeneity.csv` | Table S27: matched LiDAR pairs, pooled estimates and influence |
| `TableS28_LiDAR_pairs_sign_test.csv` | Table S28: matched LiDAR pairs, sign test |
| `TableS29_central_africa.csv` | Table S29: Central African studies |
| `TableS30_sensitivity.csv` | Table S30: sensitivity of the pooled estimate |
| `TableS30_leave_one_study_out.csv` | Table S30 (note): leave-one-study-out fits |
| `TableS30_leave_one_cluster_out.csv` | Table S30 (note): leave-one-cluster-out fits |
| `TableS31_sensitivity_additional.csv` | Table S31: additional sensitivity scenarios |
| `TableS32_clustering_levels.csv` | Table S32: clustering levels |
| `TableS33_multilevel_Egger.csv` | Table S33: multilevel Egger test |
| `TableS34_RMSE_scale_checks.csv` | Table S34: checks of the RMSE scale |
| `TableS35_model_diagnostics.csv` | Table S35: diagnostics of the primary model (per-cluster values in `TableS35_cluster_influence.csv`) |
| `FigS01_data_GEDI_test_verified.csv` | Fig. S1: GEDI test, unverified roles excluded |
| `FigS05_data_terrain_sample_size_reference.csv` | Fig. S5: estimates and tests |
| `FigS06_data_validation_resolution_region.csv` | Fig. S6: estimates |
| `All_tables.xlsx` | All tables, one sheet per table |

## Scales of analysis

Untransformed R² is the primary scale: every conclusion of the article rests on it, and every model, moderator test and regional test is fitted on it. Fisher's Z is complementary: for most values (multiple regression, machine-learning and predictive R²) atanh(√R²) is not variance-stabilising and 1/(n − 3) is not the sampling variance. A result is called robust only when it holds on both scales. Prediction intervals are reported on both scales, untransformed first: in the tables, columns `R2_PI_untransformed_*` give the primary interval and `R2_PI_lower`/`R2_PI_upper` the back-transformed Fisher's Z interval; in the figures, the thin bar or shaded band is the untransformed interval and the dotted bar or line the Fisher's Z interval.

## Documentation (`docs/`) and tests (`tests/`)

| File | Content |
|---|---|
| `docs/README.md` | This file. |
| `docs/data_dictionary.csv` | Variables of the databases and derived variables. |
| `tests/consistency_checks.R` | Automated checks (see step 4 above). |
| `tests/results_manifest.csv` | Each numerical result of the article checked by the tests, with its S2 table or figure, script, output file and source observations. |

## Reproducibility

- A complete re-execution of `run_all.R` from a copy of this folder without `output/` (Ubuntu 24.04, R 4.3.3, packages as in `renv.lock`) completed the 16 scripts without error; the regenerated tables were identical, cell by cell, to the submitted outputs (24 November 2026).
- `tests/consistency_checks.R` returned 33 PASS, 1 REVIEW and 0 FAIL, and 60 of 60 results of the article matched the outputs within half a unit of their last reported digit. The REVIEW item lists six pairs of observations with the same study, site, sensors, algorithm, training n and rounded R², to be confirmed as distinct models.
- Seeds: `set.seed()` in each script that resamples; each cluster wild bootstrap (9,999 draws) starts from `set.seed(20261123)`. Other resampling: cluster bootstrap 5,000 resamples, 500 random draws of one model per study.
- Not yet done: restoration of `renv.lock` on another operating system and verification by an independent analyst. To verify: install R 4.3.3, run `renv::restore()`, then `run_all.R` and `tests/consistency_checks.R`, and compare `output/tables/` with the archived outputs.

## Licence and citation

Please cite the article when using these data or scripts. All figures are drawn with ggplot2 from the data in `data/`; no templates or generative image tools are used.
