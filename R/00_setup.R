cat("=== AGB Meta-Analysis Setup ===\n")
cat("R version required: >= 4.2.0\n")
cat("Current R version:", R.version$major, ".", R.version$minor, "\n", sep = "")

pkgs <- c(
  "metafor",
  "robumeta",
  "dplyr",
  "tidyr",
  "ggplot2",
  "patchwork",
  "scales",
  "readr",
  "openxlsx",
  "grid",
  "ggrepel",
  "clubSandwich",
  "sf",
  "maps",
  "stringr",
  "nlme"
)

installed <- rownames(installed.packages())
to_install <- pkgs[!pkgs %in% installed]

if (length(to_install) > 0) {
  cat("\nInstalling missing packages:", paste(to_install, collapse = ", "), "\n")
  type_inst <- if (.Platform$OS.type == "windows") "binary" else getOption("pkgType")
  install.packages(to_install, repos = "https://cloud.r-project.org", type = type_inst, quiet = TRUE)
} else {
  cat("\nAll packages already installed.\n")
}

not_loaded <- c()
suppressPackageStartupMessages({
  for (pkg in pkgs) {
    ok <- suppressWarnings(require(pkg, character.only = TRUE, quietly = TRUE))
    if (!ok) not_loaded <- c(not_loaded, pkg)
  }
})
if (length(not_loaded)) {
  warning("Packages not available: ", paste(not_loaded, collapse = ", "),
          " — scripts that need them will fail.")
} else cat("All packages loaded successfully.\n\n")

set.seed(2024)

options(
  scipen    = 999,
  digits    = 4,
  warn      = 0
)

SENSOR_COLORS <- c(
  "LiDAR+Optical"       = "#2196F3",
  "LiDAR-only"          = "#009688",
  "Optical+SAR"         = "#9C27B0",
  "SAR+LiDAR"           = "#FF9800",
  "Optical-only"        = "#03A9F4",
  "SAR-only"            = "#F44336",
  "Optical+SAR+LiDAR"   = "#212121"
)

MODEL_COLORS <- c(
  "Statistical"         = "#2196F3",
  "Machine-learning"    = "#4CAF50",
  "Deep-learning"       = "#E91E63"
)

CLIMATE_COLORS <- c(
  "Mediterranean"       = "#2196F3",
  "Subtropical"         = "#FF9800",
  "Temperate"           = "#9C27B0",
  "Arid"                = "#FFEB3B",
  "Tropical"            = "#F44336",
  "Boreal"              = "#90CAF9"
)

GEDI_COLORS <- c("Yes" = "#009688", "No" = "#FF7043")

SCALE_COLORS <- c(
  "Local"         = "#4CAF50",
  "Regional"      = "#2196F3",
  "National"      = "#FF9800",
  "Multi-country" = "#E91E63"
)

FOREST_COLORS <- c(
  "Mixed forest"          = "#8D6E63",
  "Tropical moist forest" = "#F44336",
  "Plantation"            = "#FF9800",
  "Shrubland"             = "#FFEB3B",
  "Coniferous forest"     = "#1B5E20",
  "Deciduous forest"      = "#66BB6A",
  "Evergreen_forest"      = "#00796B",
  "Mangrove"              = "#0288D1",
  "Miombo woodland"       = "#CE93D8",
  "Dry forest"            = "#FF8A65",
  "Flooded forest"        = "#80DEEA"
)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)
CACHE_DIR <- file.path(tempdir(), "agb_cache"); dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

z_to_r <- function(z) tanh(z)

z_to_r2 <- function(z) tanh(z)^2

r_to_z <- function(r) 0.5 * log((1 + r) / (1 - r))

se_z <- function(n) 1 / sqrt(n - 3)

cohen_magnitude <- function(r) {
  ifelse(abs(r) >= 0.50, "LARGE",
    ifelse(abs(r) >= 0.30, "MEDIUM", "SMALL"))
}

epsilon_squared <- function(H, n) H / (n - 1)

cat("Helper functions defined.\n")
cat("Output directories: output/figures/, output/tables/\n")
cat("Global seed: set.seed(2024)\n\n")
cat("=== Setup complete. Run 01_data_import.R next. ===\n")

print(sessionInfo())
