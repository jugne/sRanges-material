library(ggplot2) 
library(treeio)
library(ggtree)
library(ape)
library(tidytree)


# Define base directories
base_dir <- getwd()
ggtree_source_dir <- "~/Documents/Source/ggtree/R"
beast_app_dir <- "'/Applications/BEAST 2.7.7/bin'"

# Define application paths
appLauncherPath <- file.path(beast_app_dir, "applauncher")
logCombinerPath <- file.path(beast_app_dir, "logcombiner")

# Define specific input file paths
fbdOriginalPath <- file.path(base_dir, "Gavryushkina_et_al/original_rerun")
fbdUpdatedPath <- file.path(base_dir, "Gavryushkina_et_al/updated_dates_rerun")
sRangesPath <- file.path(base_dir, "morph_at_start")
sRangesBothPath <- file.path(base_dir, "morph_at_both")

# Define output figure paths
figure6_dir <- file.path(base_dir, "figures/Figure6")
figure4_tableA6_dir <- file.path(base_dir, "figures/Figure4_TableA6")
figure13a_dir <- file.path(base_dir, "figures/Figure13a")
figure7_dir <- file.path(base_dir, "figures/Figure7")
table2_figure5a_dir <- file.path(base_dir, "figures/Table2_Figure5a")
table3_4_dir <- file.path(base_dir, "figures/Tables3_4")
tableA45_dir <- file.path(base_dir, "figures/TableA5")

# Create necessary directories if they do not exist
dir.create(figure6_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure4_tableA6_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure13a_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table2_figure5a_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table3_4_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tableA45_dir, recursive = TRUE, showWarnings = FALSE)

# Other constants
figure_extension <- ".pdf"
sa_percent <- 5
burnin <- 10

# Source all R files in ggtree (if needed)
if (dir.exists(ggtree_source_dir)) {
  setwd(ggtree_source_dir)
  file.sources <- list.files(pattern = "*.R")
  sapply(file.sources, source, .GlobalEnv)
}

# Set back to base_dir as working directory
setwd(base_dir)
