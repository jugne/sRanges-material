#######################################################
##### PENGUIN: helper functions and common setup  #####
#######################################################

source("../helper.R")


# Define base directories
base_dir <- getwd()

# Define specific input file paths
fbdOriginalPath <- file.path(base_dir, "Gavryushkina_et_al/original_rerun")
fbdUpdatedPath <- file.path(base_dir, "Gavryushkina_et_al/updated_dates_rerun")
sRangesPath <- file.path(base_dir, "morph_at_start")
sRangesBothPath <- file.path(base_dir, "morph_at_both")

# Define output figure paths
figure6_dir <- file.path(base_dir, "figures/Figure6")
figure4_tableA6_dir <- file.path(base_dir, "figures/Figure4_TableA6")
figureA13a_dir <- file.path(base_dir, "figures/FigureA13a")
figure7ab_dir <- file.path(base_dir, "figures/Figure7ab")
figure5a_table2_dir <- file.path(base_dir, "figures/Figure5a_Table2")
figureA15_A16dir <- file.path(base_dir, "figures/FigureA15_A16")

table3_4_dir <- file.path(base_dir, "figures/Tables3_4")
tableA5_dir <- file.path(base_dir, "figures/TableA5")

# Create necessary directories if they do not exist
dir.create(figure6_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure4_tableA6_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figureA13a_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure7ab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure5a_table2_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figureA15_A16dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table3_4_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tableA5_dir, recursive = TRUE, showWarnings = FALSE)


# Set back to base_dir as working directory
setwd(base_dir)
