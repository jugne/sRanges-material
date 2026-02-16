#####################################################
##### CANID: helper functions and common setup  #####
#####################################################

source("../helper.R")

##### SET PATHS ######

# Define base directories
base_dir <- getwd()

# Define specific input file paths
fbdOriginalPath <- file.path(base_dir, "FBD_rho_fixed_Slater_dates")
fbdUpdatedPath <- file.path(base_dir, "canids/FBD_zero_offset_new_dates_extant_with_range_till_present")
sRangesPath <- file.path(base_dir, "SRFBD_rho_new_dates")
sRangesBothPath <- file.path(base_dir, "SRFBD_rho_new_dates_morph_first_last")


# Define output figure paths
figures5b_A14_table2_dir <- file.path(base_dir, "figures/Figures5b-A14_Table2")
figure7cd_dir <- file.path(base_dir, "figures/Figure7cd")
figure8_tableA7_dir <- file.path(base_dir, "figures/Figure8_TableA7")
figure9_dir <- file.path(base_dir, "figures/Figure9")
figureA13b_dir <- file.path(base_dir, "figures/FigureA13b")
figureA14_dir <- file.path(base_dir, "figures/FigureA14")
figureA17_A18dir <- file.path(base_dir, "figures/FigureA17_A18")



# Create necessary directories if they do not exist
dir.create(figures5b_A14_table2_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure7cd_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure8_tableA7_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure9_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figureA13b_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figureA14_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figureA17_A18dir, recursive = TRUE, showWarnings = FALSE)
