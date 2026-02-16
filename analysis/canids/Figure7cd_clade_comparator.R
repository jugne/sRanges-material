rm(list = ls())
gc()

source("canid_helper.R")

cmd <- paste0(appLauncherPath, " SRangesAndSACladeSetComparator -SATree ",
              fbdUpdatedPath,"/canidae_rho_zero.combined.trees -SALog ",
              fbdUpdatedPath,"/canidae_rho_zero.combined.log -sRangesTree ",
              sRangesPath, "/canidae_rho.combined.trees  -out ",
              figure7cd_dir,"/clade_comp_morph_start.txt -png ",
              figure7cd_dir,"/clade_comp_morph_start.png ",
             "-burnin 0")

system(cmd)

cmd <- paste0(appLauncherPath, " SRangesAndSACladeSetComparator -SATree ",
              fbdUpdatedPath,"/canidae_rho_zero.combined.trees -SALog ",
              fbdUpdatedPath,"/canidae_rho_zero.combined.log -sRangesTree ",
              sRangesBothPath, "/canidae_rho.combined.trees  -out ",
              figure7cd_dir,"/clade_comp_morph_both.txt -png ",
              figure7cd_dir,"/clade_comp_morph_both.png ",
             "-burnin 0")

system(cmd)