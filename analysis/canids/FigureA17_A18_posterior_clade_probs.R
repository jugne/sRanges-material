rm(list = ls())
gc()



source("canid_helper.R")

figure_path <- paste0(figureA17_A18dir, "/")
plotPosteriorCladeAge("canid", figure7cd_dir)