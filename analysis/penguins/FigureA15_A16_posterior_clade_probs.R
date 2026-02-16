rm(list = ls())
gc()



source("penguin_helper.R")

figure_path <- paste0(figureA15_A16dir, "/")
plotPosteriorCladeAge("penguin", figure7ab_dir)