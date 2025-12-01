rm(list = ls())
gc()

library(ggplot2)
library(dplyr)
library(coda)
library(stringr)

'%notin%' <- Negate('%in%')

source("penguin_paths.R")

combine_logs <- function(path, pattern, extesion, burnin, resample=NULL){
  print(paste0(pattern, "_\\d.", extesion))
  logs <- list.files(path=path, pattern = paste0(pattern, "_[0-9]*\\.", extesion))
  cmd <- paste(logCombinerPath, "-b", burnin)
  if (!is.null(resample)){
    cmd <- paste(cmd, "-resample", resample)
  }
  for (log in logs){
    cmd <- paste(cmd, "-log", log)
  }
  cmd <- paste0(cmd, " -o ", pattern, ".combined")
  cmd <- paste0(cmd, ".", extesion)
  print(cmd)
  system(paste0("cd ",path ,";",cmd))
}

get_tree_length <- function(treePath){
  cmd <- paste(appLauncherPath, 'TreeStat2 stats.txt', treePath)
  system(cmd)
}




# combine and read logs for original FBD run
if (!file.exists(paste0(fbdOriginalPath, "/Mkv_morph+dna.combined.log"))){
  combine_logs(fbdOriginalPath, "Mkv_morph\\+dna", "log", burnin)
}
if (!file.exists(paste0(fbdOriginalPath, "/Mkv_morph+dna.combined.trees"))){
  combine_logs(fbdOriginalPath, "Mkv_morph\\+dna", "trees", burnin)
  get_tree_length(paste0(fbdOriginalPath, "/Mkv_morph+dna.combined.trees"))
}


# combine and read logs for original FBD run with updated dates
if (!file.exists(paste0(fbdUpdatedPath, "/Mkv_morph+dna.combined.log"))){
  combine_logs(fbdUpdatedPath, "Mkv_morph\\+dna", "log", burnin)
}
if (!file.exists(paste0(fbdUpdatedPath, "/Mkv_morph+dna.combined.trees"))){
  combine_logs(fbdUpdatedPath, "Mkv_morph\\+dna", "trees", burnin)
  get_tree_length(paste0(fbdUpdatedPath, "/Mkv_morph+dna.combined.trees"))
}

# combine and read logs for original sRanges run (also with updated dates)

# morph data at both ends
if (!file.exists(paste0(sRangesBothPath, "/penguins_inf_morph_at_both.combined.log"))){
  combine_logs(sRangesBothPath, "penguins_inf_morph_at_both", "log", burnin)
}
if (!file.exists(paste0(sRangesBothPath, "/penguins_inf_morph_at_both.combined.speciation.log"))){
  combine_logs(sRangesBothPath, "penguins_inf_morph_at_both", "speciation.log", burnin)
}
if (!file.exists(paste0(sRangesBothPath, "/penguins_inf_morph_at_both.combined.trees"))){
  combine_logs(sRangesBothPath, "penguins_inf_morph_at_both", "trees", burnin)
  get_tree_length(paste0(sRangesBothPath, "/penguins_inf_morph_at_both.combined.trees"))
}


# morph data at start
if (!file.exists(paste0(sRangesPath, "/penguins_inf_morph_at_start.combined.log"))){
  combine_logs(sRangesPath, "penguins_inf_morph_at_start", "log", burnin)
}
if (!file.exists(paste0(sRangesPath, "/penguins_inf_morph_at_start.combined.speciation.log"))){
  combine_logs(sRangesPath, "penguins_inf_morph_at_start", "speciation.log", burnin)
}
if (!file.exists(paste0(sRangesPath, "/penguins_inf_morph_at_start.combined.trees"))){
  combine_logs(sRangesPath, "penguins_inf_morph_at_start", "trees", burnin)
  get_tree_length(paste0(sRangesPath, "/penguins_inf_morph_at_start.combined.trees"))
}

if (!file.exists(paste0(sRangesBothPath,"/penguins_inf_morph_at_both.tipAge.log"))){
  cmd <- paste0(beast_app_dir,"beast relog_tip_age.xml")
  system(paste0("cd ",sRangesBothPath ,";",cmd))
  system()
}

if (!file.exists(paste0(sRangesPath,"/penguins_inf_morph_at_start.tipAge.log"))){
  cmd <- paste0(beast_app_dir,"beast relog_tip_age.xml")
  system(paste0("cd ",sRangesPath ,";",cmd))
  system()
}



