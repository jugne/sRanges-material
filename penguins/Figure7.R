rm(list = ls())
gc()

source("penguin_paths.R")

cmd <- paste0(appLauncherPath, " SRangesAndSACladeSetComparator -SATree ",
             fbdUpdatedPath,"/Mkv_morph+dna.combined.trees -sRangesTre",
             sRangesPath, "penguins_inf_morph_at_start.combined.trees  -out ",
             figure7_dir,"/clade_comp_morph_start.txt -png clade_comp_morph_start.png 
             -taxonsetSA Burnside_Palaeudyptes,Waimanu_tuatahi 
             -taxonsetSRanges Burnside_Palaeeudyptes,Muriwaimanu_tuatahi 
             -burnin 0")

system(cmd)

cmd <- paste0(appLauncherPath, " SRangesAndSACladeSetComparator -SATree ",
              fbdUpdatedPath,"/Mkv_morph+dna.combined.trees -sRangesTre",
              sRangesBothPath, "penguins_inf_morph_at_both.combined.trees  -out ",
              figure7_dir,"/clade_comp_morph_both.txt -png clade_comp_morph_both.png 
             -taxonsetSA Burnside_Palaeudyptes,Waimanu_tuatahi 
             -taxonsetSRanges Burnside_Palaeeudyptes,Muriwaimanu_tuatahi 
             -burnin 0")

system(cmd)
