rm(list = ls())
gc()

source("penguin_helper.R")

cmd <- paste0(appLauncherPath, " SRangesAndSACladeSetComparator -SATree ",
             fbdUpdatedPath,"/Mkv_morph+dna.combined.trees -sRangesTree ",
             sRangesPath, "/penguins_inf_morph_at_start.combined.trees  -out ",
             figure7ab_dir,"/clade_comp_morph_start.txt -png ",
             figure7ab_dir,"/clade_comp_morph_start.png ",
             "-taxonsetSA Burnside_Palaeudyptes,Waimanu_tuatahi ",
             "-taxonsetSRanges Burnside_Palaeeudyptes,Muriwaimanu_tuatahi ",
             "-burnin 0")

system(cmd)

cmd <- paste0(appLauncherPath, " SRangesAndSACladeSetComparator -SATree ",
              fbdUpdatedPath,"/Mkv_morph+dna.combined.trees -sRangesTree ",
              sRangesBothPath, "/penguins_inf_morph_at_both.combined.trees  -out ",
              figure7ab_dir,"/clade_comp_morph_both.txt -png ",
              figure7ab_dir,"/clade_comp_morph_both.png ",
             "-taxonsetSA Burnside_Palaeudyptes,Waimanu_tuatahi ",
             "-taxonsetSRanges Burnside_Palaeeudyptes,Muriwaimanu_tuatahi ",
             "-burnin 0")

system(cmd)
