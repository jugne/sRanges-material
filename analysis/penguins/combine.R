rm(list = ls())
gc()
options(scipen = 999)

source("penguin_helper.R")


##############################
############ FBD #############
##############################

# combine and read logs for original FBD run
combine_and_get_length(fbdOriginalPath, "Mkv_morph+dna", burnin)

# combine and read logs for original FBD run with updated dates
combine_and_get_length(fbdUpdatedPath, "Mkv_morph+dna", burnin)


##############################
############ SRFBD ###########
##############################

# combine and read logs for sRanges with updated dates

# morph data at both ends
combine_and_get_length(sRangesBothPath, "penguins_inf_morph_at_both", burnin, model="srfbd")


# morph data at start
combine_and_get_length(sRangesPath, "penguins_inf_morph_at_start", burnin, model="srfbd", resample=150000)



# relog the tip ages
relog_tip_ages(sRangesBothPath, "tipAgeRelog/penguins_inf_morph_at_both")

relog_tip_ages(sRangesPath, "tipAgeRelog/penguins_inf_morph_at_start")



