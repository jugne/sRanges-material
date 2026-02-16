rm(list = ls())
gc()

source("canid_helper.R")


##############################
############ FBD #############
##############################

# combine and read logs for original FBD run
combine_and_get_length(fbdOriginalPath, "canids_fbd_fixed", burnin)

# combine and read logs for original FBD run with updated dates
combine_and_get_length(fbdUpdatedPath, "canids_zero_offset_fbd", burnin)


##############################
############ SRFBD ###########
##############################

# combine and read logs for sRanges with updated dates

# morph data at both ends
combine_and_get_length(sRangesBothPath, "canids_rho", burnin, "srfbd")


# morph data at start
combine_and_get_length(sRangesPath, "canids_rho", burnin, "srfbd")


# relog the tip ages
relog_tip_ages(sRangesBothPath, "canids_rho")

relog_tip_ages(sRangesPath, "canids_rho")




