#######################################################
######### plots to compare distributions      #########
######### of parameters and some tree         #########
######### statistics between FBD and sRanges  #########
#######################################################
rm(list = ls())
gc()

source("penguin_helper.R")
figure_path <- paste0(figure4_tableA6_dir,"/")


table1 <- read.table(paste0(fbdOriginalPath, "/Mkv_morph+dna.combined.log"), sep = "\t", header = TRUE)

table2 <- read.table(paste0(fbdUpdatedPath, "/Mkv_morph+dna.combined.log"), sep = "\t", header = TRUE)

table3 <- read.table(paste0(sRangesBothPath, "/penguins_inf_morph_at_both.combined.log"), sep = "\t", header = TRUE)

table4 <- read.table(paste0(sRangesPath, "/penguins_inf_morph_at_start.combined.log"), sep = "\t", header = TRUE)


if ("TreeLength" %notin% colnames(table1)){
  table1$TreeLength <- as.numeric(unlist(read.csv(paste0(fbdOriginalPath,"/Mkv_morph+dna.combined.treestreestats.log"),
                                sep="\t")["Tree.Length"]))
}
if ("TreeLength" %notin% colnames(table2)){
  table2$TreeLength <- as.numeric(unlist(read.csv(paste0(fbdUpdatedPath,"/Mkv_morph+dna.combined.treestreestats.log"),
                                sep="\t")["Tree.Length"]))
}
if ("TreeLength" %notin% colnames(table3)){
  table3$TreeLength <- as.numeric(unlist(read.csv(paste0(sRangesBothPath, "/penguins_inf_morph_at_both.combined.treestreestats.log"),
                                sep="\t")["Tree.Length"]))
}
if ("TreeLength" %notin% colnames(table4)){
  table4$TreeLength <- as.numeric(unlist(read.csv(paste0(sRangesPath,"/penguins_inf_morph_at_start.combined.treestreestats.log"),
                                sep="\t")["Tree.Length"]))
}


# div rate prior
div.prior <- rlnorm(1000, meanlog=-3.5, sdlog=1.5)
turnover.prior <- runif(1000, min=0, max=1)
sampling.prior <- runif(1000, min=0, max=1)


colnames(table1) <- str_remove_all(colnames(table1), "SABD")
colnames(table2) <- str_remove_all(colnames(table2), "SABD")
colnames(table3) <- str_remove_all(colnames(table3), "SRFBD")
colnames(table4) <- str_remove_all(colnames(table4), "SRFBD")
table1$model <- "FBD, orig."
table2$model <- "FBD, upd."
table3$model <- "SRFBD, both"
table4$model <- "SRFBD, first"

# Combine the tables into a single data frame
# combined_df_all <- bind_rows(table1, table2, table3, .id = "Table") # all 3 runs
combined_df <- bind_rows(table1, table2, table3, table4) # only updated date

# Get the unique matching column names (params and statistics that are logged for all runs)
# matching_columns_all <- intersect(intersect(colnames(table1), colnames(table2)), colnames(table3))
matching_columns_ <- Reduce(intersect, list(colnames(table1), colnames(table2), colnames(table3), colnames(table4)))[]
matching_columns <- c(matching_columns_[6:21], matching_columns_[59])
rm(matching_columns_)
# combined_df_all$Table <- as.factor(combined_df_all$Table)
combined_df$model <- as.factor(combined_df$model)

# hpd <- data.frame(name=matching_columns_all[6:26], hpd.low.fbd=numeric(21), hpd.high.fbd=numeric(21), hpd.width.fbd=numeric(21),
#                   hpd.low.fbdUp=numeric(21), hpd.high.fbdUp=numeric(21), hpd.width.fbdUp=numeric(21),
#                   hpd.low.srfbd=numeric(21), hpd.high.srfbd=numeric(21), hpd.width.srfbd=numeric(21))

hpd_fbd_orig <- data.frame(name=matching_columns, median=numeric(17), 
                           hpd.low=numeric(17), hpd.high=numeric(17), hpd.width=numeric(17))
hpd_fbd_upd <- hpd_fbd_orig
hpd_srfbd_first <- hpd_fbd_orig
hpd_srfbd_both <- hpd_fbd_orig
hpd_frames <- list(hpd_fbd_orig, hpd_fbd_upd, hpd_srfbd_both, hpd_srfbd_first)
tables <- list(table1, table2, table3, table4)
names <- c("hpd_fbd_orig","hpd_fbd_upd" , "hpd_srfbd_both","hpd_srfbd_first")

addToHPD<- function(dataFrame, table, colName, i){
  h <- HPDinterval(as.mcmc(table[,which(colnames(table)==colName)]))
  dataFrame[i, 2]<-median(table[,which(colnames(table)==colName)])
  dataFrame[i, 3]<-h[1]
  dataFrame[i, 4]<-h[2]
  dataFrame[i, 5]<-h[2]-h[1]
  return(dataFrame)
}

i=1
#"FBD, \nUpdated Dates" = "#FF7F50", "FBD" = "#B38CB4"
# Create parameter and tree statistics plots comparing all 3 runs
for (col in matching_columns){
  
  # Create the violin plot for the matching column name
  plot <- ggplot(combined_df, aes_string(x = "model", y = col, fill = "model", color="model")) +
    geom_violin(alpha=0.8) +
    stat_summary(fun = "median", geom = "crossbar", width = 0.25,
                 position = position_dodge(width = .25), color = "#5A5A5A", alpha=0.8, show.legend = F) +
    labs(x = "", y = "") +
    ggtitle("") +
    theme_minimal() +
    scale_fill_manual(values = c("FBD, orig."="#913615","FBD, upd."="#FF7F50",
                                 "SRFBD, both"="#255957", "SRFBD, first"="#43AA8B"),
                      labels = c("FBD, original dates",  "FBD, updated dates",
                                 "SRFBD morph. data at both bounds",
                                 "SRFBD morph. data at first occurence"))+
    scale_colour_manual(values = c("FBD, orig."="#913615","FBD, upd."="#FF7F50",
                                   "SRFBD, both"="#255957", "SRFBD, first"="#43AA8B"), 
                        labels = c("FBD, original dates",  "FBD, updated dates",
                                   "SRFBD morph. data at both bounds",
                                   "SRFBD morph. data at first occurence"))+
    theme(
      legend.position = "none",
      legend.direction = "horizontal",
      legend.title = element_blank(),
      text = element_text(size = 24),
      plot.title = element_text(size = 30, face = "bold"),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 20),
    )
  if (col=="rho"){
    plot <- plot +
      scale_x_discrete(labels = c())  
    ggsave(paste0(figure_path,"penguin_",col,figure_extension), plot, height=4, width=7)
  } else{
    plot <- plot +
      scale_x_discrete(labels = c())
    ggsave(paste0(figure_path,"penguin_",col,figure_extension), plot, height=4, width=7)
  }
  
  for (j in 1:length(hpd_frames)){
    hpd_frames[[j]] <- addToHPD(hpd_frames[[j]], tables[[j]], col, i)
  }
  i=i+1
}


## to get the legend for the above plot use this:
# my_legend <- get_legend(plot)
# as_ggplot(my_legend)
# ggsave(paste0(figure_path,"penguin_legend",figure_extension), my_legend, height=0.5, width=16)

##########################################
######## Output for the HPD table ########
##########################################

for (j in 1:length(hpd_frames)){
  hpd<- hpd_frames[[j]]
  hpd[,2:ncol(hpd)] <- round(hpd[,2:ncol(hpd)], 3)
  write.csv(hpd, paste0(figure_path, names[j],".csv"))
}


# table4 <- data.frame()
# 
# # Combine the tables into a single data frame
# combined_df_all <- bind_rows(table1, table2, table3, .id = "Table") # all 3 runs
# combined_df <- bind_rows(table2, table3, .id = "Table") # only updated date
# 
# # Get the unique matching column names (params and statistics that are logged for all runs)
# matching_columns_all <- intersect(intersect(colnames(table1), colnames(table2)), colnames(table3))
# matching_columns <- intersect(colnames(table2), colnames(table3))
# 
# combined_df_all$Table <- as.factor(combined_df_all$Table)
# combined_df$Table <- as.factor(combined_df$Table)
# 
# hpd <- data.frame(name=matching_columns_all[6:26], median.fbd=numeric(21), hpd.low.fbd=numeric(21), hpd.high.fbd=numeric(21), hpd.width.fbd=numeric(21),
#                   median.fbdUp=numeric(21), hpd.low.fbdUp=numeric(21), hpd.high.fbdUp=numeric(21), hpd.width.fbdUp=numeric(21),
#                   median.srfbd=numeric(21), hpd.low.srfbd=numeric(21), hpd.high.srfbd=numeric(21), hpd.width.srfbd=numeric(21))
# i=1
# # Create parameter and tree statistics plots comparing all 3 runs
# for (col in matching_columns_all[6:26]){
#   
#   # Create the violin plot for the matching column name
#   plot <- ggplot(combined_df_all, aes_string(x = "Table", y = col, fill = "Table", color="Table")) +
#     geom_violin(alpha=0.8) +
#     stat_summary(fun = "median", geom = "crossbar", width = 0.25,
#                  position = position_dodge(width = .25), color = "#5A5A5A", alpha=0.8) +
#     labs(x = "", y = "") +
#     ggtitle("") +
#     theme_minimal() +
#     scale_x_discrete(labels = c("FBD", "FBD, Upd.",  "SRFBD")) +
#     scale_fill_manual(values = c("#B38CB4", "#FF7F50", "#43AA8B"))+
#     scale_colour_manual(values = c("#B38CB4", "#FF7F50",  "#43AA8B"))+
#     theme(
#       legend.position = "none",
#       text = element_text(size = 24),
#       plot.title = element_text(size = 30, face = "bold"),
#       axis.title = element_text(size = 24),
#       axis.text = element_text(size = 20),
#     )
#   
#   ggsave(paste0(figure_path,"penguin_",col,figure_extension), plot, width=7.5, height=6)
#   
#   h <- HPDinterval(as.mcmc(table1[,which(colnames(table1)==col)]))
#   hpd[i, 2]<-median(table1[,which(colnames(table1)==col)])
#   hpd[i, 3]<-h[1]
#   hpd[i, 4]<-h[2]
#   hpd[i, 5]<-h[2]-h[1]
#   
#   h <- HPDinterval(as.mcmc(table2[,which(colnames(table2)==col)]))
#   hpd[i, 6]<-median(table2[,which(colnames(table2)==col)])
#   hpd[i, 7]<-h[1]
#   hpd[i, 8]<-h[2]
#   hpd[i, 9]<-h[2]-h[1]
#   
#   h <- HPDinterval(as.mcmc(table3[,which(colnames(table3)==col)]))
#   hpd[i, 10]<-median(table3[,which(colnames(table3)==col)])
#   hpd[i, 11]<-h[1]
#   hpd[i, 12]<-h[2]
#   hpd[i, 13]<-h[2]-h[1]
#   i=i+1
# }
# 
# hpd[,2:ncol(hpd)] <- round(hpd[,2:ncol(hpd)], 3)
# write.csv(hpd, paste0(figure_path, "hpd.csv"))

# # Create parameter and tree statistics plots comparing only FBD and sRanges plots with updated dates
# for (col in matching_columns[6:26]){
# 
#   # Create the violin plot for the matching column name
#   plot <- ggplot(combined_df, aes_string(x = "Table", y = col, fill = "Table", color="Table")) +
#     geom_violin(alpha=0.8) +
#     stat_summary(fun = "median", geom = "crossbar", width = 0.25,
#                  position = position_dodge(width = .25), color = "#5A5A5A", alpha=0.8) +
#     labs(x = "", y = "") +
#     ggtitle("") +
#     theme_minimal() +
#     scale_x_discrete(labels = c("FBD, \nUpdated Dates", "SRFBD")) +
#     scale_fill_manual(values = c("#FF7F50", "#43AA8B"))+
#     scale_colour_manual(values = c("#FF7F50", "#43AA8B"))+
#     theme(
#       legend.position = "none",
#       text = element_text(size = 24),
#       plot.title = element_text(size = 30, face = "bold"),
#       axis.title = element_text(size = 24),
#       axis.text = element_text(size = 20),
#     )
#   
#   ggsave(paste0(figure_path,col,"_updatedDates",figure_extension), plot)
# }



 

