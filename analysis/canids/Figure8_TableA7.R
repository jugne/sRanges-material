#######################################################
######### plots to compare distributions      #########
######### of parameters and some tree         #########
######### statistics between FBD and sRanges  #########
#######################################################
rm(list = ls())
gc()

source("penguin_helper.R")
figure_path <- paste0(figure8_tableA7_dir,"/")


# read logs for original FBD run
table1 <- read.table(paste0(fbdOriginalPath, "/canids_fbd_fixed.combined.log"), sep = "\t", header = TRUE)

# read logs for original FBD run with updated dates
table2 <- read.table(paste0(fbdUpdatedPath, "/canids_zero_offset_fbd.combined.log"), sep = "\t", header = TRUE)

# read logs for original sRanges run, morph data at both ends (also with updated dates)
table3 <- read.table(paste0(sRangesBothPath, "/canids_rho.combined.log"), sep = "\t", header = TRUE)

# read logs for original sRanges run, morph data at the start (also with updated dates)
table4 <- read.table(paste0(sRangesPath, "/canids_rho.combined.log"), sep = "\t", header = TRUE)


# added the tree length that we forgot to log but got when combining trees

if ("TreeLength" %notin% colnames(table1)){
  table1$TreeLength <- as.numeric(unlist(read.csv(paste0(fbdOriginalPath,"canids_fbd_fixed.length.log"),
                                                  sep="\t")["tree.treeLength"]))
}
if ("TreeLength" %notin% colnames(table2)){
  table2$TreeLength <- as.numeric(unlist(read.csv(paste0(fbdUpdatedPath,"canids_zero_offset_fbd.length.log"),
                                                  sep="\t")["tree.treeLength"]))
}
if ("TreeLength" %notin% colnames(table3)){
  table3$TreeLength <- as.numeric(unlist(read.csv(paste0(sRangesBothPath, "canidae_rho.length.log"),
                                                  sep="\t")["tree.treeLength"]))
}
if ("TreeLength" %notin% colnames(table4)){
  table4$TreeLength <- as.numeric(unlist(read.csv(paste0(sRangesPath,"canidae_rho.length.log"),
                                                  sep="\t")["tree.treeLength"]))
}


colnames(table1) <- str_remove_all(colnames(table1), "SABD")
colnames(table2) <- str_remove_all(colnames(table2), "SABD")
colnames(table3) <- str_remove_all(colnames(table3), "SRFBD")
colnames(table4) <- str_remove_all(colnames(table4), "SRFBD")

table1$model <- "FBD, orig."
table2$model <- "FBD, upd."
table3$model <- "SRFBD, both"
table4$model <- "SRFBD, first"


# Combine the tables into a single data frame
combined_df <- bind_rows(table1, table2, table3, table4)

# Get the unique matching column names (params and statistics that are logged for all runs)
# matching_columns_all <- intersect(intersect(colnames(table1), colnames(table2)), colnames(table3))
matching_columns <- Reduce(intersect, list(colnames(table1), colnames(table2), colnames(table3), colnames(table4)))

# combined_df_all$Table <- as.factor(combined_df_all$Table)
combined_df$model <- as.factor(combined_df$model)

hpd_fbd_orig <- data.frame(name=matching_columns[6:22], median=numeric(17), 
                           hpd.low=numeric(17), hpd.high=numeric(17), hpd.width=numeric(17))
hpd_fbd_upd <- hpd_fbd_orig
hpd_srfbd_first <- hpd_fbd_orig
hpd_srfbd_both <- hpd_fbd_orig
hpd_frames <- list(hpd_fbd_orig, hpd_fbd_upd, hpd_srfbd_first, hpd_srfbd_both)
tables <- list(table1, table2, table3, table4)
names <- c("hpd_fbd_orig","hpd_fbd_upd", "hpd_srfbd_first", "hpd_srfbd_both")

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
for (col in matching_columns[6:22]){
  
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
      legend.title = element_blank(),
      text = element_text(size = 24),
      plot.title = element_text(size = 30, face = "bold"),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 20),
    )
  if (col=="rho"){
    plot <- plot +
      scale_x_discrete(labels = c())  
    ggsave(paste0(figure_path,"canid_",col,figure_extension), plot, height=4, width=7)
  } else{
    plot <- plot +
      scale_x_discrete(labels = c())
    ggsave(paste0(figure_path,"canid_",col,figure_extension), plot, height=4, width=7)
  }
  
  for (j in length(hpd_frames)){
    hpd_frames[[j]]<-addToHPD(hpd_frames[[j]], tables[[j]], col, i)
  }
  i=i+1
}

## to get the legend for the above plot use this:
# my_legend <- get_legend(plot)
# as_ggplot(my_legend)
# ggsave(paste0(figure_path,"canid_legend",figure_extension), my_legend, height=4, width=5)


for (j in length(hpd_frames)){
  hpd<- hpd_frames[[j]]
  hpd[,2:ncol(hpd)] <- round(hpd[,2:ncol(hpd)], 3)
  write.csv(hpd, paste0(figure_path, "canid_", names[j],".csv"))
}


getSATable<- function(treeFile, outFile, analyse=T){
  if (analyse){
    cmd <- paste(appLauncherPath,'SampledAncestorTreeAnalyser -file', treeFile, ">", outFile)
    system(cmd)
  }
  # Determine the total number of lines in the file
  # total_lines <- length(readLines(outFile))
  # Skip the first two lines and read until (total_lines - 2)
  data <- read.table(outFile, sep = "\t", header = TRUE, skip = 2)#, nrows = total_lines - 3)
  
  # old_names <- c("Burnside_Palaeudyptes", "Waimanu_tuatahi")
  # new_names <- c("Burnside_Palaeeudyptes", "Muriwaimanu_tuatahi")
  # for (i in 1:length(old_names)){
  #   id<-which(data$SA==old_names[i])
  #   if (length(id)!=0){
  #     data$SA[id] <- new_names[i]
  #   }
  # }
  data$SA <- gsub("_last", "", gsub("_first", "", data$SA))
  data$SA <- gsub("_", " ", data$SA)
  n <- names(which(table(data$SA)==2))
  for (nn in n){
    ids <- which(data$SA==nn) 
    data <- data[-ids[1],]
  }
  return(data)
}
