#######################################################
######### plots to compare distributions      #########
######### of parameters and some tree         #########
######### statistics between FBD and sRanges  #########
#######################################################

rm(list = ls())
library(ggplot2)
library(dplyr)
library(coda)
library(stringr)

'%notin%' <- Negate('%in%')

##########################################################
## Make sure the SA (FBD) package is installed in BEAST ##
##########################################################

##### SET PATHS ######
appLauncherPath <- "'/Applications/BEAST 2.6.7/bin/applauncher'"
logCombinerPath <- "'/Applications/BEAST 2.6.7/bin/logcombiner'"
fbdOriginalPath <- "~/Documents/Source/beast2.7/sRanges-material/canids/FBD_rho_fixed_Slater_dates/"
fbdUpdatedPath <- "~/Documents/Source/beast2.7/sRanges-material/canids/FBD_zero_offset_new_dates_extant_with_range_till_present_/"
sRangesBothPath <- "~/Documents/Source/beast2.7/sRanges-material/canids/SRFBD_rho_new_dates_morph_first_last1/"
sRangesPath <- "~/Documents/Source/beast2.7/sRanges-material/canids/SRFBD_rho_new_dates1_/"
figure_path <- "~/Documents/Source/beast2.7/sRanges-material/canids/figures/"
figure_extension <- ".pdf" # e.g., "pdf", "png", "svg"


combine_logs <- function(path, pattern, extesion, burnin, resample){
  print(paste0(pattern, "_\\d.", extesion))
  logs <- list.files(path=path, pattern = paste0(pattern, "\\.", extesion))
  cmd <- paste(logCombinerPath, "-b", burnin)
  for (log in logs){
    cmd <- paste(cmd, "-log", log)
  }
  cmd <- paste0(cmd, " -o ", pattern, ".combined")
  cmd <- paste0(cmd, ".", extesion)
  if (resample !=0){
    cmd <- paste0(cmd, " -resample ", resample)
  }
  print(cmd)
  system(paste0("cd ",path ,";",cmd))
}


###### SET THRESHOLDS #######
sa_percent <- 5
burnin <- 10



# combine and read logs for original FBD run
table1 <- read.table(paste0(fbdOriginalPath, "canids_rho_zero.log"), sep = "\t", header = TRUE)


# combine and read logs for original FBD run with updated dates
table2 <- read.table(paste0(fbdUpdatedPath, "canidae_rho_zero.log"), sep = "\t", header = TRUE)


# combine and read logs for original sRanges run (also with updated dates)
table3 <- read.table(paste0(sRangesBothPath, "canidae_rho.combined.log"), sep = "\t", header = TRUE)


# combine and read logs for original sRanges run (also with updated dates)
table4 <- read.table(paste0(sRangesPath, "canidae_rho.combined.log"), sep = "\t", header = TRUE)



if ("TreeLength" %notin% colnames(table1)){
  table1$TreeLength <- as.numeric(unlist(read.csv(paste0(fbdOriginalPath,"canids_rho_zero.length.log"),
                                                  sep="\t")["tree.treeLength"]))
}
if ("TreeLength" %notin% colnames(table2)){
  table2$TreeLength <- as.numeric(unlist(read.csv(paste0(fbdUpdatedPath,"canidae_rho_zero.length.log"),
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

table1<- table1[-1:-round(nrow(table1)*(burnin/100)),]
table2 <- table2[-1:-round(nrow(table2)*(burnin/100)),]
table3 <- table3[-1:-round(nrow(table3)*(burnin/100)),]
table4 <- table4[-1:-round(nrow(table4)*(burnin/100)),]


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
matching_columns <- Reduce(intersect, list(colnames(table1), colnames(table2), colnames(table3), colnames(table4)))

# combined_df_all$Table <- as.factor(combined_df_all$Table)
combined_df$model <- as.factor(combined_df$model)

# hpd <- data.frame(name=matching_columns_all[6:26], hpd.low.fbd=numeric(21), hpd.high.fbd=numeric(21), hpd.width.fbd=numeric(21),
#                   hpd.low.fbdUp=numeric(21), hpd.high.fbdUp=numeric(21), hpd.width.fbdUp=numeric(21),
#                   hpd.low.srfbd=numeric(21), hpd.high.srfbd=numeric(21), hpd.width.srfbd=numeric(21))

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

# library(ggpubr)
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

# Combine tree log files for FBD run and 
# read the table of from SampledAncestorTreeAnalyser (from SA package)
# combine_logs(fbdOriginalPath, "Mkv_morph\\+dna", "trees", burnin)
# table1 <- getSATable(paste0(fbdOriginalPath, "Mkv_morph+dna.combined.trees"),
#                      paste0(fbdOriginalPath, "combined_SA_analysis.txt"))

# Combine tree log files for FBD run with updated dates and 
# read the table of from SampledAncestorTreeAnalyser (from SA package)
# combine_logs(fbdUpdatedPath, "canidae_rho_zero", "trees", burnin, 0)
table2 <- getSATable(paste0(fbdUpdatedPath, "canidae_rho_zero.combined.trees"),
                     paste0(fbdUpdatedPath, "combined_SA_analysis.txt"), analyse=F)

# Combine tree log files for sRAnges run with updated dates and 
# read the table of from SampledAncestorTreeAnalyser (from SA package)
# combine_logs(sRangesPath, "canidae_rho", "trees", burnin, "500000")
# table3 <- getSATable(paste0(sRangesPath, "canidae_rho.combined.trees"),
#                      paste0(sRangesPath, "combined_SA_analysis.txt"), analyse=F)

table3 <- read.csv(paste0(figure_path, "morph_at_first_direct_ancestor_probs.csv"))
table3 <- table3[, c("species", "labels")]
table3$species <- str_replace_all(table3$species, "_", " ")
table3$labels <- table3$labels*100
colnames(table3)<-c("SA", "Percent")

# to get the SA analysis log for sRnages, you can run the SA analysis script of SA package from App Launcher in Beast2. 
# Then for all species that have a range associated with it, the sample representing the start of this range should be removed form the log 
# (it will also have an assigned percentage if 100%, but it's not a real sampled ancestors)
# id <- which(table3$Percent==100)
# table3 <- table3[-id,]

# table1$SA <- gsub("_", " ", table1$SA)
table2$SA <- gsub("_", " ", table2$SA)
table3$SA <- gsub("_", " ", table3$SA)

# table1 <- table1[which(table1$Percent>5),]
table2 <- table2[which(table2$Percent>5),]
table3 <- table3[which(table3$Percent>5),]


# Make the plot only for updated dates SA analysis
p<-ggplot(table2, aes(x = SA, y = Percent*0.01)) +
  geom_point(aes(color = "FBD", shape="FBD"), size = 4, alpha=0.8) + 
  # geom_point(data=table1,
  #            aes(x = SA, y = Percent*0.01, color = "FBD", shape="FBD"),
  #            size = 4, alpha=0.8) +
  geom_point(data=table3,
             aes(x = SA, y = Percent*0.01, color = "SRFBD", shape="SRFBD"), 
             size = 4, alpha=0.8) +
  # geom_point(data=ranges_full_s, 
  #            aes(x = range, y = median_s/median_s_prior, color = "Median Start"), 
  #            size = 3, alpha=0.7)+
  # geom_point(data=ranges_full_e, 
  #            aes(x = range, y = median_e/median_e_prior, color = "Median Last"), 
  #            size = 3, alpha=0.7) +
  ylab("Probability of being SA") + 
  xlab("Species") + 
  scale_color_manual(values = c("FBD" = "#FF7F50", "SRFBD"="#43AA8B")) +
  scale_shape_manual(values = c("FBD" = 16, "SRFBD"=18)) +
  theme_minimal()  + labs(color="", shape="")+
  theme(text=element_text(size=21), legend.position = c(0.80, 0.1), 
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"), 
        legend.title = element_blank())+coord_flip()+
  scale_x_discrete(limits=rev(c(levels(reorder(table3$SA,table3$Percent, decreasing = T)), 
                                levels(reorder(table2$SA[table2$SA %notin% table3$SA],table2$Percent[table2$SA %notin% table3$SA], decreasing = T)))))

ggsave(paste0(figure_path,"canid_SA_prob_comparison_updated_dates_full", figure_extension), p, width = 8, height=18)

name<-unique(c(table2$SA[which(table2$Percent>30)], table3$SA[which(table3$Percent>30)]))
table2_30 <- table2[which(table2$SA %in% name),]
table3_30 <- table3[which(table3$SA %in% name),]

p<-ggplot(table2_30, aes(x = SA, y = Percent*0.01)) +
  geom_point(aes(color = "FBD", shape="FBD"), size = 4, alpha=0.8) + 
  # geom_point(data=table1,
  #            aes(x = SA, y = Percent*0.01, color = "FBD", shape="FBD"),
  #            size = 4, alpha=0.8) +
  geom_point(data=table3_30,
             aes(x = SA, y = Percent*0.01, color = "SRFBD", shape="SRFBD"), 
             size = 4, alpha=0.8) +
  # geom_point(data=ranges_full_s, 
  #            aes(x = range, y = median_s/median_s_prior, color = "Median Start"), 
  #            size = 3, alpha=0.7)+
  # geom_point(data=ranges_full_e, 
  #            aes(x = range, y = median_e/median_e_prior, color = "Median Last"), 
  #            size = 3, alpha=0.7) +
  ylab("Probability of being SA") + 
  xlab("Species") + 
  scale_color_manual(values = c("FBD" = "#FF7F50", "SRFBD"="#43AA8B")) +
  scale_shape_manual(values = c("FBD" = 16, "SRFBD"=18)) +
  theme_minimal()  + labs(color="", shape="")+
  theme(text=element_text(size=21), legend.position = c(0.80, 0.08), 
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"), 
        legend.title = element_blank())+
  scale_y_continuous(breaks=c(0.2, 0.4, 0.6, 0.8))+
  scale_x_discrete(limits=rev(c(levels(reorder(table3_30$SA,table3_30$Percent, decreasing = T)), 
                                levels(reorder(table2_30$SA[table2_30$SA %notin% table3_30$SA],
                                               table2_30$Percent[table2_30$SA %notin% table3_30$SA], decreasing = T)))))+
  coord_flip()

ggsave(paste0(figure_path,"canid_SA_prob_comparison_updated_dates", figure_extension), p, width = 8, height=15)

