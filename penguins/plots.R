#######################################################
######### plots to compare distributions      #########
######### of parameters and some tree         #########
######### statistics between FBD and sRanges  #########
#######################################################


library(ggplot2)
library(dplyr)

####################################################
## Make sure the SA package is installed in BEAST ##
####################################################

##### SET PATHS ######
appLauncherPath <- "'/Applications/BEAST 2.6.7/bin/applauncher'"
logCombinerPath <- "'/Applications/BEAST 2.6.7/bin/logcombiner'"
fbdOriginalPath <- "~/Documents/Source/beast2.7/sRanges-material/penguins/Gavryushkina_et_al/original_rerun/"
fbdUpdatedPath <- "~/Documents/Source/beast2.7/sRanges-material/penguins/Gavryushkina_et_al/updated_dates_rerun/"
sRangesPath <- "~/Documents/Source/beast2.7/sRanges-material/penguins/morph_at_start/log_branch_rates/"
figure_path <- "~/Documents/Source/beast2.7/sRanges-material/penguins/figures2/"
figure_extension <- ".svg" # e.g., "pdf", "png", "svg"


combine_logs <- function(path, pattern, extesion, burnin){
  print(paste0(pattern, "_\\d.", extesion))
  logs <- list.files(path=path, pattern = paste0(pattern, "_[0-9]\\.", extesion))
  cmd <- paste(logCombinerPath, "-b", burnin)
  for (log in logs){
    cmd <- paste(cmd, "-log", log)
  }
  cmd <- paste0(cmd, " -o ", pattern, ".combined")
  cmd <- paste0(cmd, ".", extesion)
  print(cmd)
  system(paste0("cd ",path ,";",cmd))
}


###### SET THRESHOLDS #######
sa_percent <- 5
burnin <- 10



# combine and read logs for original FBD run
combine_logs(fbdOriginalPath, "Mkv_morph\\+dna", "log", burnin)
table1 <- read.table(paste0(fbdOriginalPath, "Mkv_morph+dna.combined.log"), sep = "\t", header = TRUE)

# combine and read logs for original FBD run with updated dates
combine_logs(fbdUpdatedPath, "Mkv_morph\\+dna", "log", burnin)
table2 <- read.table(paste0(fbdUpdatedPath, "Mkv_morph+dna.combined.log"), sep = "\t", header = TRUE)

# combine and read logs for original sRanges run (also with updated dates)
combine_logs(sRangesPath, "penguins_inf_morph_at_start", "log", burnin)
table3 <- read.table(paste0(sRangesPath, "penguins_inf_morph_at_start.combined.log"), sep = "\t", header = TRUE)


# Combine the tables into a single data frame
combined_df_all <- bind_rows(table1, table2, table3, .id = "Table") # all 3 runs
combined_df <- bind_rows(table2, table3, .id = "Table") # only updated dates

# Get the unique matching column names (params and statistics that are logged for all runs)
matching_columns_all <- intersect(intersect(colnames(table1), colnames(table2)), colnames(table3))
matching_columns <- intersect(colnames(table2), colnames(table3))

combined_df_all$Table <- as.factor(combined_df_all$Table)
combined_df$Table <- as.factor(combined_df$Table)

# Create parameter and tree statistics plots comparing all 3 runs
for (col in matching_columns_all[6:26]){
  
  # Create the violin plot for the matching column name
  plot <- ggplot(combined_df_all, aes_string(x = "Table", y = col, fill = "Table", color="Table")) +
    geom_violin(alpha=0.8) +
    stat_summary(fun = "median", geom = "crossbar", width = 0.25,
                 position = position_dodge(width = .25), color = "#5A5A5A", alpha=0.8) +
    labs(x = "", y = "") +
    ggtitle("") +
    theme_minimal() +
    scale_x_discrete(labels = c("FBD", "FBD, \nUpdated Dates",  "sRanges")) +
    scale_fill_manual(values = c("#B38CB4", "#FF7F50", "#43AA8B"))+
    scale_colour_manual(values = c("#B38CB4", "#FF7F50",  "#43AA8B"))+
    theme(
      legend.position = "none",
      text = element_text(size = 24),
      plot.title = element_text(size = 30, face = "bold"),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 20),
    )
  
  ggsave(paste0(figure_path,col,figure_extension), plot)
}

# Create parameter and tree statistics plots comparing only FBD and sRanges plots with updated dates
for (col in matching_columns[6:26]){

  # Create the violin plot for the matching column name
  plot <- ggplot(combined_df, aes_string(x = "Table", y = col, fill = "Table", color="Table")) +
    geom_violin(alpha=0.8) +
    stat_summary(fun = "median", geom = "crossbar", width = 0.25,
                 position = position_dodge(width = .25), color = "#5A5A5A", alpha=0.8) +
    labs(x = "", y = "") +
    ggtitle("") +
    theme_minimal() +
    scale_x_discrete(labels = c("FBD, \nUpdated Dates", "sRanges")) +
    scale_fill_manual(values = c("#FF7F50", "#43AA8B"))+
    scale_colour_manual(values = c("#FF7F50", "#43AA8B"))+
    theme(
      legend.position = "none",
      text = element_text(size = 24),
      plot.title = element_text(size = 30, face = "bold"),
      axis.title = element_text(size = 24),
      axis.text = element_text(size = 20),
    )
  
  ggsave(paste0(figure_path,col,"_updatedDates",figure_extension), plot)
}


getSATable<- function(treeFile, outFile, analise=T){
  if (analise){
    cmd <- paste(appLauncherPath,'SampledAncestorTreeAnalyser -file', treeFile, ">", outFile)
    system(cmd)
  }
  # Determine the total number of lines in the file
  # total_lines <- length(readLines(outFile))
  # Skip the first two lines and read until (total_lines - 2)
  data <- read.table(outFile, sep = "\t", header = TRUE, skip = 2)#, nrows = total_lines - 3)
  
  old_names <- c("Burnside_Palaeudyptes", "Waimanu_tuatahi")
  new_names <- c("Burnside_Palaeeudyptes", "Muriwaimanu_tuatahi")
  for (i in 1:length(old_names)){
    id<-which(data$SA==old_names[i])
    if (length(id)!=0){
      data$SA[id] <- new_names[i]
    }
  }
  data$SA <- gsub("_last", "", gsub("_first", "", data$SA))
  data$SA <- gsub("_", " ", data$SA)
  return(data)
}

# Combine tree log files for FBD run and 
# read the table of from SampledAncestorTreeAnalyser (from SA package)
combine_logs(fbdOriginalPath, "Mkv_morph\\+dna", "trees", burnin)
table1 <- getSATable(paste0(fbdOriginalPath, "Mkv_morph+dna.combined.trees"),
                     paste0(fbdOriginalPath, "combined_SA_analysis.txt"))

# Combine tree log files for FBD run with updated dates and 
# read the table of from SampledAncestorTreeAnalyser (from SA package)
combine_logs(fbdUpdatedPath, "Mkv_morph\\+dna", "trees", burnin)
table2 <- getSATable(paste0(fbdUpdatedPath, "Mkv_morph+dna.combined.trees"),
                     paste0(fbdUpdatedPath, "combined_SA_analysis.txt"))

# Combine tree log files for sRAnges run with updated dates and 
# read the table of from SampledAncestorTreeAnalyser (from SA package)
combine_logs(sRangesPath, "penguins_inf_morph_at_start", "trees", burnin)
table3 <- getSATable(paste0(sRangesPath, "penguins_inf_morph_at_start.combined.tree"),
                     paste0(sRangesPath, "combined_SA_analysis.txt"))

# to get the SA analysis log for sRnages, you can run the SA analysis script of SA package from App Launcher in Beast2. 
# Then for all species that have a range associated with it, the sample representing the start of this range should be removed form the log 
# (it will also have an assigned percentage if 100%, but it's not a real sampled ancestors)
id <- which(table3$Percent==100)
table3 <- table3[-id,]

table1$SA <- gsub("_", " ", table1$SA)
table2$SA <- gsub("_", " ", table2$SA)
table3$SA <- gsub("_", " ", table3$SA)

table1 <- table1[which(table1$Percent>5),]
table2 <- table2[which(table2$Percent>5),]
table3 <- table3[which(table3$Percent>5),]


# Make the plot only for updated dates SA analysis
p<-ggplot(table2, aes(x = reorder(SA,Percent), y = Percent*0.01)) +
  geom_point(aes(color = "FBD", shape="FBD"), size = 4, alpha=0.8) + 
  # geom_point(data=table1,
  #            aes(x = SA, y = Percent*0.01, color = "FBD", shape="FBD"),
  #            size = 4, alpha=0.8) +
  geom_point(data=table3,
             aes(x = SA, y = Percent*0.01, color = "sRanges", shape="sRanges"), 
             size = 4, alpha=0.8) +
  # geom_point(data=ranges_full_s, 
  #            aes(x = range, y = median_s/median_s_prior, color = "Median Start"), 
  #            size = 3, alpha=0.7)+
  # geom_point(data=ranges_full_e, 
  #            aes(x = range, y = median_e/median_e_prior, color = "Median Last"), 
  #            size = 3, alpha=0.7) +
  ylab("Probability of being SA") + 
  xlab("Species") + 
  scale_color_manual(values = c("FBD" = "#FF7F50", "FBD" = "#B38CB4", "sRanges"="#43AA8B")) +
  scale_shape_manual(values = c("FBD" = 16, "FBD" = 16, "sRanges"=18)) +
  theme_minimal()  + labs(color="Method", shape="Method")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),text=element_text(size=21))

ggsave(paste0(figure_path,"SA_prob_comparison_updated_dates", figure_extension), p, width = 9, height=9)
