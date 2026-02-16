#######################################################
######### Parameter and tree statistics      #########
######### comparison plots (FBD vs SRFBD)    #########
#######################################################

rm(list = ls())
source("helper.R")

# Read logs for all runs
table1 <- read.table(paste0(fbdOriginalPath, "canids_rho_zero.log"), sep = "\t", header = TRUE)
table2 <- read.table(paste0(fbdUpdatedPath, "canidae_rho_zero.log"), sep = "\t", header = TRUE)
table3 <- read.table(paste0(sRangesBothPath, "canidae_rho.combined.log"), sep = "\t", header = TRUE)
table4 <- read.table(paste0(sRangesPath, "canidae_rho.combined.log"), sep = "\t", header = TRUE)

# Add TreeLength if not present
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

# Apply burnin
table1 <- table1[-1:-round(nrow(table1)*(burnin/100)),]
table2 <- table2[-1:-round(nrow(table2)*(burnin/100)),]
table3 <- table3[-1:-round(nrow(table3)*(burnin/100)),]
table4 <- table4[-1:-round(nrow(table4)*(burnin/100)),]

# Clean column names
colnames(table1) <- str_remove_all(colnames(table1), "SABD")
colnames(table2) <- str_remove_all(colnames(table2), "SABD")
colnames(table3) <- str_remove_all(colnames(table3), "SRFBD")
colnames(table4) <- str_remove_all(colnames(table4), "SRFBD")

# Add model labels
table1$model <- "FBD, orig."
table2$model <- "FBD, upd."
table3$model <- "SRFBD, both"
table4$model <- "SRFBD, first"

# Combine tables
combined_df <- bind_rows(table1, table2, table3, table4)

# Get matching column names across all runs
matching_columns <- Reduce(intersect, list(colnames(table1), colnames(table2), colnames(table3), colnames(table4)))

combined_df$model <- as.factor(combined_df$model)

# Initialize HPD data frames
hpd_fbd_orig <- data.frame(name=matching_columns[6:22], median=numeric(17),
                  hpd.low=numeric(17), hpd.high=numeric(17), hpd.width=numeric(17))
hpd_fbd_upd <- hpd_fbd_orig
hpd_srfbd_first <- hpd_fbd_orig
hpd_srfbd_both <- hpd_fbd_orig
hpd_frames <- list(hpd_fbd_orig, hpd_fbd_upd, hpd_srfbd_first, hpd_srfbd_both)
tables <- list(table1, table2, table3, table4)
names <- c("hpd_fbd_orig","hpd_fbd_upd", "hpd_srfbd_first", "hpd_srfbd_both")

i <- 1

# Create parameter and tree statistics plots
for (col in matching_columns[6:22]){

  # Create the violin plot
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

  # Calculate HPD intervals
  for (j in 1:length(hpd_frames)){
    hpd_frames[[j]] <- addToHPD(hpd_frames[[j]], tables[[j]], col, i)
  }
  i <- i+1
}

# Save HPD tables
for (j in 1:length(hpd_frames)){
  hpd <- hpd_frames[[j]]
  hpd[,2:ncol(hpd)] <- round(hpd[,2:ncol(hpd)], 3)
  write.csv(hpd, paste0(figure_path, "canid_", names[j],".csv"))
}

print("Parameter comparison plots completed successfully!")
