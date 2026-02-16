#######################################################
######### Sampled ancestor probability       #########
######### comparison plots (FBD vs SRFBD)    #########
#######################################################

rm(list = ls())
gc()
source("canid_helper.R")

figure_path <- paste0(figures5b_A14_table2_dir,"/")

# Read sampled ancestor data for FBD run with updated dates
table2 <- getSATable(paste0(fbdUpdatedPath, "canidae_rho_zero.combined.trees"),
                     paste0(fbdUpdatedPath, "combined_SA_analysis.txt"), analyse=F)


# for SRFBD with morph data at both ends, generate the files, but don't read it - won't produce the figure
getSATableSRFBD(sRangesBothPath, "canid_morph_at_both", "morph_at_both")

table3 <- getSATableSRFBD(sRangesPath, "canid_morph_at_start", "morph_at_start")

table3 <- table3[, c("species", "labels")]
table3$species <- str_replace_all(table3$species, "_", " ")
table3$labels <- table3$labels*100
colnames(table3) <- c("SA", "Percent")

# Clean species names
table2$SA <- gsub("_", " ", table2$SA)
table3$SA <- gsub("_", " ", table3$SA)

# Filter by SA threshold
table2 <- table2[which(table2$Percent>sa_percent),]
table3 <- table3[which(table3$Percent>sa_percent),]

# Create full comparison plot (all species with SA > threshold)
p <- ggplot(table2, aes(x = SA, y = Percent*0.01)) +
  geom_point(aes(color = "FBD", shape="FBD"), size = 4, alpha=0.8) +
  geom_point(data=table3,
             aes(x = SA, y = Percent*0.01, color = "SRFBD", shape="SRFBD"),
             size = 4, alpha=0.8) +
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

# Create focused comparison plot (species with SA > 30%)
name <- unique(c(table2$SA[which(table2$Percent>30)], table3$SA[which(table3$Percent>30)]))
table2_30 <- table2[which(table2$SA %in% name),]
table3_30 <- table3[which(table3$SA %in% name),]

p <- ggplot(table2_30, aes(x = SA, y = Percent*0.01)) +
  geom_point(aes(color = "FBD", shape="FBD"), size = 4, alpha=0.8) +
  geom_point(data=table3_30,
             aes(x = SA, y = Percent*0.01, color = "SRFBD", shape="SRFBD"),
             size = 4, alpha=0.8) +
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
