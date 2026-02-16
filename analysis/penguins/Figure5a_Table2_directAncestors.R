rm(list = ls())
gc()

source("penguin_helper.R")
figure_path <- paste0(figure5a_table2_dir,"/")

table1 <- getSATable(fbdOriginalPath, "Mkv_morph+dna.combined.trees")
table2 <- getSATable(fbdUpdatedPath, "Mkv_morph+dna.combined.trees")


# for SRFBD with morph data at both ends, generate the files, but don't read it - won't produce the figure
getSATableSRFBD(sRangesBothPath, "penguins_inf_morph_at_both", "morph_at_both")

table3 <- getSATableSRFBD(sRangesPath, "penguins_inf_morph_at_start", "morph_at_start")
table3 <- table3[, c("species", "labels")]
table3$species <- str_replace_all(table3$species, "_", " ")
table3$labels <- table3$labels*100
colnames(table3)<-c("SA", "Percent")

table1$SA <- gsub("_", " ", table1$SA)
table2$SA <- gsub("_", " ", table2$SA)
table3$SA <- gsub("_", " ", table3$SA)

table1 <- table1[which(table1$Percent>sa_percent),]
table2 <- table2[which(table2$Percent>sa_percent),]
table3 <- table3[which(table3$Percent>sa_percent),]


# Make the plot only for updated dates SA analysis
p<-ggplot(table2, aes(x = reorder(SA,Percent), y = Percent*0.01)) +
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
  theme(text=element_text(size=21), legend.position = c(0.750, 0.08), 
        legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black"), legend.title = element_blank())+coord_flip()+
  scale_x_discrete(limits=rev(c(levels(reorder(table3$SA,table3$Percent, decreasing = T)), 
                                levels(reorder(table2$SA[table2$SA %notin% table3$SA],table2$Percent[table2$SA %notin% table3$SA], decreasing = T)))))

ggsave(paste0(figure_path,"penguin_SA_prob_comparison_updated_dates", figure_extension), p, width = 8, height=11)

