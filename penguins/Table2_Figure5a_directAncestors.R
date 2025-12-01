rm(list = ls())
gc()

library(coda)
library(stringr)

'%notin%' <- Negate('%in%')

source("penguin_paths.R")
figure_path <- paste0(table2_figure5a_dir,"/")


getSATable<- function(treeFile, outFile, analize=T){
  if (analize){
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
  n <- names(which(table(data$SA)==2))
  for (nn in n){
    ids <- which(data$SA==nn) 
    data <- data[-ids[1],]
  }
  return(data)
}



for (morph in c("start", "both")){
  # setwd(paste0("morph_at_",morph,"/"))
  log <-  read.table(paste0(base_dir,"/morph_at_",morph,"/penguins_inf_morph_at_",morph,".combined.speciation.log"),
                     header=TRUE, sep="\t")
  
  # burnIn <- 0.1
  # log <- tail(log, -round(nrow(log)*burnIn)) # take out 10% burnin
  # log<-log[,1:(ncol(log)-1)]
  columnsWithValue <- colnames(log)[colSums(log != "0.0,0.0") != 0 & !is.na(colSums(log != "0.0,0.0"))]
  # we don't need "Sample" column
  columnsWithValue <- columnsWithValue[2:length(columnsWithValue)]
  log <- log[, columnsWithValue]
  
  n <- length(columnsWithValue)
  d <- data.frame(from=character(n), to=character(n), prob=numeric(n), length_median=numeric(n),
                  length_hpd_l=numeric(n), length_hpd_h=numeric(n), nodes_median=numeric(n),
                  nodes_hpd_l=numeric(n), nodes_hpd_h=numeric(n))
  # Split each string by "."
  split_strings <- strsplit(columnsWithValue, "\\.")
  # Extract the first and second elements
  d$from <- sapply(split_strings, function(x) x[1])
  d$to <- sapply(split_strings, function(x) x[2])
  
  tmp <- unique(d$from)
  nn_from <- numeric(length(tmp))
  names(nn_from) <- tmp
  rm(tmp)
  
  nn <- numeric(nrow(log))
  for (j in 1:nrow(log)){
    tmp <- unique(sapply(strsplit(names(log[j,log[j,]!= "0.0,0.0"]), "\\."),
                         function(x) x[1]))
    nn_from[tmp] = nn_from[tmp] + 1
    nn[j] <- length(tmp)
  }
  
  nn_from <- nn_from/nrow(log)
  write.csv(data.frame(species=names(nn_from), probabilities=nn_from, labels=round(nn_from, 3)),
            paste0(figure_path,"morph_at_",morph,"_direct_ancestor_probs.csv"))
  write.csv(data.frame(median=median(nn), mean=mean(nn),low=HPDinterval(as.mcmc(nn))[1], high=HPDinterval(as.mcmc(nn))[2]),
            paste0(figure_path,"morph_at_",morph,"_direct_ancestor_hpd.csv"))
  
  for (i in 1:n){
    split_strings <- strsplit(log[,columnsWithValue[i]], ",")
    length <- sapply(split_strings, function(x) as.double(x[1]))
    nodes <- sapply(split_strings, function(x) as.double(x[2]))
    
    d$prob[i] <- length(length)
    idx <- which(length != 0)
    length <- length[idx] # no length means no speciation between the two species recorded.
    nodes <- nodes[idx] # 0 intermediate nodes are legit
    d$prob[i] <- length(length)/d$prob[i]
    
    if (length(idx)>1){
      d$length_median[i] <- median(length)
      length_hpd <- HPDinterval(as.mcmc(length))
      d$length_hpd_l[i] <- length_hpd[1]
      d$length_hpd_h[i] <- length_hpd[2]
      
      d$nodes_median[i] <- median(nodes)
      nodes_hpd <- HPDinterval(as.mcmc(nodes))
      d$nodes_hpd_l[i] <- nodes_hpd[1]
      d$nodes_hpd_h[i] <- nodes_hpd[2]
    }
  }
  
  d <- d[order(d$prob,decreasing = TRUE),]
  head(d, 30)
  
  d_20 <- d[which(d$prob>0.10),]
  d_20$prob <- round(d_20$prob, 3)
  write.csv(d_20, paste0(figure_path,
                         paste0("morph_at_",morph,"_speciation.csv")))
}


an <- F
if (!file.exists(paste0(fbdOriginalPath, "/combined_SA_analysis.txt"))){
  an <- T
}
table1 <- getSATable(paste0(fbdOriginalPath, "/Mkv_morph+dna.combined.trees"),
                     paste0(fbdOriginalPath, "/combined_SA_analysis.txt"), analize=an)

an <- F
if (!file.exists(paste0(fbdUpdatedPath, "/combined_SA_analysis.txt"))){
  an <- T
}
table2 <- getSATable(paste0(fbdUpdatedPath, "/Mkv_morph+dna.combined.trees"),
                     paste0(fbdUpdatedPath, "/combined_SA_analysis.txt"), analize=an)

# Combine tree log files for sRAnges run with updated dates and 
# read the table of from SampledAncestorTreeAnalyser (from SA package)
# combine_logs(sRangesPath, "penguins_inf_morph_at_start", "trees", burnin)


table3 <- read.csv(paste0(figure_path, "morph_at_start_direct_ancestor_probs.csv"))
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

