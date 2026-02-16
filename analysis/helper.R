#######################################################
######### Helper functions and common setup  #########
#######################################################

library(ggplot2) 
library(treeio)
library(ggtree)
library(ape)
library(stringr)
library(dplyr)
library(tidytree)
library(ggpubr)
library(coda)
library(deeptime)
library(ggridges)
library(cowplot)
library(forcats)
library(XML)


'%notin%' <- Negate('%in%')

#####################################################################################
## Make sure the SA (FBD) , TreeStat2 and sRanges, packages are installed in BEAST ##
#####################################################################################



# ggtree_source_dir <- "~/Documents/Source/ggtree/R"
beast_app_dir <- "'/Applications/BEAST 2.7.7/bin'"

# Define application paths
appLauncherPath <- file.path(beast_app_dir, "applauncher")
logCombinerPath <- file.path(beast_app_dir, "logcombiner")
stats_path <- getwd()

# Define figure extension
figure_extension <- ".pdf" # e.g., "pdf", "png", "svg"

###### SET THRESHOLDS #######
sa_percent <- 5
burnin <- 10


# # Source all R files in ggtree (if needed)
# if (dir.exists(ggtree_source_dir)) {
#   setwd(ggtree_source_dir)
#   file.sources <- list.files(pattern = "*.R")
#   sapply(file.sources, source, .GlobalEnv)
# }

# Set back to base_dir as working directory
# setwd(base_dir)


##### HELPER FUNCTIONS #####


combine_and_get_length <- function(path, pattern, burnin, model="fbd", resample=NULL, force=F){
  if (force || !file.exists(paste0(path,"/", pattern,".combined.log"))){
    combine_logs(path, pattern, "log", burnin, resample)
  }
  
  if (force || !file.exists(paste0(path,"/", pattern,".combined.trees"))){
    combine_logs(path, pattern, "trees", burnin, resample)
    get_tree_length(path, pattern)
  }
  
  if (model=="srfbd"){
    if (force || !file.exists(paste0(path,"/", pattern,".combined.speciation.log"))){
      combine_logs(path, pattern, "speciation.log", burnin, resample)
    }
  }
  
}


combine_logs <- function(path, pattern, extesion, burnin, resample=NULL){
  pattern <- str_replace(pattern, "\\+", "\\\\+")
  logs <- list.files(path=path, pattern = paste0(pattern, "_[0-9]*\\.", extesion))
  cmd <- paste(logCombinerPath, "-b", burnin)
  for (log in logs){
    cmd <- paste(cmd, "-log", log)
  }
  cmd <- paste0(cmd, " -o ", pattern, ".combined")
  cmd <- paste0(cmd, ".", extesion)
  if (!is.null(resample)){
    cmd <- paste0(cmd, " -resample ", resample)
  }
  system(paste0("cd ",path ,";",cmd))
}

get_tree_length <- function(path, pattern){
  cmd <- paste(appLauncherPath, 'TreeStat2', 
               paste0(stats_path,'/stats.txt'), 
               paste0(path,"/", pattern,".combined.trees"))
  system(paste0("cd ",path ,";",cmd))
}


relog_tip_ages <-function(path, pattern){
  if (!file.exists(paste0(path,"/",pattern,".tipAge.log"))){
    cmd <- paste0(beast_app_dir,"/beast relog_tip_age.xml")
    system(paste0("cd ",path,"/tipAgeRelog",";",cmd))
  }
}


addToHPD <- function(dataFrame, table, colName, i){
  h <- HPDinterval(as.mcmc(table[,which(colnames(table)==colName)]))
  dataFrame[i, 2] <- median(table[,which(colnames(table)==colName)])
  dataFrame[i, 3] <- h[1]
  dataFrame[i, 4] <- h[2]
  dataFrame[i, 5] <- h[2]-h[1]
  return(dataFrame)
}

getSATable <- function(path, treeFileName, outFile="combined_SA_analysis.txt",
                       force=F, old_names=NULL, new_names=NULL){
  outFile = paste0(path, "/", outFile)
  if (force || !file.exists(outFile)){
    treeFile = paste0(path, "/", treeFileName)
    cmd <- paste(appLauncherPath,'SampledAncestorTreeAnalyser -file', treeFile, ">", outFile)
    system(cmd)
  }
  # Determine the total number of lines in the file
  # total_lines <- length(readLines(outFile))
  # Skip the first two lines and read until (total_lines - 2)
  data <- read.table(outFile, sep = "\t", header = TRUE, skip = 2)#, nrows = total_lines - 3)
  
  if (!is.null(old_names) & !is.null(new_names)){
    for (i in 1:length(old_names)){
      id<-which(data$SA==old_names[i])
      if (length(id)!=0){
        data$SA[id] <- new_names[i]
      }
    }
  }
  
  # Clean up species names
  data$SA <- gsub("_last", "", gsub("_first", "", data$SA))
  data$SA <- gsub("_", " ", data$SA)
  n <- names(which(table(data$SA)==2))
  for (nn in n){
    ids <- which(data$SA==nn)
    data <- data[-ids[1],]
  }
  return(data)
}

getSATableSRFBD <- function(path, pattern, figPrefix, force=F){
  if (force || !file.exists(paste0(figure_path, figPrefix,"_speciation.csv"))){
    log <-  read.table(paste0(path,"/",pattern,".combined.speciation.log"),
                       header=TRUE, sep="\t")
    
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
    data_probs <- data.frame(species=names(nn_from), probabilities=nn_from, labels=round(nn_from, 3))
    write.csv(data_probs,
              paste0(figure_path, figPrefix, "_direct_ancestor_probs.csv"))
    write.csv(data.frame(median=median(nn), mean=mean(nn),low=HPDinterval(as.mcmc(nn))[1], high=HPDinterval(as.mcmc(nn))[2]),
              paste0(figure_path, figPrefix, "_direct_ancestor_hpd.csv"))
    
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
    write.csv(d_20, paste0(figure_path, figPrefix,"_speciation.csv"))
  }
  
  return(data_probs)
}

plotPosteriorCladeAge<- function(species, cladeComparePath){
  for (file_name in c("clade_comp_morph_first", "clade_comp_morph_both")){
    # Comaparison of the clades obtained by SA and sRanges 
    # The next log file is obtained by running a SRangesAndSACladeSetComparator tool on trees obtained by Sampled ancestors
    # and sRanges packages. For all clades that were found, the posterior probability of the clade, its height and HPD interval bounds for the height are recorded.
    
    log <-  read.table(paste0(cladeComparePath,"/",file_name,".txt"),
                       header=FALSE, sep=" ")
    colnames(log) <- c("Species", "PosteriorProbFBD", "PosteriorProbSRanges",
                       "HeightFBD", "HeightFBD_HPD_l", "HeightFBD_HPD_h",
                       "HeightSRFBD", "HeightSRFBD_HPD_l", "HeightSRFBD_HPD_h")
    
    # Add colums for probability difference and number of species in each clade
    log$probDiff <- abs(log$PosteriorProbFBD-log$PosteriorProbSRanges)
    log$nSpecies <- count.fields(textConnection(log$Species), sep = ",")
    b<-lapply(log$Species, function(x) str_remove_all(x, "[{}]"))
    log$cladeName <- lapply(b, function(x) str_flatten(unique(sapply(unlist(strsplit(x, ',')),
                                                                     function(y) str_split(y, "_")[[1]])[1,]), "--"))
    # Take only clades with higher than 0.2 probability in at least one of the methods for comparison
    thr <- 0.2
    log_thr<-  log[which(log$PosteriorProbFBD>=thr | log$PosteriorProbSRanges>=thr),]
    log_thr <- log_thr[which(log_thr$HeightFBD !=0 & log_thr$HeightSRFBD!=0),]
    
    
    # Arrows to mark clades where mean age estimated by SRFBD is more recent than that estimated by FBD
    l <- rep("", nrow(log_thr))
    l[which(log_thr$HeightFBD-log_thr$HeightSRFBD>0)]<-"\u2194"
    l<-l[order(log_thr$probDiff)]
    
    
    # plot clade ages, estimated by SRFBD and FBD and their means
    p<-ggplot(log_thr, aes(x = HeightFBD_HPD_l, y = fct_reorder(Species, probDiff))) +
      geom_segment(aes(xend = HeightFBD_HPD_h, yend = fct_reorder(Species, probDiff),
                       colour = "fbd"), alpha=0.4,  linewidth = 1, ) +
      geom_point(aes(x = HeightFBD, shape="fbd"), size = 1) +
      geom_segment(aes(x=HeightSRFBD_HPD_l, xend = HeightSRFBD_HPD_h,
                       yend = fct_reorder(Species, probDiff), colour = "srfbd"),
                   linewidth = 1, alpha=0.4) +
      geom_point(aes(x = HeightSRFBD, shape="srfbd"), size = 1) +
      scale_color_manual(values = c(fbd = "blue", srfbd = "red"),
                         labels = c(fbd = "FBD", srfbd = "SRFBD"),
                         name = "")+
      scale_shape_manual(values = c(fbd = 4, srfbd = 3),
                         labels = c(fbd = "FBD", srfbd = "SRFBD"),
                         name = "")+
      theme_bw() + xlab("Age (Ma)") + ylab("Clades") +
      theme(plot.margin=unit(c(0.2, -0.2, 1, 0.2),"cm"),
            axis.ticks.y = element_blank(), legend.direction="horizontal",
            legend.position = c(0.8, -0.15))  + scale_y_discrete(
              name = "",
              labels = l, position="right"
            )
    
    # Plot posterior clade probabilities
    p2 <- ggplot(log_thr, aes(x = 0, y = fct_reorder(Species, probDiff))) +
      geom_segment(aes(xend = PosteriorProbFBD, yend = fct_reorder(Species, probDiff)),
                   linewidth = 1, alpha=0.4, colour = "blue")+
      geom_segment(aes(x=0, xend = PosteriorProbSRanges, yend = fct_reorder(Species, probDiff)),
                   linewidth = 1, alpha=0.4, colour = "red")+
      theme_bw() + xlab("Posterior probability") +
      theme(legend.position = "none",axis.text.y = element_blank(),
            axis.title.y=element_blank(), axis.ticks.y = element_blank(), 
            plot.margin=unit(c(0.2, 0.2, 1, -0.2),"cm"))
    
    p3<-plot_grid(p, p2, ncol=2, rel_widths = c(3,1))
    p3
    ggsave2(paste0(figure_path,"/",species,"_",file_name,"_thr_",thr,".pdf"),
            p3, width=5, height=7, device=cairo_pdf, family="Arial Unicode MS")
    
    write.csv2(data.frame(nCladesInBoth=nrow(log_thr), nCladesFBDOlder=sum(log_thr$HeightFBD>log_thr$HeightSRFBD),
                          nCladesSRFBDhpdNarrower=sum((log_thr$HeightFBD_HPD_h-log_thr$HeightFBD_HPD_l)>(log_thr$HeightSRFBD_HPD_h-log_thr$HeightSRFBD_HPD_l))),
               paste0(figure_path,"/",species,"_",file_name,"_thr_",thr,"_nCladesFBDOlder.csv"))
    
    
    
    mm <- mean(abs(log_thr$HeightFBD-log_thr$HeightSRFBD))
    mm2 <- mean(log_thr$HeightFBD-log_thr$HeightSRFBD)
    cs <- mean(abs(log_thr$PosteriorProbFBD-log_thr$PosteriorProbSRanges))
    
    
    
    write.csv2(data.frame(mean_age_diff=mm,
                          mean_age_FBD_minus_SRFBD=mm2,
                          mean_clade_support=cs),
               paste0(figure_path,"/",species,"_",file_name,"_thr_",
                      thr,"_averages.csv"))
  
    }
}

#' Parse XML sampling dates and extract species information
#'
#' @param text XML text containing samplingDates elements
#' @return A list containing:
#'   - priors: data.frame with species, start, and end columns
#'   - ranges_single: vector of species with single samples
#'   - ranges: vector of species with paired samples (first/last)
#'   - ranges_: vector of all range names without suffixes
parse_sampling_dates <- function(text) {
  # Extract taxon names without the "@", and the lower and upper values
  sp <- gsub('@', '', regmatches(text, gregexpr('(?<=taxon="@)[^"]+', text, perl=TRUE))[[1]])
  start <- as.numeric(regmatches(text, gregexpr('(?<=upper=")[^"]+', text, perl=TRUE))[[1]])
  end <- as.numeric(regmatches(text, gregexpr('(?<=lower=")[^"]+', text, perl=TRUE))[[1]])
  
  # Remove '_first' or '_last' from taxon names
  ranges_ <- gsub('_(first|last)$', '', sp)
  
  priors <- data.frame(species=sp, start=start, end=end)
  
  all_string_listwise <- unlist(lapply(ranges_, unique))
  ranges_single <- names(which(table(all_string_listwise)==1))
  ranges <- names(which(table(all_string_listwise)>1))
  
  return(list(
    priors = priors,
    ranges_single = ranges_single,
    ranges = ranges,
    ranges_ = ranges_
  ))
}

#' Create species data frame for CSV output
#'
#' @param ranges_single Vector of single-sampled species
#' @param ranges Vector of paired-sampled species
#' @param sp Full species names with suffixes
#' @param start Start dates
#' @param end End dates
#' @return data.frame with species range information
create_species_dataframe <- function(ranges_single, ranges, sp, start, end, ranges_) {
  id_first <- which(ranges_ %in% ranges_single)
  sp_first <- data.frame(
    "Species" = ranges_[id_first],
    "First start" = start[id_first],
    "First end" = end[id_first],
    "Last start" = rep("-", length(id_first)),
    "Last end" = rep("-", length(id_first))
  )
  
  sp_both <- data.frame(
    "Species" = c(),
    "First start" = c(),
    "First end" = c(),
    "Last start" = c(),
    "Last end" = c()
  )
  
  for (r in ranges) {
    id1 <- which(sp == paste0(r, "_first"))
    id2 <- which(sp == paste0(r, "_last"))
    sp_both <- rbind(
      sp_both,
      data.frame(
        "Species" = r,
        "First start" = start[id1],
        "First end" = end[id1],
        "Last start" = start[id2],
        "Last end" = end[id2]
      )
    )
  }
  
  return(rbind(sp_first, sp_both))
}


#' Calculate IQR statistics for species ranges
#'
#' @param log Log data frame with posterior samples
#' @param priors Prior data frame with species, start, end
#' @param ranges_single Vector of single-sampled species
#' @param ranges Vector of paired-sampled species
#' @return List containing ranges_full_s, ranges_full_e, ranges_full, ranges_stat
get_sRange_full <- function(log, priors, ranges_single, ranges) {
  start_r <- c()
  start_r_low <- c()
  start_r_high <- c()
  start_r_min <- c()
  start_r_max <- c()
  start_r_iqr <- c()
  start_r_iqr_prior <- c()
  end_r <- c()
  end_r_low <- c()
  end_r_high <- c()
  end_r_min <- c()
  end_r_max <- c()
  end_r_iqr <- c()
  end_r_iqr_prior <- c()
  
  distr <- c()
  endpoint <- c()
  range <- c()
  range_s <- c()
  range_e <- c()
  median_s <- c()
  median_s_prior <- c()
  median_e <- c()
  median_e_prior <- c()
  
  # Process single-sampled species
  for (r in ranges_single) {
    if ((max(log[, paste0(r, "_first")]) - min(log[, paste0(r, "_first")])) > 0.1) {
      start_l <- length(log[, paste0(r, "_first")])
      distr <- c(distr, log[, paste0(r, "_first")])
      endpoint <- c(endpoint, rep("start", start_l))
      range <- c(range, rep(r, start_l))
      
      id <- which(priors$species == paste0(r, "_first"))
      un <- runif(5000, min=priors$end[id], max=priors$start[id])
      start_r_iqr <- c(start_r_iqr, IQR(log[, paste0(r, "_first")]))
      start_r_iqr_prior <- c(start_r_iqr_prior, IQR(un))
      
      range_s <- c(range_s, r)
      median_s <- c(median_s, median(log[, paste0(r, "_first")]))
      median_s_prior <- c(median_s_prior, median(un))
    }
    
    start_r <- c(start_r, median(log[, paste0(r, "_first")]))
    end_r <- c(end_r, NA)
    hpd_f <- HPDinterval(as.mcmc(log[, paste0(r, "_first")]))
    start_r_low <- c(start_r_low, hpd_f[1])
    start_r_high <- c(start_r_high, hpd_f[2])
    start_r_min <- c(start_r_min, min(log[, paste0(r, "_first")]))
    start_r_max <- c(start_r_max, max(log[, paste0(r, "_first")]))
    
    end_r_low <- c(end_r_low, NA)
    end_r_high <- c(end_r_high, NA)
    end_r_min <- c(end_r_min, NA)
    end_r_max <- c(end_r_max, NA)
  }
  
  # Process paired-sampled species
  for (r in ranges) {
    if ((max(log[, paste0(r, "_first")]) - min(log[, paste0(r, "_first")])) > 0.1) {
      start_l <- length(log[, paste0(r, "_first")])
      distr <- c(distr, log[, paste0(r, "_first")])
      endpoint <- c(endpoint, rep("start", start_l))
      range <- c(range, rep(r, start_l))
      
      id <- which(priors$species == paste0(r, "_first"))
      un <- runif(5000, min=priors$end[id], max=priors$start[id])
      start_r_iqr <- c(start_r_iqr, IQR(log[, paste0(r, "_first")]))
      start_r_iqr_prior <- c(start_r_iqr_prior, IQR(un))
      range_s <- c(range_s, r)
      median_s <- c(median_s, median(log[, paste0(r, "_first")]))
      median_s_prior <- c(median_s_prior, median(un))
    }
    
    if ((max(log[, paste0(r, "_last")]) - min(log[, paste0(r, "_last")])) > 0.1) {
      end_l <- length(log[, paste0(r, "_last")])
      distr <- c(distr, log[, paste0(r, "_last")])
      endpoint <- c(endpoint, rep("end", end_l))
      range <- c(range, rep(r, end_l))
      
      id <- which(priors$species == paste0(r, "_last"))
      un <- runif(5000, min=priors$end[id], max=priors$start[id])
      end_r_iqr <- c(end_r_iqr, IQR(log[, paste0(r, "_last")]))
      end_r_iqr_prior <- c(end_r_iqr_prior, IQR(un))
      range_e <- c(range_e, r)
      median_e <- c(median_e, median(log[, paste0(r, "_last")]))
      median_e_prior <- c(median_e_prior, median(un))
    }
    
    start_r <- c(start_r, median(log[, paste0(r, "_first")]))
    end_r <- c(end_r, median(log[, paste0(r, "_last")]))
    
    hpd_f <- HPDinterval(as.mcmc(log[, paste0(r, "_first")]))
    start_r_low <- c(start_r_low, hpd_f[1])
    start_r_high <- c(start_r_high, hpd_f[2])
    start_r_min <- c(start_r_min, min(log[, paste0(r, "_first")]))
    start_r_max <- c(start_r_max, max(log[, paste0(r, "_first")]))
    
    hpd_l <- HPDinterval(as.mcmc(log[, paste0(r, "_last")]))
    end_r_low <- c(end_r_low, hpd_l[1])
    end_r_high <- c(end_r_high, hpd_l[2])
    end_r_min <- c(end_r_min, min(log[, paste0(r, "_last")]))
    end_r_max <- c(end_r_max, max(log[, paste0(r, "_last")]))
  }
  
  ranges_full_s <- data.frame(
    range = gsub("_", " ", range_s),
    start_r_iqr = start_r_iqr,
    start_r_iqr_prior = start_r_iqr_prior,
    median_s_prior = median_s_prior,
    median_s = median_s
  )
  
  ranges_full_e <- data.frame(
    range = gsub("_", " ", range_e),
    end_r_iqr = end_r_iqr,
    end_r_iqr_prior = end_r_iqr_prior,
    median_e_prior = median_e_prior,
    median_e = median_e
  )
  
  ranges_full <- data.frame(distr=distr, endpoint=endpoint, range=range)
  
  length(end_r) <- length(start_r)
  ranges_stat <- data.frame(
    range = c(ranges_single, ranges),
    median_start = start_r,
    median_end = end_r,
    start_r_min = start_r_min,
    start_r_low = start_r_low,
    start_r_max = start_r_max,
    start_r_high = start_r_high,
    end_r_min = end_r_min,
    end_r_low = end_r_low,
    end_r_max = end_r_max,
    end_r_high = end_r_high
  )
  
  return(list(ranges_full_s, ranges_full_e, ranges_full, ranges_stat))
}


#' Create IQR ratio plot
#'
#' @param ranges_full_s_comb Combined start ranges dataframe
#' @param ranges_full_e_comb Combined end ranges dataframe
#' @param y_scale_params Optional list with breaks, limits, labels for y-axis
#' @param show_legend Logical, whether to show legend
#' @param text_size Text size for theme
#' @return ggplot object
create_iqr_plot <- function(ranges_full_s_comb, ranges_full_e_comb,
                            y_scale_params = NULL,
                            show_legend = FALSE,
                            text_size = 21) {
  p <- ggplot(ranges_full_s_comb, aes(x = range, y = start_r_iqr_prior/start_r_iqr)) +
    geom_point(aes(color = "First", shape=model), size = 3, alpha=0.7) +
    geom_point(data=ranges_full_e_comb,
               aes(x = range, y = end_r_iqr_prior/end_r_iqr, color = "Last", shape=model),
               size = 3, alpha=0.7) +
    geom_hline(yintercept = 1) +
    ylab("IQR ratio (prior/posterior)") +
    xlab("Species") +
    scale_color_manual(values = c("First" = "#FF7F50", "Last" = "#43AA8B")) +
    theme_minimal() +
    labs(color="Sample", shape="Morph. data attached to") +
    theme(
      legend.background = element_rect(fill="white", linewidth=0.5, linetype="solid"),
      legend.justification = c("right", "bottom"),
      legend.box = "horizontal",
      text = element_text(size=text_size),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
    ) +
    guides(color = guide_legend(override.aes = list(shape=15, size=5)))
  
  if (!show_legend) {
    p <- p + theme(legend.position = "None")
  }
  
  if (!is.null(y_scale_params)) {
    p <- p + scale_y_continuous(
      breaks = y_scale_params$breaks,
      limits = y_scale_params$limits,
      labels = y_scale_params$labels
    )
  }
  
  return(p)
}


#' Process sRange analysis workflow
#'
#' @param log_both Log data from "both" model
#' @param log_first Log data from "first" model
#' @param priors Prior dataframe
#' @param ranges_single Single-sampled species
#' @param ranges Paired-sampled species
#' @return List with ranges_full_s_comb and ranges_full_e_comb
process_srange_workflow <- function(log_both, log_first, priors, ranges_single, ranges) {
  all_both <- get_sRange_full(log_both, priors, ranges_single, ranges)
  ranges_full_s <- all_both[[1]]
  ranges_full_e <- all_both[[2]]
  
  all_first <- get_sRange_full(log_first, priors, ranges_single, ranges)
  ranges_full_s1 <- all_first[[1]]
  ranges_full_e1 <- all_first[[2]]
  
  ranges_full_s1$model <- "First"
  ranges_full_e1$model <- "First"
  
  ranges_full_s$model <- "Both"
  ranges_full_e$model <- "Both"
  
  ranges_full_s_comb <- rbind(ranges_full_s1, ranges_full_s)
  ranges_full_e_comb <- rbind(ranges_full_e1, ranges_full_e)
  
  return(list(
    ranges_full_s_comb = ranges_full_s_comb,
    ranges_full_e_comb = ranges_full_e_comb
  ))
}