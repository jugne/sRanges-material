get_index <- function(fl, bs, tree = TRUE) {
  bit <- paste0(".*/", bs)
  if (tree) {
    ans <- gsub(".trees", "", gsub(bit, "", fl))
  } else{
    ans <- gsub(".log", "", gsub(bit, "", fl))
  }
  return(as.numeric(ans))
}

burn_and_thin <- function(dt, b, k) {
  N <- nrow(dt)
  nb <- round(b * N) + 1
  if (N - nb < k)
    stop(paste0("Can't keep ", k, " samples after ", b, " burnin"))
  burned <- dt[nb:N,]
  Nstar <- nrow(burned)
  thinned <- burned[round(seq(1, Nstar, length.out = k)),]
  return(thinned)
}
organise_posterior_draws <- function(post_draws, i, bnin, nkeep) {
  out <- data.frame(origin = post_draws$origin,
                    diversificationRate = post_draws$diversificationRate,
                    turnover = post_draws$turnover,
                    samplingProportion = post_draws$samplingProportion,
                    samplingAtPresentProb = post_draws$samplingAtPresentProb)
  row.names(out) <- NULL
  out <- data.frame(out, data_set = paste0("data_set_", i))
  out <- burn_and_thin(dt = out, b = bnin, k = nkeep)
  return(out)
}

process_postfile <- function(fl, index,
                             burnin, n_keep,
                             tree = TRUE) {
  if (tree) {
    draws <- ape::read.nexus(fl)
    ans <- organise_posterior_draws_tree(
      post_draws = draws,
      i = index,
      bnin = burnin,
      nkeep = n_keep,
      distances = do_distances
    )
  } else{
    draws <- read.table(fl, header = TRUE)
    ans <- organise_posterior_draws(
      post_draws = draws,
      i = index,
      bnin = burnin,
      nkeep = n_keep
    )
  }
  return(ans)
}

folder <- c("estimate_all_params_ext_dna", "estimate_all_params_full_dna")
thisDir <- getwd()
burn <- c(0.1, 0.1)
remove <- c(182, 182)
n_sims <- c(199, 199)

for (m in 1:2){
  setwd(paste0(thisDir,"/",folder[m]))
  f.burnin <- burn[m]
  n.keep <- 500
  
  logs.post.files <- system("ls run_*/inf/*.log", intern = TRUE)
  logs.inds <- sapply(logs.post.files, get_index,
                      bs = "sRanges.", FALSE)
  
  sorted.logs.post.files <- names(sort(logs.inds))
  # sorted.logs.post.files <- sapply(which(post_ess>30), function(x) paste0("run_", x, "/inf/sRanges.",x,".log"))
  
  prior <- read.csv("true_rates.csv")
  if (!is.na(remove[m])){
    sorted.logs.post.files <- sorted.logs.post.files[-remove[m]]
    prior <- prior[-remove[m],]
  }
  
  
  names(prior)[2:5] <- c("diversificationRate", "turnover"  ,  "samplingProportion"  ,  "samplingAtPresentProb")
  ranks <- list()
  for (i in 1:length(sorted.logs.post.files)){
    continuous <- process_postfile(
      fl = sorted.logs.post.files[i],
      index = get_index(fl = sorted.logs.post.files[i],
                        bs = "sRanges.",
                        tree = FALSE),
      tree = FALSE,
      burnin = f.burnin,
      n_keep = n.keep
    )
    
    if (!is.na(remove[m])){
      names(prior)[2:5] <- names(continuous)[2:5]
    }
    
    p <-prior[i,]
    
    
    ranks<-append(ranks, SBC:::calculate_ranks_draws_matrix(variables = p[,names(continuous)[1:5]],
                                                            dm = posterior::as_draws_matrix(continuous[,1:5])))
    
  }
  
  dt <- data.frame(variable=names(ranks), rank=unlist(ranks), sim_id=rep(1:n_sims[m], each=5))
  ecdfFormatAndSave <- function(v, filePath){
    p<- SBC::plot_ecdf(subset(dt, variable==v), max_rank=n.keep)
    p+ggplot2::theme_bw()+ggplot2::xlab("Normalised rank")+
      ggplot2::ylab("Cumulative Probability")+
      ggplot2::theme(legend.position="none",
                                         text = ggplot2::element_text(size = 30),
                                         strip.text = ggplot2::element_blank(),
                                         strip.background = ggplot2::element_blank())
    ggplot2::ggsave(filePath, width=10, height=10)
  }
  
  vars<-c("diversificationRate", "turnover"  ,  "samplingProportion"  ,  "samplingAtPresentProb", "origin")
  for (v in vars){
    ecdfFormatAndSave(v, paste0("figures/", v,"_ecdf.pdf"))
  }
}


#p<-SBC::plot_rank_hist(subset(dt, variable=="origin"), max_rank=n.keep)
#p<- SBC::plot_ecdf(subset(dt, variable=="samplingProportion"), max_rank=500)
#p+ggplot2::theme_bw()+ggplot2::theme(legend.position="none",strip.text = element_blank(), text = element_text(size = 20), strip.background = element_blank())


##########################################################################
##########################################################################
##########################################################################
##########################################################################

