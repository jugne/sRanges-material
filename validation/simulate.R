rm(list = ls())
options(digits=16)

library(ape)
library(phytools)
library(ggplot2)
library(stringr)
library(beastio)
library(FossilSim)
library(R.utils)

'%notin%' <- Negate('%in%')
wd <- "/Users/jugne/Documents/Source/beast2.7/sRanges-material/validation/"
templates_dir <- "/Users/jugne/Documents/Source/beast2.7/sRanges-material/validation/templates/"


### Tree and fossil simulation from supplementary material of 
### https://doi.org/10.1098/rspb.2019.0685

set.seed(5647829)

### Simulation parameters

f_lambda <- function(diversification, turnover){
  return(diversification/(1-turnover))
}

f_mu <- function(turnover, lambda){
  return(turnover*lambda)
}

f_psi <- function(mu, sampling_prop){
  return(mu*sampling_prop/(1-sampling_prop))
}




div_rate = NULL
turnover = NULL
sampling_prop = NULL
sampl_extant_prob = NULL
lambda = NULL
mu = NULL
psi = NULL
redraws = 0
origin = NULL

redraw_params <- function(){
  div_rate <<- rlnorm(1, meanlog = -3.5, sdlog = 1.5)
  turnover <<- runif(1)
  sampling_prop <<- runif(1, 0.5, 1)
  sampl_extant_prob <<- runif(1, 0.8, 1.)
  
  lambda <<- f_lambda(div_rate, turnover)
  mu <<- f_mu(turnover, lambda)
  psi <<- f_psi(mu, sampling_prop)
  origin <<- runif(1, 40, 80)
}

min_age = 10
max_age = 100

ntrees = 1000
# nextant = 25

min_samples = 5
max_samples = 150

true_rates <- data.frame(div_rate=numeric(), turnover=numeric(),
                         sampling_prop=numeric(), rho=numeric(),
                         origin=numeric(), mrca=numeric(), tree=character(),
                         n_samples=numeric(), redraws = numeric())

### Simulate Trees

# simulating trees and fossils with parameters
trees = list()
fossils = list()
taxonomy = list()
samp_trees = list()
while(length(trees) < ntrees) {
  tryCatch({
  redraw_params()
    redraws = redraws+1
    tree_tmp <- withTimeout({
      TreeSim::sim.bd.age(origin, 1, lambda, mu)
    }, timeout = 0.005)

  
  if (length(tree_tmp[[1]])==1){
    next
  }
  # tree_tmp <- TreeSim::sim.bd.taxa(nextant, 1, lambda, mu, complete = T)
  mrca = max(ape::node.depth.edgelength(tree_tmp[[1]]))
  if (length(which((mrca - ape::node.depth.edgelength(tree_tmp[[1]])) < 1e-7)) < 1){
    next
  }
  # origin = mrca + tree_tmp[[1]]$root.edge
  # if(origin > max_age || origin < min_age) {
  #   next
  # }
  
  taxonomy_tmp = sim.taxonomy(tree_tmp[[1]])
  fossils_tmp = FossilSim::sim.fossils.poisson(psi, tree = tree_tmp[[1]], taxonomy = taxonomy_tmp)
  
  if(length(fossils_tmp$edge) < 4 ||
     length(fossils_tmp$edge) > 100) {
    next
  }
  
  
  if(length(tree_tmp[[1]]$tip.label) < min_samples ||
      length(tree_tmp[[1]]$tip.label) > max_samples) {
    next
  }
  

  trees = c(trees, tree_tmp)
  fossils[[length(trees)]] = fossils_tmp
  taxonomy = c(taxonomy, taxonomy_tmp)
  true_rates_tmp <- data.frame(div_rate, turnover, sampling_prop,
                               rho=sampl_extant_prob, origin, mrca, tree="",
                               n_samples=0, redraws=redraws)
  true_rates = rbind(true_rates, true_rates_tmp)
  redraws = 0
}, TimeoutException = function(ex) {
  message("Timeout. Skipping.")
  
})
  next}

for (i in 1:ntrees){
  
  beast_tree = beast.fbd.format(trees[[i]],fossils[[i]], rho=sampl_extant_prob, digits=16)
  true_rates$tree[i] = beast_tree
  tmp_tree = ape::read.tree(text=beast_tree)
  true_rates$mrca[i] = max(ape::node.depth.edgelength(tmp_tree))
  true_rates$n_samples[i] = length(tmp_tree$tip.label)
  sample_times = ape::node.depth.edgelength(tmp_tree)
  sample_times = max(sample_times) - sample_times
  sample_times_round = round(sample_times,15)
  sample_times_round[which(sample_times_round < 1e-10)] = 0
  
  taxon = list()
  taxon_extant = list()
  strat_ranges = list()
  strat_ranges_refs = list()
  j=0
  k=1
  for (k in 1:length(unique(sub("_.*", "", tmp_tree$tip.label)))){
    tip = unique(sub("_.*", "", tmp_tree$tip.label))[k]
    taxon <- append(taxon, paste0(tip, "_first"))
    if (length(which(sub("_.*", "", tmp_tree$tip.label) %in% tip))>1){
      taxon <- append(taxon, paste0(tip, "_last"))
      strat_ranges <- append(strat_ranges,
                             paste0('<stratigraphicRange id="r',
                                    j,'" spec="StratigraphicRange" firstOccurrence="@',
                                    paste0(tip, "_first"),'" lastOccurrence="@',
                                    paste0(tip, "_last"),'"/>'))
      if (sample_times_round[which(tmp_tree$tip.label==paste0(tip, "_last"))]==0.0){
        taxon_extant<- append(taxon_extant, paste0(tip, "_last"))
      }
    } else{
      strat_ranges <- append(strat_ranges,
                             paste0('<stratigraphicRange id="r',
                                    j,'" spec="StratigraphicRange" firstOccurrence="@',
                                    paste0(tip, "_first"),'" lastOccurrence="@',
                                    paste0(tip, "_first"),'"/>'))
      if (sample_times_round[which(tmp_tree$tip.label==paste0(tip, "_first"))]==0.0){
        taxon_extant<- append(taxon_extant, paste0(tip, "_first"))
      }
    }
    strat_ranges_refs <- append(strat_ranges_refs, paste0('<stratigraphicRange idref="r',j,'"/>'))
    j = j + 1
  }
  # while(k<=length(tmp_tree$tip.label)){
  #   taxon <- append(taxon, tmp_tree$tip.label[k])
  #   
  #   if (k+1 <= length(tmp_tree$tip.label) && str_split(tmp_tree$tip.label[k], "_")[[1]][1] == str_split(tmp_tree$tip.label[k+1], "_")[[1]][1]){
  #     taxon <- append(taxon, tmp_tree$tip.label[k+1])
  #     strat_ranges <- append(strat_ranges,
  #                            paste0('<stratigraphicRange id="r',
  #                                   j,'" spec="StratigraphicRange" firstOccurrence="@',
  #                                   tmp_tree$tip.label[k+1],'" lastOccurrence="@',
  #                                   tmp_tree$tip.label[k],'"/>'))
  #     if (sample_times_round[k+1]==0.0){
  #       taxon_extant<- append(taxon_extant, tmp_tree$tip.label[k+1])
  #     }
  #     k = k + 2
  #   } else{
  #     strat_ranges <- append(strat_ranges,
  #                            paste0('<stratigraphicRange id="r',
  #                                   j,'" spec="StratigraphicRange" firstOccurrence="@',
  #                                   tmp_tree$tip.label[k],'" lastOccurrence="@',
  #                                   tmp_tree$tip.label[k],'"/>'))
  #     if (sample_times_round[k]==0.0){
  #       taxon_extant<- append(taxon_extant, tmp_tree$tip.label[k])
  #     }
  #     k = k + 1
  #   }
  #   strat_ranges_refs <- append(strat_ranges_refs, paste0('<stratigraphicRange idref="r',j,'"/>'))
  #   j = j + 1
  # }
  
  
  sim <- readLines(paste0(templates_dir, "sRanges_simDNA_template.xml"))
  sim  <- gsub(pattern = "insertNewick",
               replace = paste0("newick='", beast_tree, "'"), x = sim)
  taxon_extant_str = c()
  for (tx in taxon){
    taxon_extant_str <- c(taxon_extant_str, paste0("<sequence spec='Sequence' taxon='",tx,"' value='?'/>"))
  }
  taxon_extant_str <- paste0(taxon_extant_str, collapse = "\n\t\t\t")
  
  taxon_str = c()
  taxon_set_str = c()
  taxa_age_str = c()
  for (tx in taxon){
    taxon_str <- c(taxon_str, paste0("<sequence spec='Sequence' taxon='",tx,"' value='?'/>"))
    taxon_set_str <- c(taxon_set_str, paste0("<taxon spec='Taxon' id='",tx,"'/>"))
    taxa_age_str <- c(taxa_age_str, paste0(tx, "=", sample_times_round[which(tmp_tree$tip.label==tx)]))
  }
  taxon_str <- paste0(taxon_str, collapse = "\n\t\t\t")
  taxon_set_str <- paste0(taxon_set_str, collapse = "\n\t\t\t\t\t\t")
  taxa_age_str <- paste0(taxa_age_str, collapse = ", ")
  
  sim  <- gsub(pattern = "<insertSequence/>",
               replace = taxon_extant_str, x = sim)

  
  sim_dir <- paste0(wd, "run_", i, "/sim")
  dir.create(sim_dir, recursive=T)
  setwd(sim_dir)
  writeLines(sim, con='sRanges_simDNA.xml')
  cmd <- "'/Applications/BEAST\ 2.7.1/bin/beast' sRanges_simDNA.xml"
  system(cmd)
  
  
  extinct <- taxon[which(taxon %notin% taxon_extant)]
  dna_sim <- readLines("simulated_dna_alignment.xml")
  idx  <- grep("<data id=", dna_sim)
  dna_sim[idx] <- "<data id='dna_alignment' spec='beast.base.evolution.alignment.Alignment'>"
  for (tax in extinct){
    idx  <- grep(tax, dna_sim)
    dna_sim[idx] <- paste0("    <sequence spec='beast.base.evolution.alignment.Sequence' taxon='",
                           tax,"' value='", paste0(rep("-", 1000),
                                                   collapse = ""), "'/>")
  }
  dna_sim  <- gsub(pattern = "id='Sequence",
               replace = "id='dna", x = dna_sim)
  
  writeLines(dna_sim, "simulated_dna_alignment.xml")
  
  
  
  ####### simulate morph sequences
  
  sim <- readLines(paste0(templates_dir,"sRanges_simMorph_template.xml"))
  sim  <- gsub(pattern = "insertNewick",
               replace = paste0("newick='", beast_tree, "'"), x = sim)
  taxon_extant_str = c()
  for (tx in taxon){
    taxon_extant_str <- c(taxon_extant_str, paste0("<sequence spec='Sequence' taxon='",tx,"' value='?'/>"))
    # taxon_extant_str<- paste0(taxon_extant_str, "<sequence spec='Sequence' taxon='",tx,"' value='?'/>", collapse = "\n")
  }
  taxon_extant_str <- paste0(taxon_extant_str, collapse = "\n\t\t\t")
  
  taxon_str = c()
  taxon_set_str = c()
  taxa_age_str = c()
  for (tx in taxon){
    taxon_str <- c(taxon_str, paste0("<sequence spec='Sequence' taxon='",tx,"' value='?'/>"))
    taxon_set_str <- c(taxon_set_str, paste0("<taxon spec='Taxon' id='",tx,"'/>"))
    taxa_age_str <- c(taxa_age_str, paste0(tx, "=", sample_times_round[which(tmp_tree$tip.label==tx)]))
  }
  taxon_str <- paste0(taxon_str, collapse = "\n\t\t\t")
  taxon_set_str <- paste0(taxon_set_str, collapse = "\n\t\t\t\t\t\t")
  taxa_age_str <- paste0(taxa_age_str, collapse = ", ")
  
  sim  <- gsub(pattern = "<insertMorphSequence/>",
               replace = taxon_str, x = sim)
  

  writeLines(sim, con='sRanges_simMorph.xml')
  cmd <- "'/Applications/BEAST\ 2.7.1/bin/beast' sRanges_simMorph.xml"
  system(cmd)
  
  morph_sim <- readLines("simulated_morph_alignment.xml")
  idx  <- grep("<data id=", morph_sim)
  morph_sim[idx] <- "<data id='morph_alignment' spec='beast.base.evolution.alignment.Alignment'>"
  morph_sim  <- gsub(pattern = "id='Sequence",
                   replace = "id='morph", x = morph_sim)
  morph_sim  <- gsub(pattern = ",",
                     replace = "", x = morph_sim)
  
  
  sim <- readLines(paste0(templates_dir,"sRanges_inference_template.xml"))
  sim  <- gsub(pattern = "<insertStartMorphData/>",
               replace = paste0(morph_sim, collapse = "\n"), x = sim)
  sim  <- gsub(pattern = "<insertStartDNAData/>",
               replace = paste0(dna_sim, collapse = "\n"), x = sim)

  sim  <- gsub(pattern = "<insertMorphSequence/>",
               replace = taxon_str, x = sim)

  sim  <- gsub(pattern = "<inputTaxa/>",
               replace = taxon_set_str, x = sim)

  sim  <- gsub(pattern = "<inputTaxaAge/>",
               replace = taxa_age_str, x = sim)

  sim  <- gsub(pattern = "<inputStratRanges/>",
               replace = paste0(strat_ranges, collapse = "\n\t\t\t\t"), x = sim)
  sim  <- gsub(pattern = "<inputStratRangesRef/>",
               replace = paste0(strat_ranges_refs, collapse = "\n\t\t\t"), x = sim)


  rnd_origin <- runif(1, mrca*2, 500)
  rnd_div_rate <- rlnorm(1, meanlog = -3.5, sdlog = 1.5)
  rnd_turnover <- runif(1)
  rnd_sampling_prop <- runif(1)
  rnd_sampl_extant_prob <- runif(1, 0.8, 1.)

  sim  <- gsub(pattern = "<initOrigin/>",
               replace = paste0("<parameter id='origin' lower='0.0'
                                name='stateNode'>",rnd_origin,"</parameter>"),
               x = sim)
  sim  <- gsub(pattern = "<initDiversificationRate/>",
               replace = paste0("<parameter id='diversificationRate' lower='0.0'
                                name='stateNode'>",rnd_div_rate,"</parameter>"),
               x = sim)
  sim  <- gsub(pattern = "<initTurnover/>",
               replace = paste0("<parameter id='turnover' lower='0.' upper = '1.'
                                name='stateNode'>",rnd_turnover,"</parameter>"),
               x = sim)
  sim  <- gsub(pattern = "<initSamplingProportion/>",
               replace = paste0("<parameter id='samplingProportion' lower='0.0'
                                name='stateNode'>",rnd_sampling_prop,"</parameter>"),
               x = sim)

  sim  <- gsub(pattern = "<initSamplingAtPresentProb/>",
               replace = paste0("<parameter id='samplingAtPresentProb' lower='0.0'
                                name='stateNode'>",sampl_extant_prob,"</parameter>"),
               x = sim)
  
  inf_dir <- paste0(wd, "run_", i, "/inf/") 
  dir.create(inf_dir, recursive=T)
  setwd(inf_dir)
  writeLines(sim, con="sRanges_inference.xml")
}

write.csv(true_rates, paste0(wd, "true_rates.csv"))

