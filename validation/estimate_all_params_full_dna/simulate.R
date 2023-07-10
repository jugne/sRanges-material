rm(list = ls())
options(digits=16)

library(ape)
library(phytools)
library(ggplot2)
library(stringr)
library(beastio)
library(FossilSim)
library(R.utils)
source("/Users/jugne/Documents/Source/beast2.7/sRanges-material/validation/estimate_all_params_full_dna/helper.R")

tmpfun <- get("sim2.bd.origin", envir = asNamespace("TreeSim"))
environment(sim2.bd.origin2) <- environment(tmpfun)
attributes(sim2.bd.origin2) <- attributes(tmpfun)  # don't know if this is really needed
assignInNamespace("sim2.bd.origin", sim2.bd.origin2, ns="TreeSim")

'%notin%' <- Negate('%in%')
wd <- "~/Documents/Source/beast2.7/sRanges-material/validation/estimate_all_params_full_dna/"
templates_dir <- "~/Documents/Source/beast2.7/sRanges-material/validation/estimate_all_params_full_dna/templates/"


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
origin = 4

redraw_params <- function(){
  div_rate <<- runif(1, 0.7, 0.9)#rlnorm(1, meanlog = -3.5, sdlog = 1.5)
  turnover <<- runif(1, 0.2, 0.8)
  sampling_prop <<- runif(1, 0.2, 0.8)
  sampl_extant_prob <<- runif(1, 0.7, 1.)
  
  lambda <<- f_lambda(div_rate, turnover)
  mu <<- f_mu(turnover, lambda)
  psi <<- f_psi(mu, sampling_prop)
  # origin <<- runif(1, 40, 80)
}

ntrees = 200

min_ext_samples = 5
max_samples = 1000

min_fossils = 0
max_fossils = 1000

true_rates <- data.frame(div_rate=numeric(), turnover=numeric(),
                         sampling_prop=numeric(), rho=numeric(),
                         origin=numeric(), mrca=numeric(), tree=character(),
                         n_samples=numeric(), n_extant=numeric(), n_ranges=numeric(), draws = numeric())

### Simulate Trees

# simulating trees and fossils with parameters
trees = list()
beast_trees = list()
fossils = list()
taxonomy = list()
samp_trees = list()
while(length(trees) < ntrees) {
  redraw_params()
  redraws = redraws+1
  tree_tmp <- TreeSim::sim.bd.age(origin, 1, lambda, mu)
  
  if (length(tree_tmp[[1]])==1 || tree_tmp[[1]]$Nnode>1000){
    next
  }
  mrca = max(ape::node.depth.edgelength(tree_tmp[[1]]))
  if (length(which((mrca - ape::node.depth.edgelength(tree_tmp[[1]])) < 1e-7)) < min_ext_samples){
    next
  }
  
  taxonomy_tmp = sim.taxonomy(tree_tmp[[1]], beta=0, lambda.a = 0)
  fossils_tmp = FossilSim::sim.fossils.poisson(psi, taxonomy = taxonomy_tmp)
  
  # if(#length(tree_tmp[[1]]$tip.label) < min_samples ||
  #   length(tree_tmp[[1]]$tip.label) > max_samples) {
  #   next
  # }
  
  beast_tree_tmp = beast.fbd.format(tree_tmp[[1]], fossils_tmp, rho=sampl_extant_prob, digits=16)
  tree_after_rho = ape::read.tree(text=beast_tree_tmp)
  mrca = max(ape::node.depth.edgelength(tree_after_rho))
  n_ext <- length(which((mrca - ape::node.depth.edgelength(tree_after_rho)) < 1e-7))
  n_fossils <- length(tree_after_rho$tip.label) - n_ext
  if (n_ext < min_ext_samples || n_fossils < min_fossils || n_fossils > max_fossils){
    next
  }
  
  trees = c(trees, tree_tmp)
  beast_trees[[length(trees)]] = beast_tree_tmp
  fossils[[length(trees)]] = fossils_tmp
  taxonomy = c(taxonomy, taxonomy_tmp)
  true_rates_tmp <- data.frame(div_rate, turnover, sampling_prop,
                               rho=sampl_extant_prob, origin, mrca, tree=beast_tree_tmp,
                               n_samples=n_fossils+n_ext, n_extant=n_ext, n_ranges=0, draws=redraws)
  true_rates = rbind(true_rates, true_rates_tmp)
  redraws = 0
  print(n_fossils+n_ext)
  
}

  for (i in 1:ntrees){
  
  beast_tree = beast_trees[[i]]#beast.fbd.format(trees[[i]],fossils[[i]], rho=true_rates$rho[i], digits=16)
  true_rates$tree[i] = beast_tree
  tmp_tree = ape::read.tree(text=beast_tree)
  # true_rates$mrca[i] = max(ape::node.depth.edgelength(tmp_tree))
  # true_rates$n_samples[i] = length(tmp_tree$tip.label)
  # true_rates$n_extant[i] = sum(true_rates$mrca[i]-ape::node.depth.edgelength(tmp_tree)<1e-10)
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
  true_rates$n_ranges[i] = j
  
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
  cmd <- "'/Applications/BEAST\ 2.7.1/bin/beast' -seed 42 sRanges_simDNA.xml"
  system(cmd)


  extinct <- taxon[which(taxon %notin% taxon_extant)]
  dna_sim <- readLines("simulated_dna_alignment.xml")
  idx  <- grep("<data id=", dna_sim)
  dna_sim[idx] <- "<data id='dna_alignment' spec='beast.base.evolution.alignment.Alignment'>"
  # for (tax in extinct){
  #   idx  <- grep(tax, dna_sim)
  #   dna_sim[idx] <- paste0("    <sequence spec='beast.base.evolution.alignment.Sequence' taxon='",
  #                          tax,"' value='", paste0(rep("-", 1000),
  #                                                  collapse = ""), "'/>")
  # }
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
  cmd <- "'/Applications/BEAST\ 2.7.1/bin/beast' -seed 42 sRanges_simMorph.xml"
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
  
  
  rnd_origin <- runif(1, mrca*2, 500) # just so origin is not smaller than mrca
  rnd_div_rate <- runif(1, 0.7, 0.9)
  rnd_turnover <- runif(1, 0.2, 0.8)
  rnd_sampling_prop <- runif(1, 0.2, 0.8)
  rnd_sampl_extant_prob <- runif(1, 0.7, 1.) 
  
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
                                name='stateNode'>",rnd_sampl_extant_prob,"</parameter>"),
               x = sim)
  
  inf_dir <- paste0(wd, "run_", i, "/inf/")
  dir.create(inf_dir, recursive=T)
  setwd(inf_dir)
  writeLines(sim, con="sRanges_inference.xml")
}

save.image(file=paste0(wd, "simulation.RData"))

write.csv(true_rates, paste0(wd, "true_rates.csv"))


