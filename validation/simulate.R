rm(list = ls())
options(digits=16)

library(ape)
library(phytools)
library(ggplot2)
library(stringr)
library(beastio)
library(FossilSim)

'%notin%' <- Negate('%in%')
wd <- "/Users/jugne/Documents/Source/beast2.7/sRanges-material/"

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

redraw_params <- function(){
  div_rate <- rlnorm(1, meanlog = -3.5, sdlog = 1.5)
  turnover <- runif(1)
  sampling_prop <- runif(1)
  sampl_extant_prob <- runif(1, 0.8, 1.)
  
  lambda <- f_lambda(div_rate, turnover)
  mu <- f_mu(turnover, lambda)
  psi <- f_psi(mu, sampling_prop)
}


div_rate = NULL
turnover = NULL
sampling_prop = NULL
sampl_extant_prob = NULL
lambda = NULL
mu = NULL
psi = NULL
redraws = 0

min_age = 40
max_age = 100

ntrees = 1
nextant = 25

min_fossils = 4
max_fossils = 125

### Simulate Trees

# simulating trees and fossils with parameters
trees = list()
beast_trees = vector(, length=10)
fossils = list()
samp_trees = list()
while(length(trees) < ntrees) {
  redraw_params()
  nsim = ntrees - length(trees)
  trees = c(trees, TreeSim::sim.bd.taxa(nextant, nsim, lambda, mu, complete = T))
  for (i in ntrees:(ntrees-nsim+1)) {
    mrca = max(ape::node.depth.edgelength(trees[[i]]))
    origin = mrca + trees[[i]]$root.edge
    if(origin > max_age || origin < min_age) {
      trees = trees[-i]
      fossils = fossils[-i]
      next
    }
    
    fossils[[i]] = FossilSim::sim.fossils.poisson(psi, tree = trees[[i]])
    if(length(fossils[[i]]$edge) < min_fossils || length(fossils[[i]]$edge) > max_fossils) {
      fossils = fossils[-i]
      trees = trees[-i]
    }
  }
}

for (i in 1:ntrees){
  
  beast_tree = beast.fbd.format(trees[[i]],fossils[[i]], rho=sampl_extant_prob, digits=9)
  tmp_tree = ape::read.tree(text=beast_tree)
  sample_times = ape::node.depth.edgelength(tmp_tree)
  sample_times = max(sample_times) - sample_times
  sample_times_round = round(sample_times,9)
  
  
  # ntaxon = length(t$tip.label)
  taxon = list()
  taxon_extant = list()
  strat_ranges = list()
  strat_ranges_refs = list()
  j=0
  k=1
  # taxon_start = list()
  # taxon_end = list()
  while(k<=length(tmp_tree$tip.label)){
    taxon <- append(taxon, tmp_tree$tip.label[k])
    
    if (k+1 <= length(tmp_tree$tip.label) && str_split(tmp_tree$tip.label[k], "_")[[1]][1] == str_split(tmp_tree$tip.label[k+1], "_")[[1]][1]){
      taxon <- append(taxon, tmp_tree$tip.label[k+1])
      strat_ranges <- append(strat_ranges,
                             paste0('<stratigraphicRange id="r',
                                    j,'" spec="StratigraphicRange" firstOccurrence="@',
                                    tmp_tree$tip.label[k+1],'" lastOccurrence="@',
                                    tmp_tree$tip.label[k],'"/>'))
      if (sample_times_round[k+1]==0.0){
        taxon_extant<- append(taxon_extant, tmp_tree$tip.label[k+1])
      }
      k = k + 2
    } else{
      strat_ranges <- append(strat_ranges,
                             paste0('<stratigraphicRange id="r',
                                    j,'" spec="StratigraphicRange" firstOccurrence="@',
                                    tmp_tree$tip.label[k],'" lastOccurrence="@',
                                    tmp_tree$tip.label[k],'"/>'))
      if (sample_times_round[k]==0.0){
        taxon_extant<- append(taxon_extant, tmp_tree$tip.label[k])
      }
      k = k + 1
    }
    strat_ranges_refs <- append(strat_ranges_refs, paste0('<stratigraphicRange idref="r',j,'"/>'))
    j = j + 1
  }
  
  
  # for (tx in tmp_tree$tip.label){
  #   coord = gregexpr(paste0(tx,"_"), beast_tree)[[1]]
  #   if (coord[1]!=-1){
  #     taxon <- append(taxon, paste0(tx, "_first"))
  #     if (sample_times_round[which(tmp_tree$tip.label==paste0(tx, "_first"))]==0.0){
  #       taxon_extant<- append(taxon_extant, paste0(tx, "_first"))
  #     }
  #     if (length(coord)==2){
  #       taxon <- append(taxon, paste0(tx, "_last"))
  #       if (sample_times_round[which(tmp_tree$tip.label==paste0(tx, "_last"))]==0.0){
  #         taxon_extant<- append(taxon_extant, paste0(tx, "_last"))
  #       }
  #       strat_ranges <- append(strat_ranges,
  #                             paste0('<stratigraphicRange id="r',
  #                                    j,'" spec="StratigraphicRange" firstOccurrence="@',
  #                                    paste0(tx, "_first"),'" lastOccurrence="@',
  #                                    paste0(tx, "_last"),'"/>'))
  #     } else {
  #       strat_ranges <- append(strat_ranges,
  #                             paste0('<stratigraphicRange id="r',
  #                                    j,'" spec="StratigraphicRange" firstOccurrence="@',
  #                                    paste0(tx, "_first"),'" lastOccurrence="@',
  #                                    paste0(tx, "_first"),'"/>'))
  #     }
  #     strat_ranges_refs <- append(strat_ranges_refs, paste0('<stratigraphicRange idref="r',j,'"/>'))
  #     j=j+1
  #   }
  # }
  
  # extinct <- taxon[which(taxon %notin% taxon_extant)]
  # tmp_tree_extant <- ape::drop.tip(tmp_tree, extinct)
  
  sim <- readLines("/Users/jugne/Documents/Source/beast2.7/sRanges-material/sRanges_simDNA_template.xml")
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
    
    
    # taxon_str <- paste(taxon_str, "<sequence spec='Sequence' taxon='",tx,"' value='?'/>", collapse = "\n")
    # taxon_set_str <- paste(taxon_set_str, "<taxon spec='Taxon' id='",tx,"'/>", sep="", collapse = "\n")
    # tmp_line <- paste0(tx, "=", sample_times_round[which(tmp_tree$tip.label==tx)])
    # taxa_age_str <- paste(taxa_age_str, tmp_line, collapse = ",")
  }
  taxon_str <- paste0(taxon_str, collapse = "\n\t\t\t")
  taxon_set_str <- paste0(taxon_set_str, collapse = "\n\t\t\t\t\t\t")
  taxa_age_str <- paste0(taxa_age_str, collapse = ", ")
  
  sim  <- gsub(pattern = "<insertSequence/>",
               replace = taxon_extant_str, x = sim)
  # 
  # sim  <- gsub(pattern = "<insertMorphSequence/>",
  #              replace = taxon_str, x = sim)
  # 
  # sim  <- gsub(pattern = "<inputTaxa/>",
  #              replace = taxon_set_str, x = sim)
  # 
  # sim  <- gsub(pattern = "<inputTaxaAge/>",
  #              replace = taxa_age_str, x = sim)
  # 
  # sim  <- gsub(pattern = "<inputStratRanges/>",
  #              replace = paste0(strat_ranges, collapse = "\n\t\t\t\t"), x = sim)
  # sim  <- gsub(pattern = "<inputStratRangesRef/>",
  #              replace = paste0(strat_ranges_refs, collapse = "\n\t\t\t"), x = sim)
  
  
  # rnd_origin <- runif(1, mrca, mrca+1000)
  # rnd_birth <- runif(1, 0.0001, 1000.0001)
  # rnd_death <- runif(1, 0.0001, 1000.0001)
  # rnd_sample <- runif(1, 0.0001, 1000.0001)
  # rnd_sampl_extant_prob <- runif(1, 0., 1.)
  # 
  # sim  <- gsub(pattern = "<initOrigin/>",
  #              replace = paste0("<parameter id='origin' lower='0.0' 
  #                               name='stateNode'>",rnd_origin,"</parameter>"),
  #              x = sim)
  # sim  <- gsub(pattern = "<initBirthRate/>",
  #              replace = paste0("<parameter id='birthRate' lower='0.0'
  #                               name='stateNode'>",rnd_birth,"</parameter>"),
  #              x = sim)
  # sim  <- gsub(pattern = "<initDeathRate/>",
  #              replace = paste0("<parameter id='deathRate' lower='0.0'
  #                               name='stateNode'>",rnd_death,"</parameter>"),
  #              x = sim)
  # sim  <- gsub(pattern = "<initSamplingRate/>",
  #              replace = paste0("<parameter id='samplingRate' lower='0.0'
  #                               name='stateNode'>",rnd_sample,"</parameter>"),
  #              x = sim)
  # 
  # sim  <- gsub(pattern = "<initSamplingAtPresentProb/>",
  #              replace = paste0("<parameter id='samplingAtPresentProb' lower='0.0'
  #                               name='stateNode'>",rnd_sampl_extant_prob,"</parameter>"),
  #              x = sim)
  
  
  
  
  setwd(wd)
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
                           tax,"' value='", paste0(rep("-", 2000),
                                                   collapse = ""), "'/>")
  }
  dna_sim  <- gsub(pattern = "id='Sequence",
               replace = "id='dna", x = dna_sim)
  
  writeLines(dna_sim, con=paste0(wd, "simulated_dna_alignment.xml"))
  
  
  
  ####### simulate morph sequences
  
  sim <- readLines("/Users/jugne/Documents/Source/beast2.7/sRanges-material/sRanges_simMorph_template.xml")
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
    
    
    # taxon_str <- paste(taxon_str, "<sequence spec='Sequence' taxon='",tx,"' value='?'/>", collapse = "\n")
    # taxon_set_str <- paste(taxon_set_str, "<taxon spec='Taxon' id='",tx,"'/>", sep="", collapse = "\n")
    # tmp_line <- paste0(tx, "=", sample_times_round[which(tmp_tree$tip.label==tx)])
    # taxa_age_str <- paste(taxa_age_str, tmp_line, collapse = ",")
  }
  taxon_str <- paste0(taxon_str, collapse = "\n\t\t\t")
  taxon_set_str <- paste0(taxon_set_str, collapse = "\n\t\t\t\t\t\t")
  taxa_age_str <- paste0(taxa_age_str, collapse = ", ")
  
  sim  <- gsub(pattern = "<insertMorphSequence/>",
               replace = taxon_str, x = sim)
  
  setwd(wd)
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
  
  
  sim <- readLines("sRanges_inference_template.xml")
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


  rnd_origin <- runif(1, mrca, mrca+1000)
  rnd_birth <- runif(1, 0.0001, 1000.0001)
  rnd_death <- runif(1, 0.0001, 1000.0001)
  rnd_sample <- runif(1, 0.0001, 1000.0001)
  rnd_sampl_extant_prob <- runif(1, 0., 1.)

  sim  <- gsub(pattern = "<initOrigin/>",
               replace = paste0("<parameter id='origin' lower='0.0'
                                name='stateNode'>",rnd_origin,"</parameter>"),
               x = sim)
  sim  <- gsub(pattern = "<initBirthRate/>",
               replace = paste0("<parameter id='birthRate' lower='0.0'
                                name='stateNode'>",rnd_birth,"</parameter>"),
               x = sim)
  sim  <- gsub(pattern = "<initDeathRate/>",
               replace = paste0("<parameter id='deathRate' lower='0.0'
                                name='stateNode'>",rnd_death,"</parameter>"),
               x = sim)
  sim  <- gsub(pattern = "<initSamplingRate/>",
               replace = paste0("<parameter id='samplingRate' lower='0.0'
                                name='stateNode'>",rnd_sample,"</parameter>"),
               x = sim)

  sim  <- gsub(pattern = "<initSamplingAtPresentProb/>",
               replace = paste0("<parameter id='samplingAtPresentProb' lower='0.0'
                                name='stateNode'>",rnd_sampl_extant_prob,"</parameter>"),
               x = sim)
  
  writeLines(sim, con="sRanges_inference.xml")
}



# for each tree
# get unique number of species
# for each species find first and last fossil
# use this for taxons 


