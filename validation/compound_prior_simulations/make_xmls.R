library(stringr)
'%notin%' <- Negate('%in%')

wd <-  "~/Documents/Source/beast2.7/sRanges-material/validation/compound_prior_simulations/"
templates_dir<- paste0(wd, "templates/") 
lines<-readLines("~/Documents/Source/beast2.7/sRanges-material/validation/compound_prior_simulations/SampleTrees.txt")
trees <- str_split(lines, ";")
param_lines<- readLines("~/Documents/Source/beast2.7/sRanges-material/validation/compound_prior_simulations/SampleTaxa.txt")
tree_taxa_id_ends <- which(grepl(";", param_lines))
param_lines[tree_taxa_id_ends]<-gsub(";", "", param_lines[tree_taxa_id_ends])

getD <- function(rho){
  return(c((log(10)/rho)/100, (log(100)/rho)/100))
}
getTurnover <- function(d, origin){
  z <- exp(d*origin)
  return(c((10+1)/(10+z), (50+1)/(50+z)))
}

rho_bounds <- c(0.5, 1)
s_bounds <- c(0.2, 0.5)



start<-1
for (t in 1:length(trees)){
  end <- tree_taxa_id_ends[t]
  true_tree <- gsub("\\[&reaction=psi\\]", "", trees[[t]][1])
  true_tree <- gsub("\\[&rho=true\\]", "", true_tree)
  true_tree <- gsub("\\[&reaction=lambda\\]", "", true_tree)
  taxon <- str_split_i(param_lines[start:end], "=", 1)
  for (k in 1:length(unique(sub("_.*", "", taxon)))){
    tip <- unique(sub("_.*", "", taxon))[k]
    if (length(which(sub("_.*", "", taxon) %in% tip))==1){
      val<- str_extract(true_tree,paste0("[^0-9]",tip, "_last"))
      true_tree <- gsub(paste0("[^0-9]",tip, "_last"), paste0(substr(val, 1,1),tip, "_first"), true_tree)
      taxon[which(sub("_.*", "", taxon) %in% tip)] <- paste0(tip, "_first")
    }
  }
  dates <- str_split_i(param_lines[start:end], "=", 2)
  taxon_extant<-taxon[which(dates=="0.0")]
  start<-end+1

  # true_tree <- gsub("[]", "", true_tree)
  
  ################# dna sim ###############
  sim <- readLines(paste0(templates_dir, "sRanges_simDNA_template.xml"))
  sim  <- gsub(pattern = "insertNewick",
               replace = paste0("newick='", true_tree, "'"), x = sim)
  taxon_extant_str = c()
  for (tx in taxon){
    taxon_extant_str <- c(taxon_extant_str, paste0("<sequence spec='Sequence' taxon='",tx,"' value='?'/>"))
  }
  taxon_extant_str <- paste0(taxon_extant_str, collapse = "\n\t\t\t")
  
  taxon_str = c()
  taxon_set_str = c()
  taxa_age_str = c()
  for (tx_i in 1:length(taxon)){
    tx<-taxon[tx_i]
    taxon_str <- c(taxon_str, paste0("<sequence spec='Sequence' taxon='",tx,"' value='?'/>"))
    taxon_set_str <- c(taxon_set_str, paste0("<taxon spec='Taxon' id='",tx,"'/>"))
    taxa_age_str <- c(taxa_age_str, paste0(tx, "=", dates[tx_i]))
  }
  taxon_str <- paste0(taxon_str, collapse = "\n\t\t\t")
  taxon_set_str <- paste0(taxon_set_str, collapse = "\n\t\t\t\t\t\t")
  taxa_age_str <- paste0(taxa_age_str, collapse = ", ")
  
  sim  <- gsub(pattern = "<insertSequence/>",
               replace = taxon_extant_str, x = sim)
  
  
  sim_dir <- paste0(wd, "run_", t, "/sim")
  dir.create(sim_dir, recursive=T)
  setwd(sim_dir)
  writeLines(sim, con='sRanges_simDNA.xml')
  cmd <- "'/Applications/BEAST\ 2.7.1/bin/beast' -seed 42 sRanges_simDNA.xml"
  system(cmd)
  
  
  extinct <- taxon[which(taxon %notin% taxon_extant)]
  dna_sim <- readLines("simulated_dna_alignment.xml")
  idx  <- grep("<data id=", dna_sim)
  dna_sim[idx] <- "<data id='dna_alignment' spec='beast.base.evolution.alignment.Alignment'>"
  for (tax in extinct){
    idx  <- grep(paste0("'",tax), dna_sim)
    dna_sim[idx] <- paste0("    <sequence spec='beast.base.evolution.alignment.Sequence' taxon='",
                           tax,"' value='", paste0(rep("-", 1000),
                                                   collapse = ""), "'/>")
  }
  dna_sim  <- gsub(pattern = "id='Sequence",
                   replace = "id='dna", x = dna_sim)
  
  writeLines(dna_sim, "simulated_dna_alignment.xml")
  
  ################# morph sim ###############
  
  sim <- readLines(paste0(templates_dir,"sRanges_simMorph_template.xml"))
  sim  <- gsub(pattern = "insertNewick",
               replace = paste0("newick='", true_tree, "'"), x = sim)
  taxon_extant_str = c()
  for (tx in taxon){
    taxon_extant_str <- c(taxon_extant_str, paste0("<sequence spec='Sequence' taxon='",tx,"' value='?'/>"))
    # taxon_extant_str<- paste0(taxon_extant_str, "<sequence spec='Sequence' taxon='",tx,"' value='?'/>", collapse = "\n")
  }
  taxon_extant_str <- paste0(taxon_extant_str, collapse = "\n\t\t\t")
  
  taxon_str = c()
  taxon_set_str = c()
  taxa_age_str = c()
  for (tx_i in 1:length(taxon)){
    tx<-taxon[tx_i]
    taxon_str <- c(taxon_str, paste0("<sequence spec='Sequence' taxon='",tx,"' value='?'/>"))
    taxon_set_str <- c(taxon_set_str, paste0("<taxon spec='Taxon' id='",tx,"'/>"))
    taxa_age_str <- c(taxa_age_str, paste0(tx, "=", dates[tx_i]))
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
  
  
  ###################################################################
  strat_ranges = list()
  strat_ranges_refs = list()
  j=0
  k=1
  for (k in 1:length(unique(sub("_.*", "", taxon)))){
    tip = unique(sub("_.*", "", taxon))[k]
    if (length(which(sub("_.*", "", taxon) %in% tip))>1){
      strat_ranges <- append(strat_ranges,
                             paste0('<stratigraphicRange id="r',
                                    j,'" spec="StratigraphicRange" firstOccurrence="@',
                                    paste0(tip, "_first"),'" lastOccurrence="@',
                                    paste0(tip, "_last"),'"/>'))
    } else{
      strat_ranges <- append(strat_ranges,
                             paste0('<stratigraphicRange id="r',
                                    j,'" spec="StratigraphicRange" firstOccurrence="@',
                                    paste0(tip, "_first"),'" lastOccurrence="@',
                                    paste0(tip, "_first"),'"/>'))
    }
    strat_ranges_refs <- append(strat_ranges_refs, paste0('<stratigraphicRange idref="r',j,'"/>'))
    j = j + 1
  }
  

  
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
  
  
  rnd_origin <- runif(1, length(taxon)+50, 600) # just so origin is not smaller than mrca
  rnd_sampling_prop <- runif(1, s_bounds)
  rnd_sampl_extant_prob <- runif(1, rho_bounds) 
  rnd_div_rate <- runif(1, getD(rnd_sampl_extant_prob))
  rnd_turnover <- runif(1, getTurnover(rnd_div_rate, origin))
  # rnd_sampling_prop <- runif(1, 0.2, 0.5)
  # rnd_sampl_extant_prob <- runif(1, 0.5, 1.) 
  
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
  
  inf_dir <- paste0(wd, "run_", t, "/inf/")
  dir.create(inf_dir, recursive=T)
  setwd(inf_dir)
  writeLines(sim, con="sRanges_inference.xml")
  
}


