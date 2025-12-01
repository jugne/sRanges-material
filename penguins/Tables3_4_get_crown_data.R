library(beastio)
library(treestats)
library(phytools)
library(dispRity)
library(stringr)
library(coda)

source("penguin_paths.R")
ranges_trees <- paste0(sRangesPath, "/penguins_inf_morph_at_start.combined.trees")
fbd_trees <- paste0(fbdUpdatedPath, "/Mkv_morph+dna.combined.trees")


tt <-readTreeLog(ranges_trees,
                 burnin = 0)
ntrees<- length(tt)
species<- c("Madrynornis_mirandus", "Marplesornis_novaezealandiae",
            "Eretiscus_tonnii", "Palaeospheniscus_bergi",
            "Palaeospheniscus_biloculata", "Pygoscelis_grandis",
            "Palaeospheniscus_patagonicus", "Spheniscus_megaramphus",
            "Spheniscus_muizoni", "Spheniscus_urbinai")
tips <- tt[[1]]$tip.label
ll <- sapply(species, str_c, "_last")[sapply(species, str_c, "_last") %in% tips]
ss <- c(sapply(species, str_c, "_first"), ll)

crown_age<-numeric(ntrees)
crown_taxa<-character(ntrees)
dd <- data.frame(matrix(NA, nrow = ntrees, ncol = length(ss)))
colnames(dd) <- ss
for (i in 1:ntrees){
  cr<-crown.stem(tt[[i]], inc.nodes = F, output.names = TRUE)$crown
  crown_age[i] <- max(nodeHeights(tt[[i]]))-nodeheight(tt[[i]],
                                                       getMRCA(tt[[i]],tip=cr))
  dd[i,]<- as.numeric(ss %in% cr)
}

hpd<-HPDinterval(as.mcmc(crown_age))
dd_sort <- sort(round(colMeans(dd),2), decreasing = T)
write.csv(data.frame(crown_age_and_probs=c(mean(crown_age), hpd[1], hpd[2], dd_sort), 
                     row.names = c("crown_age_mean", "crown_age_hpd_l", "crown_age_hpd_h",
                                   names(dd_sort))),
          paste0(table3_4_dir_dir,"/ranges_crown.csv"))


tt <-readTreeLog(fbd_trees,
                 burnin = 0)
ntrees<- length(tt)
species<- c("Madrynornis_mirandus", "Marplesornis_novaezealandiae",
            "Eretiscus_tonnii", "Palaeospheniscus_bergi",
            "Palaeospheniscus_biloculata", "Pygoscelis_grandis",
            "Palaeospheniscus_patagonicus", "Spheniscus_megaramphus",
            "Spheniscus_muizoni", "Spheniscus_urbinai")

dd <- data.frame(matrix(NA, nrow = ntrees, ncol = length(species)))
colnames(dd) <- species
for (i in 1:ntrees){
  cr<-crown.stem(tt[[i]], inc.nodes = F, output.names = TRUE)$crown
  crown_age[i] <- max(nodeHeights(tt[[i]]))-nodeheight(tt[[i]],
                                                       getMRCA(tt[[i]],tip=cr))
  dd[i,]<- as.numeric(species %in% cr)
}

hpd<-HPDinterval(as.mcmc(crown_age))
dd_sort <- sort(round(colMeans(dd),2), decreasing = T)
write.csv(data.frame(crown_age_and_probs=c(mean(crown_age),hpd[1], hpd[2],  dd_sort),
                     row.names = c("crown_age_mean", "crown_age_hpd_l", "crown_age_hpd_h",
                                   names(dd_sort))),
          paste0(table3_4_dir_dir,"/fbd_crown.csv"))

