rm(list = ls())
gc()

library(stringr)
library(deeptime)
library(coda)
library(ggridges)
library(ggplot2)

source("penguin_paths.R")

# Define the text containing the sampling dates information
text <- '<samplingDates id="samplingDate1.0" spec="sa.evolution.tree.SamplingDate" taxon="@Anthropornis_grandis_last" lower="34" upper="39.5"/>
        <samplingDates id="samplingDate2.0" spec="sa.evolution.tree.SamplingDate" taxon="@Anthropornis_nordenskjoeldi_last" lower="34" upper="39"/>
        <samplingDates id="samplingDate3.0" spec="sa.evolution.tree.SamplingDate" taxon="@Archaeospheniscus_lopdelli_first" lower="25.2" upper="27.3"/>
        <samplingDates id="samplingDate4.0" spec="sa.evolution.tree.SamplingDate" taxon="@Archaeospheniscus_lowei_first" lower="25.2" upper="27.3"/>
        <samplingDates id="samplingDate5.0" spec="sa.evolution.tree.SamplingDate" taxon="@Burnside_Palaeeudyptes_first" lower="36" upper="37.2"/>
        <samplingDates id="samplingDate6.0" spec="sa.evolution.tree.SamplingDate" taxon="@Delphinornis_arctowskii_first" lower="34" upper="39.5"/>
        <samplingDates id="samplingDate7.0" spec="sa.evolution.tree.SamplingDate" taxon="@Delphinornis_gracilis_first" lower="34" upper="39.5"/>
        <samplingDates id="samplingDate8.0" spec="sa.evolution.tree.SamplingDate" taxon="@Delphinornis_larseni_first" lower="39.5" upper="41"/>
        <samplingDates id="samplingDate8.1" spec="sa.evolution.tree.SamplingDate" taxon="@Delphinornis_larseni_last" lower="34" upper="39.5"/>
        <samplingDates id="samplingDate9.0" spec="sa.evolution.tree.SamplingDate" taxon="@Delphinornis_wimani_last" lower="34" upper="39.5"/>
        <samplingDates id="samplingDate10.0" spec="sa.evolution.tree.SamplingDate" taxon="@Duntroonornis_parvus_first" lower="25.2" upper="27.3"/>
        <samplingDates id="samplingDate11.0" spec="sa.evolution.tree.SamplingDate" taxon="@Eretiscus_tonnii_first" lower="18.2" upper="19.8"/>
        <samplingDates id="samplingDate12.0" spec="sa.evolution.tree.SamplingDate" taxon="@Icadyptes_salasi_first" lower="35.7" upper="37.2"/>
        <samplingDates id="samplingDate13.0" spec="sa.evolution.tree.SamplingDate" taxon="@Inkayacu_paracasensis_first" lower="35.7" upper="37.2"/>
        <samplingDates id="samplingDate14.0" spec="sa.evolution.tree.SamplingDate" taxon="@Kairuku_grebneffi_first" lower="25.2" upper="27.3"/>
        <samplingDates id="samplingDate15.0" spec="sa.evolution.tree.SamplingDate" taxon="@Kairuku_waitaki_first" lower="25.2" upper="28.1"/>
        <samplingDates id="samplingDate16.0" spec="sa.evolution.tree.SamplingDate" taxon="@Madrynornis_mirandus_first" lower="9.7" upper="10.3"/>
        <samplingDates id="samplingDate17.0" spec="sa.evolution.tree.SamplingDate" taxon="@Marambiornis_exilis_last" lower="34" upper="39.5"/>
        <samplingDates id="samplingDate18.0" spec="sa.evolution.tree.SamplingDate" taxon="@Marplesornis_novaezealandiae_first" lower="1.8" upper="13.8"/>
        <samplingDates id="samplingDate19.0" spec="sa.evolution.tree.SamplingDate" taxon="@Mesetaornis_polaris_last" lower="34" upper="39.5"/>
        <samplingDates id="samplingDate20.0" spec="sa.evolution.tree.SamplingDate" taxon="@Pachydyptes_ponderosus_first" lower="34.3" upper="36"/>
        <samplingDates id="samplingDate21.0" spec="sa.evolution.tree.SamplingDate" taxon="@Palaeeudyptes_antarcticus_first" lower="36" upper="37.2"/>
        <samplingDates id="samplingDate21.1" spec="sa.evolution.tree.SamplingDate" taxon="@Palaeeudyptes_antarcticus_last" lower="21.7" upper="25.2"/>
        <samplingDates id="samplingDate22.0" spec="sa.evolution.tree.SamplingDate" taxon="@Palaeeudyptes_gunnari_first" lower="41.5" upper="44"/>
        <samplingDates id="samplingDate22.1" spec="sa.evolution.tree.SamplingDate" taxon="@Palaeeudyptes_gunnari_last" lower="34" upper="39.5"/>
        <samplingDates id="samplingDate23.0" spec="sa.evolution.tree.SamplingDate" taxon="@Palaeeudyptes_klekowskii_last" lower="34" upper="39.5"/>
        <samplingDates id="samplingDate24.0" spec="sa.evolution.tree.SamplingDate" taxon="@Palaeospheniscus_bergi_last" lower="11.61" upper="15.97"/>
        <samplingDates id="samplingDate24.1" spec="sa.evolution.tree.SamplingDate" taxon="@Palaeospheniscus_bergi_first" lower="18.2" upper="24.01"/>
        <samplingDates id="samplingDate25.0" spec="sa.evolution.tree.SamplingDate" taxon="@Palaeospheniscus_biloculata_first" lower="18.2" upper="19.8"/>
        <samplingDates id="samplingDate26.0" spec="sa.evolution.tree.SamplingDate" taxon="@Palaeospheniscus_patagonicus_first" lower="18.2" upper="19.8"/>
        <samplingDates id="samplingDate27.0" spec="sa.evolution.tree.SamplingDate" taxon="@Paraptenodytes_antarcticus_first" lower="18.2" upper="19.8"/>
        <samplingDates id="samplingDate27.1" spec="sa.evolution.tree.SamplingDate" taxon="@Paraptenodytes_antarcticus_last" lower="16.3" upper="17.5"/>
        <samplingDates id="samplingDate28.0" spec="sa.evolution.tree.SamplingDate" taxon="@Perudyptes_devriesi_first" lower="41.3" upper="47.8"/>
        <samplingDates id="samplingDate29.0" spec="sa.evolution.tree.SamplingDate" taxon="@Platydyptes_novaezealandiae_last" lower="25.2" upper="27.3"/>
        <samplingDates id="samplingDate29.1" spec="sa.evolution.tree.SamplingDate" taxon="@Platydyptes_novaezealandiae_first" lower="34.3" upper="36"/>
        <samplingDates id="samplingDate30.0" spec="sa.evolution.tree.SamplingDate" taxon="@Platydyptes_marplesi_first" lower="25.2" upper="27.3"/>
        <samplingDates id="samplingDate31.0" spec="sa.evolution.tree.SamplingDate" taxon="@Pygoscelis_grandis_first" lower="4.5" upper="8.6"/>
        <samplingDates id="samplingDate31.1" spec="sa.evolution.tree.SamplingDate" taxon="@Pygoscelis_grandis_last" lower="2.6" upper="4.5"/>
        <samplingDates id="samplingDate32.0" spec="sa.evolution.tree.SamplingDate" taxon="@Spheniscus_megaramphus_first" lower="5.3" upper="11.6"/>
        <samplingDates id="samplingDate32.1" spec="sa.evolution.tree.SamplingDate" taxon="@Spheniscus_megaramphus_last" lower="5.3" upper="7.2"/>
        <samplingDates id="samplingDate33.0" spec="sa.evolution.tree.SamplingDate" taxon="@Spheniscus_muizoni_first" lower="9" upper="11.6"/>
        <samplingDates id="samplingDate33.1" spec="sa.evolution.tree.SamplingDate" taxon="@Spheniscus_muizoni_last" lower="3.6" upper="9"/>
        <samplingDates id="samplingDate34.0" spec="sa.evolution.tree.SamplingDate" taxon="@Spheniscus_urbinai_last" lower="5.3" upper="9"/>
        <samplingDates id="samplingDate34.1" spec="sa.evolution.tree.SamplingDate" taxon="@Spheniscus_urbinai_first" lower="9" upper="11.6"/>
        <samplingDates id="samplingDate35.0" spec="sa.evolution.tree.SamplingDate" taxon="@Waimanu_manneringi_first" lower="60.5" upper="61.6"/>
        <samplingDates id="samplingDate36.0" spec="sa.evolution.tree.SamplingDate" taxon="@Muriwaimanu_tuatahi_first" lower="58" upper="60.5"/>
        <samplingDates id="samplingDate36.1" spec="sa.evolution.tree.SamplingDate" taxon="@Muriwaimanu_tuatahi_last" lower="55.8" upper="58.7"/>
        <samplingDates id="samplingDate37.0" spec="sa.evolution.tree.SamplingDate" taxon="@Eudyptula_minor_first" lower="0.0117" upper="0.126"/>
        <samplingDates id="samplingDate38.0" spec="sa.evolution.tree.SamplingDate" taxon="@Spheniscus_humboldti_first" lower="5.3" upper="7.2"/>
        <samplingDates id="samplingDate39.0" spec="sa.evolution.tree.SamplingDate" taxon="@Spheniscus_demersus_first" lower="0.126" upper="0.781"/>
        <samplingDates id="samplingDate40.0" spec="sa.evolution.tree.SamplingDate" taxon="@Megadyptes_antipodes_first" lower="0.0117" upper="0.126"/>
        <samplingDates id="samplingDate41.0" spec="sa.evolution.tree.SamplingDate" taxon="@Eudyptes_pachyrhynchus_first" lower="0.0117" upper="0.126"/>
'
# Use regular expressions to extract taxon names without the "@", and the lower and upper values
sp <- gsub('@', '', regmatches(text, gregexpr('(?<=taxon="@)[^"]+', text, perl=TRUE))[[1]])
start<- as.numeric(regmatches(text, gregexpr('(?<=upper=")[^"]+', text, perl=TRUE))[[1]])
end<- as.numeric(regmatches(text, gregexpr('(?<=lower=")[^"]+', text, perl=TRUE))[[1]])



# Removing '_first' or '_last' from taxon names if needed
ranges_ <- gsub('_(first|last)$', '', sp)

priors <- data.frame(species=sp, start=start, end=end)

all_string_listwise <- unlist(lapply(ranges_, unique))
ranges_single <- names(which(table(all_string_listwise)==1))
ranges <- names(which(table(all_string_listwise)>1))

id_first <- which(ranges_ %in% ranges_single)
sp_first <- data.frame("Species"=ranges_[id_first], "First start"=start[id_first], "First end"=end[id_first],
                       "Last start"=rep("-", length(id_first)), "Last end"=rep("-", length(id_first)))

sp_both <- data.frame("Species"=c(), "First start"=c(), "First end"=c(),
                       "Last start"=c(), "Last end"=c())

for (r in ranges){
  id1 <- which(sp==paste0(r,"_first"))
  id2 <- which(sp==paste0(r,"_last"))
  sp_both <- rbind(sp_both, 
                   data.frame("Species"=r, "First start"=start[id1], "First end"=end[id1],
                              "Last start"=start[id2], "Last end"=end[id2]))
}

write.csv(rbind(sp_first, sp_both),paste0(figure13a_dir,"/penguin_range_dates.csv"))

get_sRange_full <- function(log){
  start_r <-c()
  start_r_low <- c()
  start_r_high <- c()
  start_r_min <- c()
  start_r_max <- c()
  start_r_iqr <- c()
  start_r_iqr_prior <- c()
  end_r <-c()
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
  median_s_prior<-c()
  median_e <- c()
  median_e_prior<-c()
  
  for (r in ranges_single){
    if ((max(log[,paste0(r,"_first")])-min(log[,paste0(r,"_first")]))>0.1){
      start_l <- length(log[,paste0(r,"_first")])
      distr <- c(distr, log[,paste0(r,"_first")])
      endpoint <- c(endpoint, rep("start", start_l))
      range <- c(range, rep(r, start_l))
      
      id <- which(priors$species==paste0(r,"_first"))
      un <- runif(5000, min=priors$end[id], max=priors$start[id])
      start_r_iqr <- c(start_r_iqr, IQR(log[,paste0(r,"_first")]))
      start_r_iqr_prior <- c(start_r_iqr_prior, IQR(un))
      
      range_s <- c(range_s, r)
      median_s<- c(median_s, median(log[,paste0(r,"_first")]))
      median_s_prior<- c(median_s_prior, median(un))
      
    }
    start_r <- c(start_r, median(log[,paste0(r,"_first")])) 
    end_r <- c(end_r, NA)
    hpd_f<- HPDinterval(as.mcmc(log[,paste0(r,"_first")]))
    start_r_low <- c(start_r_low, hpd_f[1])
    start_r_high <- c(start_r_high, hpd_f[2])
    start_r_min <- c(start_r_min, min(log[,paste0(r,"_first")]))
    start_r_max <- c(start_r_max, max(log[,paste0(r,"_first")]))
    
    
    
    # hpd_l<- HPDinterval(as.mcmc(log[,paste0(r,"_last")]))
    end_r_low<- c(end_r_low, NA)
    end_r_high<- c(end_r_high, NA)
    end_r_min <- c(end_r_min, NA)
    end_r_max <- c(end_r_max, NA)
    
  }
  
  for (r in ranges){
    if ((max(log[,paste0(r,"_first")])-min(log[,paste0(r,"_first")]))>0.1){
      start_l <- length(log[,paste0(r,"_first")])
      distr <- c(distr, log[,paste0(r,"_first")])
      endpoint <- c(endpoint, rep("start", start_l))
      range <- c(range, rep(r, start_l))
      
      id <- which(priors$species==paste0(r,"_first"))
      un <- runif(5000, min=priors$end[id], max=priors$start[id])
      start_r_iqr <- c(start_r_iqr, IQR(log[,paste0(r,"_first")]))
      start_r_iqr_prior <- c(start_r_iqr_prior, IQR(un))
      range_s <- c(range_s, r)
      median_s<- c(median_s, median(log[,paste0(r,"_first")]))
      median_s_prior<- c(median_s_prior, median(un))
      
    }
    if ((max(log[,paste0(r,"_last")])-min(log[,paste0(r,"_last")]))>0.1){
      end_l <- length(log[,paste0(r,"_last")])
      distr <- c(distr, log[,paste0(r,"_last")])
      endpoint <- c(endpoint, rep("end", end_l))
      range <- c(range, rep(r, end_l))
      
      id <- which(priors$species==paste0(r,"_last"))
      un <- runif(5000, min=priors$end[id], max=priors$start[id])
      end_r_iqr <- c(end_r_iqr, IQR(log[,paste0(r,"_last")]))
      end_r_iqr_prior <- c(end_r_iqr_prior, IQR(un))
      range_e <- c(range_e, r)
      median_e<- c(median_e, median(log[,paste0(r,"_last")]))
      median_e_prior<- c(median_e_prior, median(un))
    }
    
    start_r <- c(start_r, median(log[,paste0(r,"_first")]))
    end_r <- c(end_r, median(log[,paste0(r,"_last")]))
    
    hpd_f<- HPDinterval(as.mcmc(log[,paste0(r,"_first")]))
    start_r_low <- c(start_r_low, hpd_f[1])
    start_r_high <- c(start_r_high, hpd_f[2])
    start_r_min <- c(start_r_min, min(log[,paste0(r,"_first")]))
    start_r_max <- c(start_r_max, max(log[,paste0(r,"_first")]))
    
    
    hpd_l<- HPDinterval(as.mcmc(log[,paste0(r,"_last")]))
    end_r_low<- c(end_r_low, hpd_l[1])
    end_r_high<- c(end_r_high, hpd_l[2])
    end_r_min <- c(end_r_min, min(log[,paste0(r,"_last")]))
    end_r_max <- c(end_r_max, max(log[,paste0(r,"_last")]))
    
  }
  
  
  
  ranges_full_s <- data.frame(range=gsub("_", " ", range_s),
                              start_r_iqr=start_r_iqr,
                              start_r_iqr_prior=start_r_iqr_prior,
                              median_s_prior=median_s_prior,
                              median_s=median_s)
  ranges_full_e <- data.frame(range=gsub("_", " ", range_e),
                              end_r_iqr=end_r_iqr,
                              end_r_iqr_prior=end_r_iqr_prior,
                              median_e_prior=median_e_prior,
                              median_e=median_e)
  
  
  
  ranges_full <- data.frame(distr=distr, endpoint=endpoint, range=range)
  
  
  length(end_r)<-length(start_r)
  ranges_stat <- data.frame(range=c(ranges_single,ranges),
                            median_start=start_r,
                            median_end=end_r,
                            start_r_min=start_r_min,
                            start_r_low=start_r_low,
                            start_r_max=start_r_max,
                            start_r_high=start_r_high,
                            end_r_min=end_r_min,
                            end_r_low=end_r_low,
                            end_r_max=end_r_max,
                            end_r_high=end_r_high)
  
  return(list(ranges_full_s, ranges_full_e, ranges_full, ranges_stat))
  
}



log <- read.table(paste0(sRangesBothPath,"/penguins_inf_morph_at_both.tipAge.log"), header = T, )
log1 <- read.table(paste0(sRangesPath,"/penguins_inf_morph_at_start.tipAge.log"), header = T, )
name <- "combined"
all <- get_sRange_full(log1)
ranges_full_s <- all[[1]]
ranges_full_e <- all[[2]]

all <- get_sRange_full(log1)
ranges_full_s1 <- all[[1]]
ranges_full_e1 <- all[[2]]
rm(all)

ranges_full_s1$model <- "First"
ranges_full_e1$model <- "First"

ranges_full_s$model <- "Both"
ranges_full_e$model <- "Both"

ranges_full_s_comb <- rbind(ranges_full_s1, ranges_full_s)
ranges_full_e_comb <- rbind(ranges_full_e1, ranges_full_e)

p<-ggplot(ranges_full_s_comb, aes(x = range, y = start_r_iqr_prior/start_r_iqr)) +
  geom_point(aes(color = "First", shape=model), size = 3, alpha=0.7) + 
  geom_point(data=ranges_full_e_comb, 
             aes(x = range, y = end_r_iqr_prior/end_r_iqr, color = "Last", shape=model), 
             size = 3, alpha=0.7) +
  geom_hline(yintercept = 1)+ 
  ylab("IQR ratio (prior/posterior)") +
  # scale_y_continuous(breaks = seq(0, 12.5, by = 2.5), limits=c(0,12.7),
  #                    labels = as.character(seq(0, 12.5, by = 2.5))) +
  xlab("Species") + 
  scale_color_manual(values = c("First" = "#FF7F50", "Last" = "#43AA8B")) +
  theme_minimal()  + labs(color="Sample", shape="Morph. data attached to") +
  theme(#legend.position = c(0.80, 0.68), 
    legend.background = element_rect(fill="white", 
                                     linewidth=0.5, linetype="solid"),
    legend.justification = c("right", "bottom"),text=element_text(size=21),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ 
  guides(color = guide_legend(override.aes = list(shape=15, size=5) ) )

ggsave(paste0(figure13a_dir,"/penguin_SRFBD_",name,"_IQR.pdf"), width=14, height=7)
