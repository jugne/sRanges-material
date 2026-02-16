rm(list = ls())
gc()

source("penguin_helper.R")

figure_path <- paste0(figureA13a_dir, "/")

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

# Parse sampling dates using shared function
parsed_data <- parse_sampling_dates(text)
priors <- parsed_data$priors
ranges_single <- parsed_data$ranges_single
ranges <- parsed_data$ranges
ranges_ <- parsed_data$ranges_
sp <- priors$species
start <- priors$start
end <- priors$end

# Create species dataframe for CSV output
species_df <- create_species_dataframe(ranges_single, ranges, sp, start, end, ranges_)
write.csv(species_df, paste0(figure_path, "penguin_range_dates.csv"), row.names = FALSE)

# Load log data
log_both <- read.table(paste0(sRangesBothPath, "/tipAgeRelog/penguins_inf_morph_at_both.tipAge.log"), header = TRUE)
log_first <- read.table(paste0(sRangesPath, "/tipAgeRelog/penguins_inf_morph_at_start.tipAge.log"), header = TRUE)

# Process sRange workflow using shared function
result <- process_srange_workflow(log_both, log_first, priors, ranges_single, ranges)
ranges_full_s_comb <- result$ranges_full_s_comb
ranges_full_e_comb <- result$ranges_full_e_comb

# Create plot using shared function without custom y-axis scaling
p <- create_iqr_plot(
  ranges_full_s_comb,
  ranges_full_e_comb,
  y_scale_params = NULL,  # No custom y-axis parameters for penguins
  show_legend = TRUE,
  text_size = 21
)

ggsave(paste0(figure_path, "penguin_SRFBD_combined_IQR.pdf"), width=14, height=7)
