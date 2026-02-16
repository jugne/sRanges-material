rm(list = ls())
gc()

source("canid_helper.R")

figure_path <- paste0(figureA13b_dir, "/")



# Define the text containing the sampling dates information
text <- '<samplingDates id="samplingDate0" spec="sa.evolution.tree.SamplingDate" taxon="@Aelurodon_asthenostylus_first" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate1" spec="sa.evolution.tree.SamplingDate" taxon="@Aelurodon_ferox_first" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate2" spec="sa.evolution.tree.SamplingDate" taxon="@Aelurodon_mcgrewi_first" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate3" spec="sa.evolution.tree.SamplingDate" taxon="@Aelurodon_stirtoni_first" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate4" spec="sa.evolution.tree.SamplingDate" taxon="@Aelurodon_taxoides_first" lower="13.6" upper="23.03"/>
<samplingDates id="samplingDate5" spec="sa.evolution.tree.SamplingDate" taxon="@Aenocyon_dirus_first" lower="0.781" upper="3.6"/>
<samplingDates id="samplingDate6" spec="sa.evolution.tree.SamplingDate" taxon="@Archaeocyon_leptodus_first" lower="30.8" upper="33.3"/>
<samplingDates id="samplingDate7" spec="sa.evolution.tree.SamplingDate" taxon="@Archaeocyon_pavidus_first" lower="30.8" upper="33.3"/>
<samplingDates id="samplingDate8" spec="sa.evolution.tree.SamplingDate" taxon="@Borophagus_diversidens_first" lower="1.806" upper="4.9"/>
<samplingDates id="samplingDate9" spec="sa.evolution.tree.SamplingDate" taxon="@Borophagus_dudleyi_first" lower="4.9" upper="10.3"/>
<samplingDates id="samplingDate10" spec="sa.evolution.tree.SamplingDate" taxon="@Borophagus_hilli_first" lower="4.9" upper="10.3"/>
<samplingDates id="samplingDate11" spec="sa.evolution.tree.SamplingDate" taxon="@Borophagus_littoralis_first" lower="10.3" upper="13.6"/>
<samplingDates id="samplingDate12" spec="sa.evolution.tree.SamplingDate" taxon="@Borophagus_orc_first" lower="4.9" upper="10.3"/>
<samplingDates id="samplingDate13" spec="sa.evolution.tree.SamplingDate" taxon="@Borophagus_parvus_first" lower="4.9" upper="10.3"/>
<samplingDates id="samplingDate14" spec="sa.evolution.tree.SamplingDate" taxon="@Borophagus_pugnator_first" lower="10.3" upper="13.6"/>
<samplingDates id="samplingDate15" spec="sa.evolution.tree.SamplingDate" taxon="@Borophagus_secundus_first" lower="10.3" upper="13.6"/>
<samplingDates id="samplingDate16" spec="sa.evolution.tree.SamplingDate" taxon="@Canis_armbrusteri_first" lower="0.781" upper="1.806"/>
<samplingDates id="samplingDate17" spec="sa.evolution.tree.SamplingDate" taxon="@Canis_edwardii_first" lower="1.806" upper="4.9"/>
<samplingDates id="samplingDate18" spec="sa.evolution.tree.SamplingDate" taxon="@Canis_latrans_first" lower="1.806" upper="4.9"/>
<samplingDates id="samplingDate19" spec="sa.evolution.tree.SamplingDate" taxon="@Canis_lepophagus_first" lower="1.8" upper="4.9"/>
<samplingDates id="samplingDate20" spec="sa.evolution.tree.SamplingDate" taxon="@Canis_lupus_first" lower="2.588" upper="3.6"/>
<samplingDates id="samplingDate21" spec="sa.evolution.tree.SamplingDate" taxon="@Carpocyon_compressus_first" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate22" spec="sa.evolution.tree.SamplingDate" taxon="@Carpocyon_limosus_first" lower="4.9" upper="10.3"/>
<samplingDates id="samplingDate23" spec="sa.evolution.tree.SamplingDate" taxon="@Carpocyon_robustus_first" lower="10.3" upper="13.6"/>
<samplingDates id="samplingDate24" spec="sa.evolution.tree.SamplingDate" taxon="@Carpocyon_webbi_first" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate25" spec="sa.evolution.tree.SamplingDate" taxon="@Cerdocyon_texanus_first" lower="4.9" upper="10.3"/>
<samplingDates id="samplingDate26" spec="sa.evolution.tree.SamplingDate" taxon="@Cormocyon_copei_first" lower="26.3" upper="30.8"/>
<samplingDates id="samplingDate27" spec="sa.evolution.tree.SamplingDate" taxon="@Cormocyon_haydeni_first" lower="20.43" upper="30.8"/>
<samplingDates id="samplingDate28" spec="sa.evolution.tree.SamplingDate" taxon="@Cuon_alpinus_first" lower="0.3" upper="2.588"/>
<samplingDates id="samplingDate29" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctoides_acridens_first" lower="20.43" upper="30.8"/>
<samplingDates id="samplingDate30" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctoides_emryi_first" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate31" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctoides_gawnae_first" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate32" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctoides_lemur_first" lower="26.3" upper="30.8"/>
<samplingDates id="samplingDate33" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctoides_luskensis_first" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate34" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctoides_roii_first" lower="26.3" upper="30.8"/>
<samplingDates id="samplingDate35" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctus_crucidens_first" lower="11.62" upper="13.82"/>
<samplingDates id="samplingDate36" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctus_galushai_first" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate37" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctus_marylandica_first" lower="13.82" upper="15.97"/>
<samplingDates id="samplingDate38" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctus_saxatilis_first" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate39" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctus_voorhiesi_first" lower="10.3" upper="13.6"/>
<samplingDates id="samplingDate40" spec="sa.evolution.tree.SamplingDate" taxon="@Cynodesmus_thooides_first" lower="30.8" upper="33.3"/>
<samplingDates id="samplingDate41" spec="sa.evolution.tree.SamplingDate" taxon="@Desmocyon_matthewi_first" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate42" spec="sa.evolution.tree.SamplingDate" taxon="@Desmocyon_thomsoni_first" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate43" spec="sa.evolution.tree.SamplingDate" taxon="@Ectopocynus_antiquus_first" lower="30.8" upper="33.3"/>
<samplingDates id="samplingDate44" spec="sa.evolution.tree.SamplingDate" taxon="@Ectopocynus_simplicidens_first" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate45" spec="sa.evolution.tree.SamplingDate" taxon="@Enhydrocyon_basilatus_first" lower="20.43" upper="30.8"/>
<samplingDates id="samplingDate46" spec="sa.evolution.tree.SamplingDate" taxon="@Enhydrocyon_crassidens_first" lower="26.3" upper="30.8"/>
<samplingDates id="samplingDate47" spec="sa.evolution.tree.SamplingDate" taxon="@Enhydrocyon_pahinsintewakpa_first" lower="26.3" upper="30.8"/>
<samplingDates id="samplingDate48" spec="sa.evolution.tree.SamplingDate" taxon="@Enhydrocyon_stenocephalus_first" lower="26.3" upper="30.8"/>
<samplingDates id="samplingDate49" spec="sa.evolution.tree.SamplingDate" taxon="@Epicyon_haydeni_first" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate50" spec="sa.evolution.tree.SamplingDate" taxon="@Epicyon_saevus_first" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate51" spec="sa.evolution.tree.SamplingDate" taxon="@Eucyon_davisi_first" lower="5.333" upper="11.608"/>
<samplingDates id="samplingDate52" spec="sa.evolution.tree.SamplingDate" taxon="@Eucyon_ferox_first" lower="4.9" upper="10.3"/>
<samplingDates id="samplingDate53" spec="sa.evolution.tree.SamplingDate" taxon="@Euoplocyon_brachygnathus_first" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate54" spec="sa.evolution.tree.SamplingDate" taxon="@Hesperocyon_coloradensis_first" lower="33.3" upper="33.9"/>
<samplingDates id="samplingDate55" spec="sa.evolution.tree.SamplingDate" taxon="@Hesperocyon_gregarius_first" lower="33.9" upper="37.2"/>
<samplingDates id="samplingDate56" spec="sa.evolution.tree.SamplingDate" taxon="@Leptocyon_douglassi_first" lower="26.3" upper="30.8"/>
<samplingDates id="samplingDate57" spec="sa.evolution.tree.SamplingDate" taxon="@Leptocyon_gregorii_first" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate58" spec="sa.evolution.tree.SamplingDate" taxon="@Leptocyon_leidyi_first" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate59" spec="sa.evolution.tree.SamplingDate" taxon="@Leptocyon_matthewi_first" lower="10.3" upper="13.6"/>
<samplingDates id="samplingDate60" spec="sa.evolution.tree.SamplingDate" taxon="@Leptocyon_vafer_first" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate61" spec="sa.evolution.tree.SamplingDate" taxon="@Leptocyon_vulpinus_first" lower="20.43" upper="30.8"/>
<samplingDates id="samplingDate62" spec="sa.evolution.tree.SamplingDate" taxon="@Lycaon_pictus_first" lower="0.0117" upper="2.588"/>
<samplingDates id="samplingDate63" spec="sa.evolution.tree.SamplingDate" taxon="@Mesocyon_brachyops_first" lower="26.3" upper="30.8"/>
<samplingDates id="samplingDate64" spec="sa.evolution.tree.SamplingDate" taxon="@Mesocyon_coryphaeus_first" lower="26.3" upper="30.8"/>
<samplingDates id="samplingDate65" spec="sa.evolution.tree.SamplingDate" taxon="@Mesocyon_temnodon_first" lower="33.3" upper="33.9"/>
<samplingDates id="samplingDate66" spec="sa.evolution.tree.SamplingDate" taxon="@Metalopex_macconnelli_first" lower="10.3" upper="13.6"/>
<samplingDates id="samplingDate67" spec="sa.evolution.tree.SamplingDate" taxon="@Metalopex_merriami_first" lower="10.3" upper="13.6"/>
<samplingDates id="samplingDate68" spec="sa.evolution.tree.SamplingDate" taxon="@Metatomarctus_canavus_first" lower="15.97" upper="20.44"/>
<samplingDates id="samplingDate69" spec="sa.evolution.tree.SamplingDate" taxon="@Microtomarctus_conferta_first" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate70" spec="sa.evolution.tree.SamplingDate" taxon="@Osbornodon_brachypus_first" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate71" spec="sa.evolution.tree.SamplingDate" taxon="@Osbornodon_fricki_first" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate72" spec="sa.evolution.tree.SamplingDate" taxon="@Osbornodon_iamonensis_first" lower="15.97" upper="23.03"/>
<samplingDates id="samplingDate73" spec="sa.evolution.tree.SamplingDate" taxon="@Osbornodon_renjiei_first" lower="33.3" upper="33.9"/>
<samplingDates id="samplingDate74" spec="sa.evolution.tree.SamplingDate" taxon="@Osbornodon_sesnoni_first" lower="30.8" upper="33.3"/>
<samplingDates id="samplingDate75" spec="sa.evolution.tree.SamplingDate" taxon="@Otarocyon_cooki_first" lower="26.3" upper="30.8"/>
<samplingDates id="samplingDate76" spec="sa.evolution.tree.SamplingDate" taxon="@Otarocyon_macdonaldi_first" lower="33.3" upper="37.2"/>
<samplingDates id="samplingDate77" spec="sa.evolution.tree.SamplingDate" taxon="@Oxetocyon_cuspidatus_first" lower="30.8" upper="33.3"/>
<samplingDates id="samplingDate78" spec="sa.evolution.tree.SamplingDate" taxon="@Paracynarctus_kelloggi_first" lower="15.97" upper="23.03"/>
<samplingDates id="samplingDate79" spec="sa.evolution.tree.SamplingDate" taxon="@Paracynarctus_sinclairi_first" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate80" spec="sa.evolution.tree.SamplingDate" taxon="@Paraenhydrocyon_josephi_first" lower="26.3" upper="30.8"/>
<samplingDates id="samplingDate81" spec="sa.evolution.tree.SamplingDate" taxon="@Paraenhydrocyon_robustus_first" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate82" spec="sa.evolution.tree.SamplingDate" taxon="@Paraenhydrocyon_wallovianus_first" lower="23.03" upper="28.1"/>
<samplingDates id="samplingDate83" spec="sa.evolution.tree.SamplingDate" taxon="@Paratomarctus_euthos_first" lower="10.3" upper="23.03"/>
<samplingDates id="samplingDate84" spec="sa.evolution.tree.SamplingDate" taxon="@Paratomarctus_temerarius_first" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate85" spec="sa.evolution.tree.SamplingDate" taxon="@Philotrox_condoni_first" lower="26.3" upper="30.8"/>
<samplingDates id="samplingDate86" spec="sa.evolution.tree.SamplingDate" taxon="@Phlaocyon_annectens_first" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate87" spec="sa.evolution.tree.SamplingDate" taxon="@Phlaocyon_latidens_first" lower="30.8" upper="33.3"/>
<samplingDates id="samplingDate88" spec="sa.evolution.tree.SamplingDate" taxon="@Phlaocyon_leucosteus_first" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate89" spec="sa.evolution.tree.SamplingDate" taxon="@Phlaocyon_marslandensis_first" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate90" spec="sa.evolution.tree.SamplingDate" taxon="@Phlaocyon_minor_first" lower="24.8" upper="30.8"/>
<samplingDates id="samplingDate91" spec="sa.evolution.tree.SamplingDate" taxon="@Protepicyon_raki_first" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate92" spec="sa.evolution.tree.SamplingDate" taxon="@Protomarctus_optatus_first" lower="15.97" upper="20.44"/>
<samplingDates id="samplingDate93" spec="sa.evolution.tree.SamplingDate" taxon="@Psalidocyon_marianae_first" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate94" spec="sa.evolution.tree.SamplingDate" taxon="@Rhizocyon_oregonensis_first" lower="26.3" upper="30.8"/>
<samplingDates id="samplingDate95" spec="sa.evolution.tree.SamplingDate" taxon="@Sunkahetanka_geringensis_first" lower="26.3" upper="30.8"/>
<samplingDates id="samplingDate96" spec="sa.evolution.tree.SamplingDate" taxon="@Tephrocyon_rurestris_first" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate97" spec="sa.evolution.tree.SamplingDate" taxon="@Tomarctus_brevirostris_first" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate98" spec="sa.evolution.tree.SamplingDate" taxon="@Tomarctus_hippophaga_first" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate99" spec="sa.evolution.tree.SamplingDate" taxon="@Urocyon_cinereoargenteus_first" lower="1.806" upper="2.588"/>
<samplingDates id="samplingDate100" spec="sa.evolution.tree.SamplingDate" taxon="@Urocyon_galushai_first" lower="1.8" upper="4.9"/>
<samplingDates id="samplingDate101" spec="sa.evolution.tree.SamplingDate" taxon="@Urocyon_minicephalus_first" lower="0.3" upper="1.8"/>
<samplingDates id="samplingDate102" spec="sa.evolution.tree.SamplingDate" taxon="@Vulpes_stenognathus_first" lower="4.9" upper="13.6"/>
<samplingDates id="samplingDate103" spec="sa.evolution.tree.SamplingDate" taxon="@Xenocyon_falconeri_first" lower="1.806" upper="2.588"/>
<samplingDates id="samplingDate104" spec="sa.evolution.tree.SamplingDate" taxon="@Aelurodon_asthenostylus_last" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate105" spec="sa.evolution.tree.SamplingDate" taxon="@Aelurodon_ferox_last" lower="10.3" upper="13.6"/>
<samplingDates id="samplingDate106" spec="sa.evolution.tree.SamplingDate" taxon="@Aelurodon_mcgrewi_last" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate107" spec="sa.evolution.tree.SamplingDate" taxon="@Aelurodon_stirtoni_last" lower="10.3" upper="13.6"/>
<samplingDates id="samplingDate108" spec="sa.evolution.tree.SamplingDate" taxon="@Aelurodon_taxoides_last" lower="3.6" upper="10.3"/>
<samplingDates id="samplingDate109" spec="sa.evolution.tree.SamplingDate" taxon="@Aenocyon_dirus_last" lower="0.0117" upper="0.126"/>
<samplingDates id="samplingDate110" spec="sa.evolution.tree.SamplingDate" taxon="@Archaeocyon_leptodus_last" lower="20.43" upper="30.8"/>
<samplingDates id="samplingDate111" spec="sa.evolution.tree.SamplingDate" taxon="@Archaeocyon_pavidus_last" lower="23.03" upper="28.1"/>
<samplingDates id="samplingDate112" spec="sa.evolution.tree.SamplingDate" taxon="@Borophagus_diversidens_last" lower="1.8" upper="2.588"/>
<samplingDates id="samplingDate113" spec="sa.evolution.tree.SamplingDate" taxon="@Borophagus_dudleyi_last" lower="3.6" upper="5.333"/>
<samplingDates id="samplingDate114" spec="sa.evolution.tree.SamplingDate" taxon="@Borophagus_hilli_last" lower="1.8" upper="4.9"/>
<samplingDates id="samplingDate115" spec="sa.evolution.tree.SamplingDate" taxon="@Borophagus_littoralis_last" lower="5.333" upper="11.608"/>
<samplingDates id="samplingDate116" spec="sa.evolution.tree.SamplingDate" taxon="@Borophagus_orc_last" lower="3.6" upper="5.333"/>
<samplingDates id="samplingDate117" spec="sa.evolution.tree.SamplingDate" taxon="@Borophagus_parvus_last" lower="4.9" upper="10.3"/>
<samplingDates id="samplingDate118" spec="sa.evolution.tree.SamplingDate" taxon="@Borophagus_pugnator_last" lower="4.9" upper="10.3"/>
<samplingDates id="samplingDate119" spec="sa.evolution.tree.SamplingDate" taxon="@Borophagus_secundus_last" lower="4.9" upper="10.3"/>
<samplingDates id="samplingDate120" spec="sa.evolution.tree.SamplingDate" taxon="@Canis_armbrusteri_last" lower="0.0117" upper="0.126"/>
<samplingDates id="samplingDate121" spec="sa.evolution.tree.SamplingDate" taxon="@Canis_edwardii_last" lower="0.3" upper="1.8"/>
<samplingDates id="samplingDate122" spec="sa.evolution.tree.SamplingDate" taxon="@Canis_lepophagus_last" lower="0.781" upper="2.588"/>
<samplingDates id="samplingDate123" spec="sa.evolution.tree.SamplingDate" taxon="@Carpocyon_compressus_last" lower="10.3" upper="15.97"/>
<samplingDates id="samplingDate124" spec="sa.evolution.tree.SamplingDate" taxon="@Carpocyon_limosus_last" lower="4.9" upper="10.3"/>
<samplingDates id="samplingDate125" spec="sa.evolution.tree.SamplingDate" taxon="@Carpocyon_robustus_last" lower="10.3" upper="13.6"/>
<samplingDates id="samplingDate126" spec="sa.evolution.tree.SamplingDate" taxon="@Carpocyon_webbi_last" lower="4.9" upper="10.3"/>
<samplingDates id="samplingDate127" spec="sa.evolution.tree.SamplingDate" taxon="@Cerdocyon_texanus_last" lower="1.8" upper="4.9"/>
<samplingDates id="samplingDate128" spec="sa.evolution.tree.SamplingDate" taxon="@Cormocyon_copei_last" lower="11" upper="23.03"/>
<samplingDates id="samplingDate129" spec="sa.evolution.tree.SamplingDate" taxon="@Cormocyon_haydeni_last" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate130" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctoides_acridens_last" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate131" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctoides_emryi_last" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate132" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctoides_gawnae_last" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate133" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctoides_lemur_last" lower="20.43" upper="28.1"/>
<samplingDates id="samplingDate134" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctoides_luskensis_last" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate135" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctoides_roii_last" lower="26.3" upper="30.8"/>
<samplingDates id="samplingDate136" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctus_crucidens_last" lower="4.9" upper="11.608"/>
<samplingDates id="samplingDate137" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctus_galushai_last" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate138" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctus_marylandica_last" lower="13.82" upper="15.97"/>
<samplingDates id="samplingDate139" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctus_saxatilis_last" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate140" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctus_voorhiesi_last" lower="10.3" upper="13.6"/>
<samplingDates id="samplingDate141" spec="sa.evolution.tree.SamplingDate" taxon="@Cynodesmus_thooides_last" lower="20.43" upper="30.8"/>
<samplingDates id="samplingDate142" spec="sa.evolution.tree.SamplingDate" taxon="@Desmocyon_matthewi_last" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate143" spec="sa.evolution.tree.SamplingDate" taxon="@Desmocyon_thomsoni_last" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate144" spec="sa.evolution.tree.SamplingDate" taxon="@Ectopocynus_antiquus_last" lower="20.43" upper="30.8"/>
<samplingDates id="samplingDate145" spec="sa.evolution.tree.SamplingDate" taxon="@Ectopocynus_simplicidens_last" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate146" spec="sa.evolution.tree.SamplingDate" taxon="@Enhydrocyon_basilatus_last" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate147" spec="sa.evolution.tree.SamplingDate" taxon="@Enhydrocyon_crassidens_last" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate148" spec="sa.evolution.tree.SamplingDate" taxon="@Enhydrocyon_pahinsintewakpa_last" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate149" spec="sa.evolution.tree.SamplingDate" taxon="@Enhydrocyon_stenocephalus_last" lower="20.43" upper="30.8"/>
<samplingDates id="samplingDate150" spec="sa.evolution.tree.SamplingDate" taxon="@Epicyon_haydeni_last" lower="4.9" upper="10.3"/>
<samplingDates id="samplingDate151" spec="sa.evolution.tree.SamplingDate" taxon="@Epicyon_saevus_last" lower="4.9" upper="10.3"/>
<samplingDates id="samplingDate152" spec="sa.evolution.tree.SamplingDate" taxon="@Eucyon_davisi_last" lower="1.8" upper="3.6"/>
<samplingDates id="samplingDate153" spec="sa.evolution.tree.SamplingDate" taxon="@Eucyon_ferox_last" lower="1.8" upper="4.9"/>
<samplingDates id="samplingDate154" spec="sa.evolution.tree.SamplingDate" taxon="@Euoplocyon_brachygnathus_last" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate155" spec="sa.evolution.tree.SamplingDate" taxon="@Hesperocyon_coloradensis_last" lower="33.3" upper="33.9"/>
<samplingDates id="samplingDate156" spec="sa.evolution.tree.SamplingDate" taxon="@Hesperocyon_gregarius_last" lower="20.43" upper="30.8"/>
<samplingDates id="samplingDate157" spec="sa.evolution.tree.SamplingDate" taxon="@Leptocyon_douglassi_last" lower="20.43" upper="30.8"/>
<samplingDates id="samplingDate158" spec="sa.evolution.tree.SamplingDate" taxon="@Leptocyon_gregorii_last" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate159" spec="sa.evolution.tree.SamplingDate" taxon="@Leptocyon_leidyi_last" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate160" spec="sa.evolution.tree.SamplingDate" taxon="@Leptocyon_matthewi_last" lower="4.9" upper="13.6"/>
<samplingDates id="samplingDate161" spec="sa.evolution.tree.SamplingDate" taxon="@Leptocyon_vafer_last" lower="10.3" upper="13.6"/>
<samplingDates id="samplingDate162" spec="sa.evolution.tree.SamplingDate" taxon="@Leptocyon_vulpinus_last" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate163" spec="sa.evolution.tree.SamplingDate" taxon="@Mesocyon_brachyops_last" lower="23.03" upper="28.1"/>
<samplingDates id="samplingDate164" spec="sa.evolution.tree.SamplingDate" taxon="@Mesocyon_coryphaeus_last" lower="20.43" upper="30.8"/>
<samplingDates id="samplingDate165" spec="sa.evolution.tree.SamplingDate" taxon="@Mesocyon_temnodon_last" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate166" spec="sa.evolution.tree.SamplingDate" taxon="@Metalopex_macconnelli_last" lower="4.9" upper="10.3"/>
<samplingDates id="samplingDate167" spec="sa.evolution.tree.SamplingDate" taxon="@Metalopex_merriami_last" lower="4.9" upper="10.3"/>
<samplingDates id="samplingDate168" spec="sa.evolution.tree.SamplingDate" taxon="@Metatomarctus_canavus_last" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate169" spec="sa.evolution.tree.SamplingDate" taxon="@Microtomarctus_conferta_last" lower="10.3" upper="15.97"/>
<samplingDates id="samplingDate170" spec="sa.evolution.tree.SamplingDate" taxon="@Osbornodon_brachypus_last" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate171" spec="sa.evolution.tree.SamplingDate" taxon="@Osbornodon_fricki_last" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate172" spec="sa.evolution.tree.SamplingDate" taxon="@Osbornodon_iamonensis_last" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate173" spec="sa.evolution.tree.SamplingDate" taxon="@Osbornodon_renjiei_last" lower="30.8" upper="33.3"/>
<samplingDates id="samplingDate174" spec="sa.evolution.tree.SamplingDate" taxon="@Osbornodon_sesnoni_last" lower="26.3" upper="30.8"/>
<samplingDates id="samplingDate175" spec="sa.evolution.tree.SamplingDate" taxon="@Otarocyon_cooki_last" lower="20.43" upper="30.8"/>
<samplingDates id="samplingDate176" spec="sa.evolution.tree.SamplingDate" taxon="@Otarocyon_macdonaldi_last" lower="33.3" upper="33.9"/>
<samplingDates id="samplingDate177" spec="sa.evolution.tree.SamplingDate" taxon="@Oxetocyon_cuspidatus_last" lower="30.8" upper="33.3"/>
<samplingDates id="samplingDate178" spec="sa.evolution.tree.SamplingDate" taxon="@Paracynarctus_kelloggi_last" lower="5.333" upper="15.97"/>
<samplingDates id="samplingDate179" spec="sa.evolution.tree.SamplingDate" taxon="@Paracynarctus_sinclairi_last" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate180" spec="sa.evolution.tree.SamplingDate" taxon="@Paraenhydrocyon_josephi_last" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate181" spec="sa.evolution.tree.SamplingDate" taxon="@Paraenhydrocyon_robustus_last" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate182" spec="sa.evolution.tree.SamplingDate" taxon="@Paraenhydrocyon_wallovianus_last" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate183" spec="sa.evolution.tree.SamplingDate" taxon="@Paratomarctus_euthos_last" lower="4.9" upper="13.6"/>
<samplingDates id="samplingDate184" spec="sa.evolution.tree.SamplingDate" taxon="@Paratomarctus_temerarius_last" lower="4.9" upper="13.6"/>
<samplingDates id="samplingDate185" spec="sa.evolution.tree.SamplingDate" taxon="@Philotrox_condoni_last" lower="20.43" upper="30.8"/>
<samplingDates id="samplingDate186" spec="sa.evolution.tree.SamplingDate" taxon="@Phlaocyon_annectens_last" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate187" spec="sa.evolution.tree.SamplingDate" taxon="@Phlaocyon_latidens_last" lower="26.3" upper="30.8"/>
<samplingDates id="samplingDate188" spec="sa.evolution.tree.SamplingDate" taxon="@Phlaocyon_leucosteus_last" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate189" spec="sa.evolution.tree.SamplingDate" taxon="@Phlaocyon_marslandensis_last" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate190" spec="sa.evolution.tree.SamplingDate" taxon="@Phlaocyon_minor_last" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate191" spec="sa.evolution.tree.SamplingDate" taxon="@Protepicyon_raki_last" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate192" spec="sa.evolution.tree.SamplingDate" taxon="@Protomarctus_optatus_last" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate193" spec="sa.evolution.tree.SamplingDate" taxon="@Psalidocyon_marianae_last" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate194" spec="sa.evolution.tree.SamplingDate" taxon="@Rhizocyon_oregonensis_last" lower="26.3" upper="30.8"/>
<samplingDates id="samplingDate195" spec="sa.evolution.tree.SamplingDate" taxon="@Sunkahetanka_geringensis_last" lower="20.43" upper="30.8"/>
<samplingDates id="samplingDate196" spec="sa.evolution.tree.SamplingDate" taxon="@Tephrocyon_rurestris_last" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate197" spec="sa.evolution.tree.SamplingDate" taxon="@Tomarctus_brevirostris_last" lower="11.608" upper="15.97"/>
<samplingDates id="samplingDate198" spec="sa.evolution.tree.SamplingDate" taxon="@Tomarctus_hippophaga_last" lower="13.6" upper="15.97"/>
<samplingDates id="samplingDate199" spec="sa.evolution.tree.SamplingDate" taxon="@Urocyon_galushai_last" lower="1.8" upper="4.9"/>
<samplingDates id="samplingDate200" spec="sa.evolution.tree.SamplingDate" taxon="@Urocyon_minicephalus_last" lower="0.3" upper="1.8"/>
<samplingDates id="samplingDate201" spec="sa.evolution.tree.SamplingDate" taxon="@Vulpes_stenognathus_last" lower="0.3" upper="1.8"/>
<samplingDates id="samplingDate202" spec="sa.evolution.tree.SamplingDate" taxon="@Xenocyon_falconeri_last" lower="0.0117" upper="2.588"/>
<samplingDates id="samplingDate203" spec="sa.evolution.tree.SamplingDate" taxon="@Archaeocyon_falkenbachi_first" lower="20.43" upper="30.8"/>
<samplingDates id="samplingDate204" spec="sa.evolution.tree.SamplingDate" taxon="@Caedocyon_tedfordi_first" lower="20.43" upper="30.8"/>
<samplingDates id="samplingDate205" spec="sa.evolution.tree.SamplingDate" taxon="@Canis_thooides_first" lower="1.8" upper="4.9"/>
<samplingDates id="samplingDate206" spec="sa.evolution.tree.SamplingDate" taxon="@Cynarctoides_harlowi_first" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate207" spec="sa.evolution.tree.SamplingDate" taxon="@Cynodesmus_martini_first" lower="20.43" upper="30.8"/>
<samplingDates id="samplingDate208" spec="sa.evolution.tree.SamplingDate" taxon="@Ectopocynus_intermedius_first" lower="20.43" upper="30.8"/>
<samplingDates id="samplingDate209" spec="sa.evolution.tree.SamplingDate" taxon="@Epicyon_aelurodontoides_first" lower="4.9" upper="10.3"/>
<samplingDates id="samplingDate210" spec="sa.evolution.tree.SamplingDate" taxon="@Euoplocyon_spissidens_first" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate211" spec="sa.evolution.tree.SamplingDate" taxon="@Leptocyon_mollis_first" lower="26.3" upper="30.8"/>
<samplingDates id="samplingDate212" spec="sa.evolution.tree.SamplingDate" taxon="@Phlaocyon_achoros_first" lower="23.03" upper="25"/>
<samplingDates id="samplingDate213" spec="sa.evolution.tree.SamplingDate" taxon="@Phlaocyon_mariae_first" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate214" spec="sa.evolution.tree.SamplingDate" taxon="@Phlaocyon_multicuspus_first" lower="20.43" upper="24.8"/>
<samplingDates id="samplingDate215" spec="sa.evolution.tree.SamplingDate" taxon="@Phlaocyon_yatkolai_first" lower="15.97" upper="20.43"/>
<samplingDates id="samplingDate216" spec="sa.evolution.tree.SamplingDate" taxon="@Prohesperocyon_wilsoni_first" lower="33.9" upper="37.2"/>
<samplingDates id="samplingDate217" spec="sa.evolution.tree.SamplingDate" taxon="@Urocyon_citrinus_first" lower="1.806" upper="2.588"/>
<samplingDates id="samplingDate218" spec="sa.evolution.tree.SamplingDate" taxon="@Urocyon_webbi_first" lower="4.9" upper="10.3"/>
<samplingDates id="samplingDate219" spec="sa.evolution.tree.SamplingDate" taxon="@Xenocyon_texanus_first" lower="0.63" upper="1.8"/>
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
write.csv(species_df, paste0(figure_path, "canid_range_dates.csv"), row.names = FALSE)

# Load log data
log_both <- read.table(paste0(sRangesBothPath, "/tipAgeRelog/canidae_rho.tipAge.log"), header = TRUE)
log_first <- read.table(paste0(sRangesPath, "/tipAgeRelog/canidae_rho.tipAge.log"), header = TRUE)

# Process sRange workflow using shared function
result <- process_srange_workflow(log_both, log_first, priors, ranges_single, ranges)
ranges_full_s_comb <- result$ranges_full_s_comb
ranges_full_e_comb <- result$ranges_full_e_comb

# Create plot using shared function with custom y-axis scaling
y_scale_params <- list(
  breaks = seq(0, 12.5, by = 2.5),
  limits = c(0, 12.7),
  labels = as.character(seq(0, 12.5, by = 2.5))
)

p <- create_iqr_plot(
  ranges_full_s_comb,
  ranges_full_e_comb,
  y_scale_params = y_scale_params,
  show_legend = FALSE,
  text_size = 21
)

ggsave(paste0(figure_path, "canid_SRFBD_combined_IQR_noLegend.pdf"), width=25, height=8)
