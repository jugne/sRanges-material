#!/bin/bash

# TreeStat2 needs to be installed!
# Run an R script to combine logs
Rscript combine.R

# Run R scripts to produce figures
echo "##########################################"
echo "##########################################"
echo "##########################################"
echo "START Figure4_TablesA6.R"
Rscript Figure4_TablesA6.R

echo "##########################################"
echo "##########################################"
echo "##########################################"
echo "START FigureA13a_ranges_iqr.R"
Rscript FigureA13a_ranges_iqr.R

echo "##########################################"
echo "##########################################"
echo "##########################################"
echo "START Figure5a_Table2_directAncestors.R"
Rscript Figure5a_Table2_directAncestors.R

echo "##########################################"
echo "##########################################"
echo "##########################################"
echo "START Figure7ab_clade_comparator.R"
Rscript Figure7ab_clade_comparator.R

echo "##########################################"
echo "##########################################"
echo "##########################################"
echo "START Tables3_4_get_crown_data.R"
Rscript Tables3_4_get_crown_data.R

echo "##########################################"
echo "##########################################"
echo "##########################################"
echo "START TableA5_missingChars.R"
Rscript TableA5_missingChars.R

echo "##########################################"
echo "##########################################"
echo "##########################################"
echo "START FigureA15_A16_posterior_clade_probs.R"
Rscript FigureA15_A16_posterior_clade_probs.R

echo "##########################################"
echo "##########################################"
echo "##########################################"
echo "START Figure6_tree_plotting_script.R"
Rscript Figure6_tree_plotting_script.R
