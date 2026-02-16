#!/bin/bash

# TreeStat2 needs to be installed!
# Run an R script to combine logs
Rscript combine.R

# Run R scripts to produce canid figures

echo "##########################################"
echo "##########################################"
echo "##########################################"
echo "START Figure5b-A14_Table2_directAncestors.R"
Rscript Figure5b-A14_Table2_directAncestors.R

echo "##########################################"
echo "##########################################"
echo "##########################################"
echo "START Figure7cd_clade_comparator.R"
Rscript Figure7cd_clade_comparator.R

echo "##########################################"
echo "##########################################"
echo "##########################################"
echo "START Figure8_TableA7.R"
Rscript Figure8_TableA7.R

echo "##########################################"
echo "##########################################"
echo "##########################################"
echo "START Figure9_tree_plotting_script.R"
Rscript Figure9_tree_plotting_script.R

echo "##########################################"
echo "##########################################"
echo "##########################################"
echo "START FigureA13b_parameter_comparison.R"
Rscript FigureA13b_parameter_comparison.R

echo "##########################################"
echo "##########################################"
echo "##########################################"
echo "START FigureA14_sampled_ancestors.R"
Rscript FigureA14_sampled_ancestors.R

echo "##########################################"
echo "##########################################"
echo "##########################################"
echo "START Rscript FigureA17_A18_posterior_clade_probs.R"
Rscript Rscript FigureA17_A18_posterior_clade_probs.R

echo "##########################################"
echo "All canid figures completed!"
echo "##########################################"
