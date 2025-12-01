### Simulation study description

We use the following parameter priors for simulation and inference:
  - diversification rate ~Uniform(0.7, 0.9)
  - turnover ~Uniform(0.2, 0.8)
  - sampling proportion ~Uniform(0.2, 0.8)
  - sampling at present ~Uniform(0.7, 1)
  - origin was fixed at 4 in simulations but estimated with a Uniform(1, 1000) prior.


We rejected simulations with trees that had less than 5 extant leaves or more than 1000 internal nodes. The simulation was done as follows:
1. Simulate 200 birth-death-sampling trees and associated taxonomy with FossilSim package.
2. For species with more than 1 sampled occurrence, we record the first and last occurrence. 
  These are start and end of stratigraphic range, associated with that species.
3. For each occurrence, we simulate morphological data with 7 states according to LewisMK substitution model. 
  Strict clock value was fixed at 0.01.
4. For EACH occurrence we simulate DNA data with JukesCantor substitution model. Strict clock value was fixed at 0.01.


Then we infer the trees and associated parameters. Results are summarised in figures folder. 
The HPD_test.csv is a table of 95% HPD coverage for parameters. The qq plots for every parameter show HPD coverage at different percentage.
Ideally, you want z% coverage for z% HPD width (in other words, values in these plots should fall at x=y diagonal).

### Reproducing the analysis

To reproduce our analysis, you should have R (https://www.r-project.org) and Beast2 (http://www.beast2.org) installed. This is enough to reproduce simulated trees and sequences. In order to run the inference, you should install the sRanges package in BEAST2.7 (see Readme.md at https://github.com/jugne/stratigraphic-ranges or linked Zenodo repository).

First, run the `simulate.R` script. This should create 200 folders with `sim` and `inf` subfolders. Folder `sim` stores simulation xmls, with trees and simulated morphological and DNA sequences. Folder `inf` has the inference xml for sRanges package in BEAST2.7.  The script will also create the `true_rates.csv` file that contains the parameter values under which true tree was simulated and the true sampled tree. 

Then, you can run the `inf/inference.xml` files with Beast2.7 for all 200 simulations. The results can then be summarised using the `FigureA10_validationPlots.R` script.
