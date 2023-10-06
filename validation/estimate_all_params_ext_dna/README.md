We use the following parameter priors for simualtion and inference:

- diversification rate ~Uniform(0.7, 0.9)
- turnover ~Uniform(0.2, 0.8)
- sampling proportion ~Uniform(0.2, 0.8)
- sampling at present ~Uniform(0.7, 1)
- origin was fixed at 4 in simulations but estimated with a Uniform(1, 1000) prior.


We rejected simulations with trees that had less than 5 or more than 1000 nodes. The simulation was done as follows:

1. Simulate 200 birth-death-sampling trees and associated taxonomy with FossilSim package.
2. For species with more than 1 sampled occurence, we record the first and last occurence.
   These are start and end of stratigraphic range, associated with that species.
3. For each occurence, we simulate morphological data with 7 states according to LewisMK substitution model.
   Strict clock value was fixed at 0.01.
4. Only for EXTANT occurences, we simulate DNA data with JukesCantor substitution model. Strict clock value was fixed at 0.01.


Then we infer the trees and associated parameters. Results are sumarrised in figures folder. 
The HPD_test.csv is a table of 95% HPD coverage for parameters. The qq plots for every parameter show HPD coverage at different percentages.
Ideally, you want z% coverage for z% HPD width (in other words, values in these plots should fall at x=y diagonal).
