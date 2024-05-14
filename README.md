# Trends-Model

This repository contains 5 files for the paper 'Detection of Temporal Trends in the Short-Term Health Effects of Air Pollution' by Rai, Brown, Morariu, and Shin. 

1. Simulation-Study-RW2-FE.Rmd: This script produces the simulation study results (plots and tables), except for the coverage probabilities in Table 1. 

2. Simulation-Study-Coverage.R: This script produces the coverage (percentages) in Table 1.

3. Simulation-Study-Helper-Functions.R: This script contains helper functions used by Simulation-Study-RW2-FE.Rmd and Simulation-Study-Coverage.R.

4. Momentum.cpp: This script contains a function used to calculate the momentum statistic.

5. sim_data_O3_Temp.RDS: A data file containing the temperature and air pollution data used for the simulation study. The data sources are noted in the paper.

These files assume that all 5 files are in R's working directory. It also assumes a 'Figures' folder exists in the working directory, and will save figures to that folder.
