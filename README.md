# s-vcom

This repository contains code to simluate trajectories from the discretized and stochasticified VCOM (S-VCOM) model. For infomation about the original deterministic ODE model, please see [our publication](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0187680).

The file "S_VCOM.pdf" contains a description of the mathematical approach used for this simulation code. The script "example.R" runs some examples and produces output.

## sim-src
  * mosquito-both.cpp: contains R-compatible C++14 code to run the stochastic and deterministic S-VCOM model, including ITN/IRS interventions as simple switches.
  * mosquito-deterministic.R: R code to simulate the deterministic S-VCOM model.
  * mosquito-equilibrium.R: R code to calculate discrete time equilibrium values.
  * mosquito-stochastic.R: R code to simulate the stochastic S-VCOM model.
  
