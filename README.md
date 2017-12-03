## Overview

This repository contains MATLAB, Python and R implementations of the Marginal Path Likelihood (MPL) algorithm, which infers the selection coefficients from observed mutant frequencies, assuming these frequences evolved according to a Wright Fisher (WF) process. For each implementation, a function generating the WF mutant frequency trajectories (`WF_sim_traj.m/R/py`) and estimating the selection coefficients using MPL (`estimate_MPL.m/R/py`) is provided.

Also included is an example script (`main_WF_MPL.m/R/py`) which demonstrates how these two functions are run, and generates visual plots to demonstrate the performance of the MPL algorithm, using ggplot2 (Python) and matplotlib (R).

For each of the MATLAB, Python and R implementations, the two functions implementing the WF process and MPL algorithm have the same inputs and outputs, described as follows:

## `WF_sim_traj`

 This function generates the WF single and double mutant frequency trajectories, and described by

` WF_sim_traj(s,mu,L,N,p_init,dt_array)`

The inputs are:

`s` : selection coefficients 

`mu`: mutation probability

`L` : number of loci

`N` : population size

`p_init` : initial genotype frequencies

`dt_array` : an array containing the generation number when the frequencies are observed

The outputs are:

`single_mut`: single mutant frequencies, stored in a T x L matrix where T is the number of observed generations

`double_mut`: double mutant frequencies, stored in a T x L x L matrix 


## `estimate_MPL`

The function estimates the selection coefficients from the single and double mutant frequency WF trajectories, and described by

`estimate_MPL(mu,dt_array,single_mut,double_mut)`

The inputs are:

` mu`: mutation probability

` dt_array` : an array containing the generation number when the frequencies are observed

`single_mut`: single mutant frequencies, stored in a T x L matrix where T is the number of observed generations

`double_mut`: double mutant frequencies, stored in a T x L x L matrix 

The outputs are:

`s_MPL`: estimate of the selection coefficients using the MPL algorithm
