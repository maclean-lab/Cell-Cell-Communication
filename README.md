# Cell-Cell-Communication

### Summary

Multiscale model of single cell dynamics of the granulocyte-monocyte vs. megakaryocyte-erythrocyte fate decision with cell-cell communication. This hybrid model is implemented in the [Julia](https://julialang.org/) programming language and consists of continuous and discrete processes:

- Continuous intracellular GATA1-PU.1 gene regulatory network dynamics, modeled using the bistable system of ODEs defined in Chickarmane et al. 2009. 
- Discrete signals sent between cells following a Poisson process

Cells can send consensus signals, recruiting neighbors to converge to the same lineage, or dissensus signals, pushing neighbors to converge to the opposite lineage. Cell-cell signaling topologies are encoded as a matrix where 1 corresponds to a consensus signal, -1 corresponds to a dissensus signal, and 0 corresponds to no communication. 

Single trajectory simulations can be performed using Topology_Trajectories.jl. Running large numbers of simulations over a range of parameters can be done without noise using Topology_Probabilities.jl or with noise (either intrinsic or extrinsic) using Simple_Topologies_w_Noise_Probabilities.jl. Large numbers of iterations over many parameter values should not be run locally. The notebook Plotting_Trajectories_and_Probabilities.ipynb plots both single trajectory simulations and approximate probabilities. 

### Paper

**A single-cell resolved cell-cell communication model explains lineage commitment in hematopoiesis**

Megan K Franke and Adam L MacLean

[https://www.biorxiv.org/content/10.1101/2021.03.31.437948v1]

### Requirements

1. Julia 1.6.0
2. DifferentialEquations, Distributions, Plots, DelimitedFiles, LsqFit

### References

- V. Chickarmane, T. Enver, C. Peterson, Computational modeling of the hematopoietic erythroid- myeloid switch reveals insights into cooperativity, priming, and irreversibility, *PLoS Comput Biol* 5 (1) (2009) e1000268.
