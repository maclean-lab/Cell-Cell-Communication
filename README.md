# Cell-Cell-Communication

### Summary

We provide code associated with the development of hybrid multiscale models of cell-cell communication. The models are developed in [Julia](https://julialang.org/) for a specific gene regulatory network in hematopoiesis that controls myeloid cell fates; however the framework is general and can be applied in various contexts to gene regulatory networks linked to single cell phenotypes. The hybrid model consists of coupled continuous and discrete dynamical processes:

- Continuous gene regulatory  dynamics inside cells (GATA1-PU.1 network), modeled with nonlinear differential equations as described previously (Chickarmane et al., 2009). 
- Discrete cell-cell communication between cells, where signals sent between cells are modeled by a Poisson process.

In general, the form of the external signals can be defined by the user. In the paper, we consider two signals: consensus signals ("be like me") that recruit neighboring cells to commit to the same lineage, or dissensus signals ("be unlike me") that push neighbors to commit to the alternative lineage. Cell-cell signaling network topologies are defined by a matrix representation
    - 1 corresponds to a consensus signal
    - -1 corresponds to a dissensus signal
    - 0 corresponds to no signaling between cells. 

Simulation of individual model runs (trajectories) is performed in `Topology_Trajectories.jl`. Running sets of simulations over a range of parameters is performed in`Topology_Probabilities.jl` if no noise is modeled, in  or in `Simple_Topologies_w_Noise_Probabilities.jl` with noise (either intrinsic or extrinsic) added. To run large numbers of simulations over many parameter values, we recommend using hpc. The notebook `Plotting_Trajectories_and_Probabilities.ipynb` gives the code used to produce the plots in the paper, both for single trajectory simulations and approximate probability distributions for large sets of parameters.


### Paper

**A single-cell resolved cell-cell communication model explains lineage commitment in hematopoiesis**

Megan K Franke and Adam L MacLean

bioRxiv: https://www.biorxiv.org/content/10.1101/2021.03.31.437948v1

doi: https://doi.org/10.1101/2021.03.31.437948


### Dependencies

- [`DifferentialEquations`](https://github.com/SciML/DifferentialEquations.jl) - for differential equation modeling 
- [`Distributions`](https://github.com/JuliaStats/Distributions.jl) - for probability distributions
- [`DelimitedFiles`]() - for I/O 
- [`Plots`](https://github.com/JuliaPlots/Plots.jl), [`LsqFit`](https://github.com/JuliaNLSolvers/LsqFit.jl) for plotting and curve fitting


### References

- V. Chickarmane, T. Enver, C. Peterson, Computational modeling of the hematopoietic erythroid- myeloid switch reveals insights into cooperativity, priming, and irreversibility, *PLoS Comput Biol* 5 (1) (2009) e1000268.
- J. Bezanson, A. Edelman, S. Karpinski, V. B. Shah, Julia: A fresh approach to numerical computing, *SIAM Review* 59 (1) (2017) 65–98.
- C. Rackauckas, Q. Nie, Differentialequations.jl–a performant and feature-rich ecosystem for solving differential equations in julia, *Journal of Open Research Software* 5 (1) (2017).
