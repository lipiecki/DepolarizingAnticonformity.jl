# Depolarizing Power of Anticonformity
`DepolarizingAnticonformity.jl` is a replication package for the paper [*Depolarizing Power of Anticonformity*, A. Lipiecki & K. Weron (2025)](https://doi.org/10.1016/j.eswa.2025.127879).

Preprint version is openly available in [WORMS](https://worms.pwr.edu.pl/RePEc/ahh/wpaper/WORMS_25_03.pdf)

If you find this repository useful in your research, please consider adding the following citation:
```bibtex
@article{lipiecki:sznajdweron:2025,
    title = {Depolarizing power of anticonformity},
    journal = {Expert Systems with Applications},
    volume = {285},
    pages = {127879},
    year = {2025},
    issn = {0957-4174},
    doi = {10.1016/j.eswa.2025.127879},
    author = {Arkadiusz Lipiecki and Katarzyna Sznajd-Weron}
}
```

## Installation
To install the package, use the following lines in Julia REPL:
```julia
julia> using Pkg
julia> Pkg.add(url = "https://github.com/lipiecki/DepolarizingAnticonformity.jl")
```
there is no need to download the repository manually, you can install it directly to your Julia environment.

## Numerical evolution of the system
Performing the evolution of the dynamical systems corresponding to different model specifications is done by a single function call. `study(q::Int, Q::Int, type::Symbol)` performs the the entire study for specified `q` and `Q` (as described in the paper) and model specification `type`, which can take one of the following values:
- `type=:dynamic1` for Conformity with BC and Anticonformity in the Dynamic (Annealed) approach
- `type=:dynamic2` for Conformity without BC and Anticonformity in the Dynamic (Annealed) approach
- `type=:dynamic3` for Conformity with and without BC in the Dynamic (Annealed) approach
- `type=:static1` for Conformity with BC and Anticonformity in the Static (Quenched) approach
- `type=:static2` for Conformity without BC and Anticonformity in the Static (Quenched) approach
- `type=:static3` for Conformity with and without BC in the Static (Quenched) approach.

The function prints the summary of phase classification and saves the output file in `DepolarizingAnticonformityResults/OutputFiles` (the directory is automatically created if it does not exists). The file is in the .jld2 format, and stores the data as a dictionary with the following key:value pairs:
- `intervention_strength`: the vector of values describing the strength of intervention (`0.01:0.0005:0.5` by default)
- `probability_outgroup`: the vector of values describing the probability of outgroup interaction (`0.01:0.0005:0.5` by default)
- `polarization_index`: the matrix of polarization index values, where the first index corresponds to the intervention strength, and the second to probability of outgroup interaction
- `opinion_concentration`: the three-dimensional array of stationary opinion concentrations, where the first index corresponds to the intervention strength, the second to probability of outgroup interaction and the third specifies agent types: 
    - `1` for agents in faction `A` holding opinion `-1`
    - `2` for agents in faction `A` holding opinion `1`
    - `3` for agents in faction `B` holding opinion `-1`
    - `4` for agents in faction `B` holding opinion `1`
- `phase`: the matrix identifying the phase, where the first index corresponds to the intervention strength, and the second to probability of outgroup interaction. The phases are denoted using the following notation: 
    - `2` for In-Group Polarization
    - `1` for Between-Group Polarization
    - `-1` for Pole Consensus
    - `-2` for Middle-Ground Consensus
    - `0` if the phase has not been classified

Below you can find an example of how to run a single study:
```julia
using DepolarizingAnticonformity # loads the package
study(3, 4, :dynamic1)
#=
q=3, Q=4, type=dynamic1, Δ=0.0
------------------------------
phase	| %
------------------------------
IGP	| 42.32757
BGP	| 4.20715
PC	| 0.0
MGC	| 53.44564
------------------------------
unclassified states: 0.01964%
=#
```
## Plots
In order to generate a polarization map from the results stored in `DepolarizingAnticonformityResults/OutputFiles`, call `polarizationmap(q, Q, type)`, for example:
```julia
polarizationmap(3, 4, :dynamic1)
```
and to generate a phase map for the same results, call `phasemap(q, Q, type)`:
```julia
phasemap(3, 4, :dynamic1)
```
the figures are saved in `DepolarizingAnticonformityResults/Figures` (the directory is automatically created if it does not exists).

## Sensitivity analysis
In addition to the code replicating the results presented in the paper, `DepolarizingAnticonformity.jl` provides the tool for sensitivity analysis. It runs the numerical evolution of systems with shifted initial conditions and compares the phase classification against the unperturbed system. The method checks the perturbations of initial conditions in the form of:
$$c_{0B}(0) = \varepsilon + \Delta, \quad c_{1B}(0) = \frac{1}{2} - \varepsilon - \Delta,$$
where $\Delta\in$ { $10^{-5}, 10^{-4}, 10^{-3}, 10^{-2}, 10^{-1}$ } and $\varepsilon = 10^{-6}$ is introduced to avoid instabilities in systems that exhibit symmetry breaking. The method creates a plot which colors the phase space by the $-\log_{10}(\Delta)$ required to introduce change: the darker the color, the higher the sensitivity of phase classification to the inital conditions.

To run the sensitivity analysis, call `sensitivitystudy(q, Q, type)`, for example:
```julia
sensitivitystudy(3, 4, :dynamic1)
```
and plot the corresponding sensitivity map with `sensitivitymap(q, Q, type)`:
```julia
sensitivitymap(3, 4, :dynamic1)
```

## Notes
For convenience, the notation for opinion values used in the source code differes from the one in the paper. In the implementation, different values of the three-state ranked opinion are denoted by 1, 2 and 3.
