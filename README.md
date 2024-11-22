# Depolarizing Power of Anticonformity

`DepolarizingAnticonformity.jl` is a replication package for the preprint A. Lipiecki, K. Weron, Depolarizing Power of Anticonformity (2024).

## Instalation
To install the package, use the following lines in Julia REPL:

```julia
julia> using Pkg
julia> Pkg.add(url = "https://github.com/lipiecki/DepolarizingAnticonformity.jl")
```
There is no need to download the repository manually, you can install it directly to your Julia environment.

## Numerical evolution of the system
Performing the evolution of the dynamical systems corresponding to different model specifications is done by a single function call. `runstudy(q::Int, Q::Int, type::Symbol)` performs the the entire study for specified values of `q` and `Q` (as described in the paper) and of the six examined model specifications, specified by `type`:
- `type=:dynamic1` for Conformity with BC and Anticonformity in the Dynamic (Annealed) approach
- `type=:dynamic2` for Conformity without BC and Anticonformity in the Dynamic (Annealed) approach
- `type=:dynamic3` for Conformity with and without BC in the Dynamic (Annealed) approach
- `type=:static1` for Conformity with BC and Anticonformity in the Static (Quenched) approach
- `type=:static2` for Conformity without BC and Anticonformity in the Static (Quenched) approach
- `type=:static3` for Conformity with and without BC in the Static (Quenched) approach.
The function does not return anything, but saves the output file in `DepolarizingAnticonformityResults/Figures/OutputFiles` (the directory is automatically created if it does not exists). The file is in the `.jld2` format, and stores the data as a dictionary with the following key:value pairs:
- `intervention_strength`: the vector of values describing the strength of intervention (defaults to `0.01:0.0005:0.5`)
- `probability_outgroup`: the vector of values describing the probability of outgroup interaction (defaults to `0.01:0.0005:0.5`)
- `polarization_index`: the matrix of polarization index values, where the first index corresponds to the intervention strength, and the second to probability of outgroup interaction
- `opinion_concentration`: the three-dimensional array of stationary opinion concentrations, where the first index corresponds to the intervention strength, the second to probability of outgroup interaction and the third specifies agent types: 
    - `1` for agents with opinion `-1` in faction `A`
    - `2` for agents with opinion `1` in faction `A`
    - `3` for agents with opinion `-1` in faction `B`
    - `4` for agents with opinion `1` in faction `B`
- `phase`: the matrix identifying the phase, where the first index corresponds to the intervention strength, and the second to probability of outgroup interaction. The phases are denoted using the following notation: 
    - `0` for Middle-Ground Consensus
    - `1` for Pole Consensus
    - `2` for Compromise
    - `3` for Between-Group Polarization
    - `5` for In-Group Polarization
    - (`-1` if the phase has not been classified, which will throw a warning)

Below you can find and example of how to run a single study:

```julia
using DepolarizingAnticonformity
q = 3
Q = 4
type = :dynamic1
DepolarizingAnticonformity.runstudy(q, Q, type)
```

In order to create a polarization map from the results stored in `DepolarizingAnticonformityResults/OutputFiles`, use `polarizationmap(q, Q, type)`, e.g.

```julia
DepolarizingAnticonformity.polarizationmap(3, 4, :dynamic1)
```

To obtain the same results, call:

```julia
DepolarizingAnticonformity.phasemap(3, 4, :dynamic1)
```
The figures are saved in `DepolarizingAnticonformityResults/Figures/Figures` (the directory is automatically created if it does not exists).

## Sensitivity analysis
In addition to the code replicating the results presented in the paper, `DepolarizingAnticonformity.jl` provides the tool for sensitivity analysis, in the form of `sensitivitymap(q, Q, type)` function. It runs the numerical evolution of systems with shifted initial conditions and compares the phase classification against the unperturbed system. The method checks the perturbations of initial conditions in the form of:
$$c_{0B}(0) = \varepsilon + \Delta, \quad c_{1B}(0) = \frac{1}{2} - \varepsilon - \Delta,$$
where $\Delta\in \{10^{-5}, 10^{-4}, 10^{-3}, 10^{-2}, 10^{-1}\}$ and $\varepsilon = 10^{-5}$ is introduced to systems that exhibit symmetry breaking. The method creates a plot which colors the phase sapce by the $-\log_{10}(\Delta)$ required to introduce change: the darker the color, the higher the sensitivity of phase classification to the inital conditions.