module DepolarizingAnticonformity

using DifferentialEquations
using JLD2
using Plots

include("model.jl")
include("optimalbeta.jl")
include("plotmaps.jl")
include("sensitivity.jl")
include("study.jl")

export
    model,
    study,
    runstudy,
    optimalbeta,
    phasemap,
    polarizationmap,
    sensitivitymap,
end
