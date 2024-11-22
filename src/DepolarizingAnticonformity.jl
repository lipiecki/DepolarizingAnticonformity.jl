module DepolarizingAnticonformity

using DifferentialEquations
using JLD2
using Plots

include("model.jl")
include("plotmaps.jl")
include("sensitivity.jl")
include("study.jl")

export
    model
    study
    runstudy
    phasemap
    polarizationmap
    sensitivitymap

end
