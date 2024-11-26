module DepolarizingAnticonformity

using DifferentialEquations
using DiffEqCallbacks
using JLD2
using Plots

const statictypes = (:static1, :static2, :static3)
const dynamictypes = (:dynamic1, :dynamic2, :dynamic3)

# generate the ODE systems for every q, Q and model type, using definitions from `model.jl`
include("model.jl")
const Qmax = 10 # maximum Q for which the ODE system is generated
const ODEs = Dict((q, Q, type) => model(q, Q, Val(type)) for Q in 1:Qmax for q in Int(floor(Q/2)+1):Q, 
                type in (statictypes..., dynamictypes...))

include("plotmaps.jl")
include("sensitivity.jl")
include("study.jl")

export
    study,
    runstudy,
    phasemap,
    polarizationmap,
    sensitivitymap
end
