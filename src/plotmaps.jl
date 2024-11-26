function polarizationmap(q::Int, Q::Int, type::Symbol)
    data = load(joinpath("DepolarizingAnticonformityResults", "OutputFiles", "q$(q)_Q$(Q)_$(type).jld2"))
    plt = heatmap(data["intervention_strength"], data["probability_outgroup"], data["polarization_index"]', c=cgrad(:RdBu, rev=true), clims=(0, 1))
    plot!(plt, framestyle=:grid, colorbar=false, size=(300, 300), xlabel="p", ylabel="β")
    savefig(plt, joinpath(mkpath(joinpath("DepolarizingAnticonformityResults", "Figures")), "polarization_q$(q)_Q$(Q)_$(type).pdf"))
end

function phasemap(q::Int, Q::Int, type::Symbol)
    data = load(joinpath("DepolarizingAnticonformityResults", "OutputFiles", "q$(q)_Q$(Q)_$(type).jld2"))
    plt = heatmap(data["intervention_strength"], data["probability_outgroup"], data["phase"]', c=cgrad(:RdBu, 5, categorical=true, rev=true), clims=(-2, 2))
    plot!(plt, framestyle=:grid, colorbar=false, size=(300, 300), xlabel="p", ylabel="β")
    mkpath(joinpath("DepolarizingAnticonformityResults", "Figures"))
    savefig(plt, joinpath(mkpath(joinpath("DepolarizingAnticonformityResults", "Figures")), "phase_q$(q)_Q$(Q)_$(type).pdf"))
end
