function sensitivitymap(q::Int, Q::Int, type::Symbol)
    _, _, unperturbed_phase, prange, βrange = study(q, Q, type)
    sensitivity = zeros(Int, size(unperturbed_phase))
    for logΔ in 1:1:5
        _, _, perturbed_phase, prange, βrange = study(q, Q, type, 10.0^(-logΔ))
        for i in axes(unperturbed_phase, 1)
            for j in axes(unperturbed_phase, 2)
                if unperturbed_phase[i, j] != perturbed_phase[i, j]
                    sensitivity[i, j] = logΔ
                end
            end
        end
    end
    save(joinpath(mkpath(joinpath("DepolarizingAnticonformityResults", "OutputFiles")), "q$(q)_Q$(Q)_$(type).jld2"), "sensitivity", sensitivity, "prange", prange, "βrange", βrange)
    plt = heatmap(prange, βrange, sensitivity', c = cgrad(:Blues, 6, categorical=true), clims = (0, 5))
    plot!(plt, framestyle=:grid, colorbar=true, size=(330, 300))
    savefig(plt, joinpath(mkpath(joinpath("DepolarizingAnticonformityResults", "Figures")), "sensitivity_q$(q)_Q$(Q)_$(type).pdf"))
end
