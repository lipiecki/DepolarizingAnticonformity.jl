function optimalbeta(q::Int, Q::Int)
    optimalβ = ones(3)
    for (no, type) in enumerate((:dynamic1, :dynamic2, :static2))
        data = load(joinpath("DepolarizingAnticonformityResults", "OutputFiles", "q$(q)_Q$(Q)_$(type).jld2"))
        p = 1.0
        for i in axes(data["polarization_index"], 1)
            for j in axes(data["polarization_index"], 2)
                if data["phase"][i, j] == -2 && data["intervention_strength"][i] < p + 1e-9
                    optimalβ[no] = data["probability_outgroup"][j]
                    p = data["intervention_strength"][i]
                end
            end
        end
    end
    return optimalβ
end
