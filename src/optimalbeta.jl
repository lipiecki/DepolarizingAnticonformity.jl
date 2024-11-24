function optimalbeta(q::Int, Q::Int)
    optimalβl = ones(3)
    optimalβu = ones(3)
    for (no, type) in enumerate((:dynamic1, :dynamic2, :static2))
        data = load(joinpath("DepolarizingAnticonformityResults", "OutputFiles", "q$(q)_Q$(Q)_$(type).jld2"))
        p = 1.0
        for i in axes(data["polarization_index"], 1)
            for j in axes(data["polarization_index"], 2)
                if data["phase"][i, j] == -2
                    if data["intervention_strength"][i] < p - 1e-9
                        optimalβu[no] = data["probability_outgroup"][j]
                        p = data["intervention_strength"][i]
                    elseif data["intervention_strength"][i] ≈ p
                        optimalβl[no] = data["probability_outgroup"][j]
                    end
                end
            end
        end
    end
    return optimalβl, optimalβu
end
