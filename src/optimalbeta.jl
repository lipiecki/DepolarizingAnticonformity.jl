function optimalbeta(q::Int, Q::Int)
    for (no, type) in enumerate((:dynamic2, :static2))
        data = load(joinpath("DepolarizingAnticonformityResults", "OutputFiles", "q$(q)_Q$(Q)_$(type).jld2"))
        βc = 1.01
        p = 1.01
        for i in axes(data["polarization_index"], 1)
            for j in axes(data["polarization_index"], 2)
                if data["phase"][i, j] == -2
                    if !(data["intervention_strength"][i] ≈ p) && data["intervention_strength"][i] < p 
                        βc = data["probability_outgroup"][j]
                        p = data["intervention_strength"][i]
                    end
                end
            end
        end
        println("Critical β for $(type): $(βc)")
    end
    type = :dynamic1
    data = load(joinpath("DepolarizingAnticonformityResults", "OutputFiles", "q$(q)_Q$(Q)_$(type).jld2"))
    βl, βu = -0.01, 1.01
    p = 1.01
    for i in axes(data["polarization_index"], 1)
        for j in axes(data["polarization_index"], 2)
            if data["phase"][i, j] == -2
                if (data["intervention_strength"][i] < p)
                    βl = min(βl, data["probability_outgroup"][j])
                    βu = max(βu, data["probability_outgroup"][j])
                end
            end
        end
    end
    println("Optimal β for $(type): [$(βl), $(βu)]")
end
