function study(q::Int, Q::Int, type::Symbol, Δ::Float64=0.0, prange = 0.01:0.0005:0.5, βrange = 0.01:0.0005:0.5)
    Q > 0 || error("incorrect value of Q")
    (q <= Q && q > floor(Q/2)) || error("incorrect value of q")
    type ∈ (:dynamic1, :dynamic2, :dynamic3, :static1, :static2, :static3)
    de = model(q, Q, Val(type))
    c = zeros(length(prange), length(βrange), 4)
    μ = zeros(length(prange), length(βrange))
    phase = zeros(length(prange), length(βrange))
    ε = (type ∈ (:static1, :static3, :dynamic3)) ? 1e-5 : 0.0
    c0 = [0.5, 0.0, 0.0, 0.5 - ε - Δ]
    if type ∈ (:dynamic1, :dynamic2, :dynamic3)
        @Threads.threads for i in eachindex(prange)
            p = prange[i]
            for j in eachindex(βrange)
                β = βrange[j]
                prob = ODEProblem(de, [0.5, 0.0, 0.0, 0.5 - ε], (0, 2*10^12), [p, β])
                sol = solve(prob, Rosenbrock23(), saveat=[10^12, 2*10^12])
                
                for var_index in eachindex(sol(1))
                    if abs(sol(2)[var_index] - sol(1)[var_index]) > 10^(-12)
                        @warn "System did not converge for p = $(p), and frac = $(β)" 
                    end
                end
                c[i, j, :] .= @view(sol(2)[1:4])
            end
        end
    else
        @Threads.threads for i in eachindex(prange)
            p = prange[i]
            for j in eachindex(βrange)
                β = βrange[j]
                prob = ODEProblem(de, [c0[1]*p, c0[2]*p, c0[1]*(1-p), c0[2]*(1-p), c0[3]*p, c0[4]*p, c0[3]*(1-p), c0[4]*(1-p)], (0, 2*10^12), [p, β])
                sol = solve(prob, Rosenbrock23(), saveat=[10^12, 2*10^12])
                for var_index in eachindex(sol(1))
                    if abs(sol(2)[var_index] - sol(1)[var_index]) > 10^(-12)
                        @warn "System did not converge for p = $(p), and frac = $(β)" 
                    end
                end
                c[i, j, 1] = sol(2)[1] + sol(2)[3] 
                c[i, j, 2] = sol(2)[2] + sol(2)[4]
                c[i, j, 3] = sol(2)[5] + sol(2)[7]
                c[i, j, 4] = sol(2)[6] + sol(2)[8]
            end
        end
    end
    for i in eachindex(prange)
        for j in eachindex(βrange)
            c1A, c3A, c1B, c3B = @view(c[i, j, :])
            c2A = 0.5 - c1A - c3A
            c2B = 0.5 - c1B - c3B
            c1 = c1A + c1B
            c3 = c3A + c3B
            c2 = 1.0 - c1 - c3

            # calculate polarization index
            tol = 1e-9
            μ[i, j] = (1 - abs(c1 - c3))*0.5*((c1)/(c1 + c2 + tol) + (c3)/(c3 + c2 + tol)) # `tol` allows to stabilize the expression when the denominator approaches 0

            μG = μ[i, j]
            μA = (1 - 2*abs(c1A - c3A))*0.5*((c1A)/(c1A + c2A + tol) + (c3A)/(c3A + c2A + tol))
            μB = (1 - 2*abs(c1B - c3B))*0.5*((c1B)/(c1B + c2B + tol) + (c3B)/(c3B + c2B + tol))
                    
            # classify phases
            if μA > 0.5 && μB > 0.5 && μG > 0.5
                phase[i, j] = 4 # in-group polarization
            elseif μA < 0.5 && μB < 0.5 && μG > 0.5
                phase[i, j] = 3 # between-group polarization
            elseif μA > 0.5 && μB > 0.5 && μG < 0.5 && c2 > (c1-tol) && c2 > (c3-tol)
                phase[i, j] = 2 # compromise
            elseif μA < 0.5 && μB < 0.5 && μG < 0.5
                if (c1 > (c2-tol) && c1 > (c3-tol)) || (c3 > (c2-tol) && c3 > (c1-tol))
                    phase[i, j] = 1 # pole consensus
                else
                    phase[i, j] = 0 # middle-ground consensus
            else
                phase[i, j] = -1
                # warn if the phase is unclassified, supress the warning for sensitivity analysis
                (Δ ≈ 0.0) && @warn "unclassified phase at (p=$(prange[i]), β=$(βrange[j])) with cA=$((c1A, c2A, c3A)), cB=$((c1B, c2B, c3B)))"
            end
        end
    end
    return c, μ, phase, prange, βrange
end

function runstudy(q::Int, Q::Int, type::Symbol)
    c, μ, phase, prange, βrange = study(q, Q, type)
    save(joinpath(mkpath(joinpath("DepolarizingAnticonformityResults", "OutputFiles")), "q$(q)_Q$(Q)_$(type).jld2"), "opinion_concentration", c, "polarization_index", μ, "phase", phase, "intervention_strength", prange, "probability_outgroup", βrange)
end
