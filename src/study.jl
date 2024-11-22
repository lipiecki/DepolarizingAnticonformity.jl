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

            # tolerance for polarization index calculations and phase classification
            tol = 1e-8

            # calculate polarization index
            μ[i, j] = (1 - abs(c1 - c3))*0.5*((c1)/(c1 + c2 + tol) + (c3)/(c3 + c2 + tol)) # `tol` allows to stabilize the expression when the denominator approaches 0

            # classify phases
            if (abs(c1A - c3A) < tol && c2A < (c1A+tol)) && (abs(c1B - c3B) < tol && c2B < (c1B+tol))
                phase[i, j] = 4 # in-group polarization
            elseif (c1A > (c2A+tol) && c2A > (c3A+tol)) && (c3B > (c2B+tol) && c2B > (c1B+tol))
                if abs(c1 - c3) < tol && c2 < (c1+tol)
                    phase[i, j] = 3 # between-group polarization
                else
                    phase[i, j] = 2 # compromise
                end
            elseif (c1A > (c2A+tol) && c2A > (c3A+tol)) && (c1B > (c2B+tol) && c2B > (c3B+tol)) ||
                (c3A > (c2A+tol) && c2A > (c1A+tol)) && (c3B > (c2B+tol) && c2B > (c1B+tol))
                phase[i, j] = 1 # pole consensus
            elseif (c2A > (c1A+tol) && c2A > (c3A+tol)) && (c2B > (c1B+tol) && c2B > (c1B+tol))
                phase[i, j] = 0 # middle-ground consensus
            else
                phase[i, j] = -1
                # warn if the phase is unclassified, supress the warning for sensitivity analysis
                @warn "unclassified phase at (p=$(prange[i]), β=$(βrange[j])) with cA=$((c1A, c2A, c3A)), cB=$((c1B, c2B, c3B)))"
            end
        end
    end
    return c, μ, phase, prange, βrange
end

function runstudy(q::Int, Q::Int, type::Symbol)
    println("v2")
    c, μ, phase, prange, βrange = study(q, Q, type)
    save(joinpath(mkpath(joinpath("DepolarizingAnticonformityResults", "OutputFiles")), "q$(q)_Q$(Q)_$(type).jld2"), "opinion_concentration", c, "polarization_index", μ, "phase", phase, "intervention_strength", prange, "probability_outgroup", βrange)
end
