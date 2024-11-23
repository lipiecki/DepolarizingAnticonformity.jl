function study(q::Int, Q::Int, type::Symbol, Δ::Float64=0.0, prange = 0.01:0.0005:0.5, βrange = 0.01:0.0005:0.5)
    Q > 0 || error("incorrect value of Q")
    (q <= Q && q > floor(Q/2)) || error("incorrect value of q")
    type ∈ (:dynamic1, :dynamic2, :dynamic3, :static1, :static2, :static3)
    de = model(q, Q, Val(type))
    c = zeros(length(prange), length(βrange), 4)
    μ = zeros(length(prange), length(βrange))
    phase = zeros(length(prange), length(βrange))
    ε = 1e-6
    c0 = [0.5, 0.0, 0.0, 0.5 - ε - Δ]
    tol = 1e-12
    T = 1e12
    if type ∈ (:dynamic1, :dynamic2, :dynamic3)
        @Threads.threads for i in eachindex(prange)
            p = prange[i]
            for j in eachindex(βrange)
                β = βrange[j]
                prob = ODEProblem(de, [0.5, 0.0, 0.0, 0.5 - ε], (0, 2T), [p, β])
                sol = solve(prob, Rosenbrock23(), saveat=[T, 2T])
                
                for var_index in eachindex(sol(1))
                    if abs(sol(2)[var_index] - sol(1)[var_index]) > tol
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
                prob = ODEProblem(de, [c0[1]*p, c0[2]*p, c0[1]*(1-p), c0[2]*(1-p), c0[3]*p, c0[4]*p, c0[3]*(1-p), c0[4]*(1-p)], (0, 2T), [p, β])
                sol = solve(prob, Rosenbrock23(), saveat=[T, 2T])
                for var_index in eachindex(sol(1))
                    if abs(sol(2)[var_index] - sol(1)[var_index]) > tol
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
            μ[i, j] = (1 - abs(c1 - c3))*0.5*((c1)/(c1 + c2 + tol) + (c3)/(c3 + c2 + tol)) # `tol` allows to stabilize the expression when the denominator approaches 0
            μG = μ[i, j]
            μA = (1 - 2*abs(c1A - c3A))*0.5*((c1A)/(c1A + c2A + tol) + (c3A)/(c3A + c2A + tol))
            μB = (1 - 2*abs(c1B - c3B))*0.5*((c1B)/(c1B + c2B + tol) + (c3B)/(c3B + c2B + tol))
                    
            # classify phases
            if μA > 0.5 && μB > 0.5 && μG > 0.5
                phase[i, j] = 4 # in-group polarization
            elseif μA < 0.5 && μB < 0.5 && μG > 0.5
                phase[i, j] = 3 # between-group polarization
            elseif μA < 0.5 && μB < 0.5 && μG < 0.5
                if (c1 > (c2-tol) && c1 > (c3-tol)) || (c3 > (c2-tol) && c3 > (c1-tol))
                    phase[i, j] = 2 # pole consensus
                else
                    phase[i, j] = 1 # middle-ground consensus
                end
            end
        end
    end
    println("q=$(q), Q=$(Q), type=$(type), Δ=$(Δ)")
    println("-"^30)
    println("phase\t| %")
    println("-"^30)
    println("IGP\t| ", round(100*count(phase .== 5)/length(phase), digits=9))
    println("BGP\t| ", round(100*count(phase .== 4)/length(phase), digits=9))
    println("Comp\t| ", round(100*count(phase .== 3)/length(phase), digits=9))
    println("PC\t| ", round(100*count(phase .== 2)/length(phase), digits=9))
    println("MGC\t| ", round(100*count(phase .== 1)/length(phase), digits=9))
    println("unknown\t| ", round(100*count(phase .== 0)/length(phase), digits=9))
    println("-"^30)
    println()
    return c, μ, phase, prange, βrange
end

function runstudy(q::Int, Q::Int, type::Symbol)
    c, μ, phase, prange, βrange = study(q, Q, type)
    save(joinpath(mkpath(joinpath("DepolarizingAnticonformityResults", "OutputFiles")), "q$(q)_Q$(Q)_$(type).jld2"), "opinion_concentration", c, "polarization_index", μ, "phase", phase, "intervention_strength", prange, "probability_outgroup", βrange)
end
