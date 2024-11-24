function study(q::Int, Q::Int, type::Symbol; Δ::Float64=0.0, intervention_strength::StepRangeLen{Float64}=0.01:0.0005:0.5, probability_outgroup::StepRangeLen{Float64}=0.01:0.0005:0.5)
    if !haskey(ODEs, (q, Q, type))
        if Q > Qmax
            error("Q (=$(Q)) exceeds the predefined maximum value ($(Qmax))")
        elseif q > Q
            error("q (=$(q)) exceeds Q (=$(Q))")
        elseif q < Int(floor(Q/2)+1)
            error("for Q=$(Q) q cannot be smaller than $(Int(floor(Q/2)+1))")
        elseif type ∉ (statictypes..., dynamictypes...)
            error("unknown model type")
        else
            error("undefined error")
        end
    end
    ode = ODEs[(q, Q, type)] # load the appropriate ODE system from the dictionary
    c = zeros(length(intervention_strength), length(probability_outgroup), 4) # array for storing stationary opinion concentrations
    μ = zeros(length(intervention_strength), length(probability_outgroup)) # array for storing stationary polarization index
    phase = zeros(length(intervention_strength), length(probability_outgroup)) # array for storing phase signature
    T = 10^12 # time for which the system is evolved
    tol = 1e-8 # tolerance for testing stationarity (also used in polarization index calculations and phase classification)
    ε = 1e-6 # shift in initial conditions allowing to avoid instability in the case of symmetry breaking
    c0 = [0.5, 0.0, 0.0, 0.5 - ε - Δ]
    if type ∈ dynamictypes # calculations for the dynamic approach
        @Threads.threads for i in eachindex(intervention_strength)
            p = intervention_strength[i]
            for j in eachindex(probability_outgroup)
                β = probability_outgroup[j]
                problem = ODEProblem(ode, c0, (0, T), [p, β]) 
                sol = solve(problem, Rosenbrock23(), save_everystep=false)
                c[i, j, :] .= sol.u[2]
                du = ones(4)
                ode(du, sol.u[2], [p, β], Inf)
                all(abs.(du) .< tol) || @warn "nonstationary solution detected at p=$(p) and β=$(β)" 
            end
        end
    else # calculations for the static approach
        @Threads.threads for i in eachindex(intervention_strength)
            p = intervention_strength[i]
            for j in eachindex(probability_outgroup)
                β = probability_outgroup[j]
                problem = ODEProblem(ode, [c0[1]*p, c0[2]*p, c0[1]*(1-p), c0[2]*(1-p), c0[3]*p, c0[4]*p, c0[3]*(1-p), c0[4]*(1-p)], (0, T), [p, β])
                sol = solve(problem, Rosenbrock23(), save_everystep=false)
                c[i, j, 1] = sol.u[2][1] + sol.u[2][3] 
                c[i, j, 2] = sol.u[2][2] + sol.u[2][4]
                c[i, j, 3] = sol.u[2][5] + sol.u[2][7]
                c[i, j, 4] = sol.u[2][6] + sol.u[2][8]
                du = ones(8)
                ode(du, sol.u[2], [p, β], Inf)
                all(abs.(du) .< tol) || @warn "nonstationary solution detected at p=$(p) and β=$(β)" 
            end
        end
    end
    for i in eachindex(intervention_strength)
        for j in eachindex(probability_outgroup)
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
                phase[i, j] = 2 # in-group polarization
            elseif μA < 0.5 && μB < 0.5 && μG > 0.5
                phase[i, j] = 1 # between-group polarization
            elseif μA < 0.5 && μB < 0.5 && μG < 0.5
                if (c1 > (c2-tol) && c1 > (c3-tol)) || (c3 > (c2-tol) && c3 > (c1-tol))
                    phase[i, j] = -1 # pole consensus
                else
                    phase[i, j] = -2 # middle-ground consensus
                end
            end
        end
    end
    println("\nq=$(q), Q=$(Q), type=$(type), Δ=$(Δ)")
    println("-"^30)
    println("phase\t| %")
    println("-"^30)
    println("IGP\t| ", round(100*count(phase .== 2)/length(phase), digits=5))
    println("BGP\t| ", round(100*count(phase .== 1)/length(phase), digits=5))
    println("PC\t| ", round(100*count(phase .== -1)/length(phase), digits=5))
    println("MGC\t| ", round(100*count(phase .== -2)/length(phase), digits=5))
    println("-"^30)
    println("unclassified states: ", round(100*count(phase .== 0)/length(phase), digits=5), "%")
    return c, μ, phase, intervention_strength, probability_outgroup
end

function runstudy(q::Int, Q::Int, type::Symbol)
    c, μ, phase, intervention_strength, probability_outgroup = study(q, Q, type)
    save(joinpath(mkpath(joinpath("DepolarizingAnticonformityResults", "OutputFiles")), "q$(q)_Q$(Q)_$(type).jld2"), "opinion_concentration", c, "polarization_index", μ, "phase", phase, "intervention_strength", intervention_strength, "probability_outgroup", probability_outgroup)
end
