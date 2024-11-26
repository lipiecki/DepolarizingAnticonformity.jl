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
    abstol = 1e-18 # absolute tolerance for testing stationarity
    ε = 1e-6 # shift in initial conditions to avoid instability in systems that exhibit symmetry breaking
    c0 = [0.5, 0.0, 0.0, 0.5 - ε - Δ]
    if type ∈ dynamictypes # calculations for the dynamic approach
        for i in eachindex(intervention_strength)
            p = intervention_strength[i]
            for j in eachindex(probability_outgroup)
                β = probability_outgroup[j]
                problem = ODEProblem(ode, c0, (0, Inf), [p, β]) 
                sol = solve(problem, Rodas5P(), save_everystep=false, callback=TerminateSteadyState(abstol, 0.0))
                c[i, j, :] .= sol.u[end]
            end
        end
    else # calculations for the static approach
        u0 = zeros(8)
        for i in eachindex(intervention_strength)
            p = intervention_strength[i]
            u0 .= c0[1]*p, c0[2]*p, c0[1]*(1-p), c0[2]*(1-p), c0[3]*p, c0[4]*p, c0[3]*(1-p), c0[4]*(1-p)
            for j in eachindex(probability_outgroup)
                β = probability_outgroup[j]
                problem = ODEProblem(ode, u0, (0, Inf), [p, β])
                sol = solve(problem, Rodas5P(), save_everystep=false, callback=TerminateSteadyState(abstol, 0.0))
                c[i, j, 1] = sol.u[end][1] + sol.u[end][3] 
                c[i, j, 2] = sol.u[end][2] + sol.u[end][4]
                c[i, j, 3] = sol.u[end][5] + sol.u[end][7]
                c[i, j, 4] = sol.u[end][6] + sol.u[end][8]
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
            μ[i, j] = (1 - abs(c1 - c3))*0.5*((c1)/(c1 + c2 + 1e-12) + (c3)/(c3 + c2 + 1e-12)) # `1e-12` allows to stabilize the expression when the denominator approaches 0
            μG = μ[i, j]
            μA = (1 - 2*abs(c1A - c3A))*0.5*((c1A)/(c1A + c2A + 1e-12) + (c3A)/(c3A + c2A + 1e-12))
            μB = (1 - 2*abs(c1B - c3B))*0.5*((c1B)/(c1B + c2B + 1e-12) + (c3B)/(c3B + c2B + 1e-12))
                    
            # classify phases
            if μA > 0.5 && μB > 0.5 && μG > 0.5
                phase[i, j] = 2 # in-group polarization
            elseif μA < 0.5 && μB < 0.5 && μG > 0.5
                phase[i, j] = 1 # between-group polarization
            elseif μA < 0.5 && μB < 0.5 && μG < 0.5
                if (c2 > c1 && c2 > c3)
                    phase[i, j] = -2 # middle-ground consensus
                else
                    phase[i, j] = -1 # pole consensus
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
