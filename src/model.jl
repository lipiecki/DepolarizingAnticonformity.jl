function model(q::Int, Q::Int, ::Val{:dynamic1})
    SA1 = "(0"
    SA2 = "(0"
    SA3 = "(0"
    SB1 = "(0"
    SB2 = "(0"
    SB3 = "(0"
    for i in q:Q
        SA1 = SA1*" + "*"$(binomial(Q, i))*((1-β)*2*c1A + β*2*c1B)^$(i)*(1 - (1-β)*2*c1A - β*2*c1B)^$(Q-i)"
        SA2 = SA2*" + "*"$(binomial(Q, i))*((1-β)*2*c2A + β*2*c2B)^$(i)*(1 - (1-β)*2*c2A - β*2*c2B)^$(Q-i)"
        SA3 = SA3*" + "*"$(binomial(Q, i))*((1-β)*2*c3A + β*2*c3B)^$(i)*(1 - (1-β)*2*c3A - β*2*c3B)^$(Q-i)"
        SB1 = SB1*" + "*"$(binomial(Q, i))*((1-β)*2*c1B + β*2*c1A)^$(i)*(1 - (1-β)*2*c1B - β*2*c1A)^$(Q-i)"
        SB2 = SB2*" + "*"$(binomial(Q, i))*((1-β)*2*c2B + β*2*c2A)^$(i)*(1 - (1-β)*2*c2B - β*2*c2A)^$(Q-i)"
        SB3 = SB3*" + "*"$(binomial(Q, i))*((1-β)*2*c3B + β*2*c3A)^$(i)*(1 - (1-β)*2*c3B - β*2*c3A)^$(Q-i)"
    end
    SA1 = SA1*")"
    SA2 = SA2*")"
    SA3 = SA3*")"
    SB1 = SB1*")"
    SB2 = SB2*")"
    SB3 = SB3*")"

    def = """
    \tfunction (du, u, params, t)
    \t\tc1A = u[1]
    \t\tc3A = u[2]
    \t\tc1B = u[3]
    \t\tc3B = u[4]
    \t\tc2A = 0.5 - c1A - c3A
    \t\tc2B = 0.5 - c1B - c3B

    \t\tβ = params[2]
    \t\tp = params[1]

    \t\tdu[1] = p*(0.5*c2A*$(SA2) - c1A*$(SA1)) + (1-p)*(c2A*$(SA1) - c1A*$(SA2))
    \t\tdu[2] = p*(0.5*c2A*$(SA2) - c3A*$(SA3)) + (1-p)*(c2A*$(SA3) - c3A*$(SA2))
    \t\tdu[3] = p*(0.5*c2B*$(SB2) - c1B*$(SB1)) + (1-p)*(c2B*$(SB1) - c1B*$(SB2))
    \t\tdu[4] = p*(0.5*c2B*$(SB2) - c3B*$(SB3)) + (1-p)*(c2B*$(SB3) - c3B*$(SB2))

    \t\tnothing
    \tend"""
    return eval(Meta.parse(def))
end

function model(q::Int, Q::Int, ::Val{:dynamic2})
    SA1 = "(0"
    SA2 = "(0"
    SA3 = "(0"
    SB1 = "(0"
    SB2 = "(0"
    SB3 = "(0"
    SAleft = "(0"
    SAright = "(0"
    SBleft = "(0"
    SBright = "(0"
    for i in q:Q
        SA1 = SA1*" + "*"$(binomial(Q, i))*((1-β)*2*c1A + β*2*c1B)^$(i)*(1 - (1-β)*2*c1A - β*2*c1B)^$(Q-i)"
        SA2 = SA2*" + "*"$(binomial(Q, i))*((1-β)*2*c2A + β*2*c2B)^$(i)*(1 - (1-β)*2*c2A - β*2*c2B)^$(Q-i)"
        SA3 = SA3*" + "*"$(binomial(Q, i))*((1-β)*2*c3A + β*2*c3B)^$(i)*(1 - (1-β)*2*c3A - β*2*c3B)^$(Q-i)"
        SB1 = SB1*" + "*"$(binomial(Q, i))*((1-β)*2*c1B + β*2*c1A)^$(i)*(1.0 - (1-β)*2*c1B - β*2*c1A)^$(Q-i)"
        SB2 = SB2*" + "*"$(binomial(Q, i))*((1-β)*2*c2B + β*2*c2A)^$(i)*(1.0 - (1-β)*2*c2B - β*2*c2A)^$(Q-i)"
        SB3 = SB3*" + "*"$(binomial(Q, i))*((1-β)*2*c3B + β*2*c3A)^$(i)*(1.0 - (1-β)*2*c3B - β*2*c3A)^$(Q-i)"

        SAleft = SAleft*" + "*"$(binomial(Q, i))*((1-β)*2*(0.5 - c3A) + β*2*(0.5 - c3B))^$(i)*(1 - (1-β)*2*(0.5 - c3A) - β*2*(0.5 - c3B))^$(Q-i)"
        SAright = SAright*" + "*"$(binomial(Q, i))*((1-β)*2*(0.5 - c1A) + β*2*(0.5 - c1B))^$(i)*(1 - (1-β)*2*(0.5 - c1A) - β*2*(0.5 - c1B))^$(Q-i)"
        SBleft = SBleft*" + "*"$(binomial(Q, i))*((1-β)*2*(0.5 - c3B) + β*2*(0.5 - c3A))^$(i)*(1 - (1-β)*2*(0.5 - c3B) - β*2*(0.5 - c3A))^$(Q-i)"
        SBright = SBright*" + "*"$(binomial(Q, i))*((1-β)*2*(0.5 - c1B) + β*2*(0.5 - c1A))^$(i)*(1 - (1-β)*2*(0.5 - c1B) - β*2*(0.5 - c1A))^$(Q-i)"
    end
    SA1 = SA1*")"
    SA2 = SA2*")"
    SA3 = SA3*")"
    SB1 = SB1*")"
    SB2 = SB2*")"
    SB3 = SB3*")"
    SAleft = SAleft*")"
    SAright = SAright*")"
    SBleft = SBleft*")"
    SBright = SBright*")"

    def = """
    \tfunction (du, u, params, t)
    \t\tc1A = u[1]
    \t\tc3A = u[2]
    \t\tc1B = u[3]
    \t\tc3B = u[4]
    \t\tc2A = 0.5 - c1A - c3A
    \t\tc2B = 0.5 - c1B - c3B

    \t\tβ = params[2]
    \t\tp = params[1]

    \t\tdu[1] = p*(0.5*c2A*$(SA2) - c1A*$(SA1)) + (1-p)*(c2A*$(SA1) - c1A*$(SAright))
    \t\tdu[2] = p*(0.5*c2A*$(SA2) - c3A*$(SA3)) + (1-p)*(c2A*$(SA3) - c3A*$(SAleft))
    \t\tdu[3] = p*(0.5*c2B*$(SB2) - c1B*$(SB1)) + (1-p)*(c2B*$(SB1) - c1B*$(SBright))
    \t\tdu[4] = p*(0.5*c2B*$(SB2) - c3B*$(SB3)) + (1-p)*(c2B*$(SB3) - c3B*$(SBleft))

    \t\tnothing
    \tend"""
    return eval(Meta.parse(def))
end

function model(q::Int, Q::Int, ::Val{:dynamic3})
    SA1 = "(0"
    SA2 = "(0"
    SA3 = "(0"
    SB1 = "(0"
    SB2 = "(0"
    SB3 = "(0"
    SAleft = "(0"
    SAright = "(0"
    SBleft = "(0"
    SBright = "(0"
    for i in q:Q
        SA1 = SA1*" + "*"$(binomial(Q, i))*((1-β)*2*c1A + β*2*c1B)^$(i)*(1 - (1-β)*2*c1A - β*2*c1B)^$(Q-i)"
        SA2 = SA2*" + "*"$(binomial(Q, i))*((1-β)*2*c2A + β*2*c2B)^$(i)*(1 - (1-β)*2*c2A - β*2*c2B)^$(Q-i)"
        SA3 = SA3*" + "*"$(binomial(Q, i))*((1-β)*2*c3A + β*2*c3B)^$(i)*(1 - (1-β)*2*c3A - β*2*c3B)^$(Q-i)"
        SB1 = SB1*" + "*"$(binomial(Q, i))*((1-β)*2*c1B + β*2*c1A)^$(i)*(1 - (1-β)*2*c1B - β*2*c1A)^$(Q-i)"
        SB2 = SB2*" + "*"$(binomial(Q, i))*((1-β)*2*c2B + β*2*c2A)^$(i)*(1 - (1-β)*2*c2B - β*2*c2A)^$(Q-i)"
        SB3 = SB3*" + "*"$(binomial(Q, i))*((1-β)*2*c3B + β*2*c3A)^$(i)*(1 - (1-β)*2*c3B - β*2*c3A)^$(Q-i)"

        SAleft = SAleft*" + "*"$(binomial(Q, i))*((1-β)*2*(0.5 - c3A) + β*2*(0.5 - c3B))^$(i)*(1 - (1-β)*2*(0.5 - c3A) - β*2*(0.5 - c3B))^$(Q-i)"
        SAright = SAright*" + "*"$(binomial(Q, i))*((1-β)*2*(0.5 - c1A) + β*2*(0.5 - c1B))^$(i)*(1 - (1-β)*2*(0.5 - c1A) - β*2*(0.5 - c1B))^$(Q-i)"
        SBleft = SBleft*" + "*"$(binomial(Q, i))*((1-β)*2*(0.5 - c3B) + β*2*(0.5 - c3A))^$(i)*(1 - (1-β)*2*(0.5 - c3B) - β*2*(0.5 - c3A))^$(Q-i)"
        SBright = SBright*" + "*"$(binomial(Q, i))*((1-β)*2*(0.5 - c1B) + β*2*(0.5 - c1A))^$(i)*(1 - (1-β)*2*(0.5 - c1B) - β*2*(0.5 - c1A))^$(Q-i)"
    end
    SA1 = SA1*")"
    SA2 = SA2*")"
    SA3 = SA3*")"
    SB1 = SB1*")"
    SB2 = SB2*")"
    SB3 = SB3*")"
    SAleft = SAleft*")"
    SAright = SAright*")"
    SBleft = SBleft*")"
    SBright = SBright*")"

    def = """
    \tfunction (du, u, params, t)
    \t\tc1A = u[1]
    \t\tc3A = u[2]
    \t\tc1B = u[3]
    \t\tc3B = u[4]
    \t\tc2A = 0.5 - c1A - c3A
    \t\tc2B = 0.5 - c1B - c3B
    
    \t\tβ = params[2]
    \t\tp = params[1]

    \t\tdu[1] = p*(0.5*c2A*$(SA1) - c1A*$(SAright)) + (1-p)*(c2A*$(SA1) - c1A*$(SA2))
    \t\tdu[2] = p*(0.5*c2A*$(SA3) - c3A*$(SAleft)) + (1-p)*(c2A*$(SA3) - c3A*$(SA2))
    \t\tdu[3] = p*(0.5*c2B*$(SB1) - c1B*$(SBright)) + (1-p)*(c2B*$(SB1) - c1B*$(SB2))
    \t\tdu[4] = p*(0.5*c2B*$(SB3) - c3B*$(SBleft)) + (1-p)*(c2B*$(SB3) - c3B*$(SB2))

    \t\tnothing
    \tend"""
    return eval(Meta.parse(def))
end

function model(q::Int, Q::Int, ::Val{:static1})
    SA1 = "(0"
    SA2 = "(0"
    SA3 = "(0"
    SB1 = "(0"
    SB2 = "(0"
    SB3 = "(0"
    for i in q:Q
        SA1 = SA1*" + "*"$(binomial(Q, i))*((1-β)*2*c1A + β*2*c1B)^$(i)*(1.0 - (1-β)*2*c1A - β*2*c1B)^$(Q-i)"
        SA2 = SA2*" + "*"$(binomial(Q, i))*((1-β)*2*(0.5 - c1A - c3A) + β*2*(0.5 - c1B - c3B))^$(i)*(1.0 - (1-β)*2*(0.5 - c1A - c3A) - β*2*(0.5 - c1B - c3B))^$(Q-i)"
        SA3 = SA3*" + "*"$(binomial(Q, i))*((1-β)*2*c3A + β*2*c3B)^$(i)*(1.0 - (1-β)*2*c3A - β*2*c3B)^$(Q-i)"
        SB1 = SB1*" + "*"$(binomial(Q, i))*((1-β)*2*c1B + β*2*c1A)^$(i)*(1.0 - (1-β)*2*c1B - β*2*c1A)^$(Q-i)"
        SB2 = SB2*" + "*"$(binomial(Q, i))*((1-β)*2*(0.5 - c1B - c3B) + β*2*(0.5 - c1A - c3A))^$(i)*(1.0 - (1-β)*2*(0.5 - c1B - c3B) - β*2*(0.5 - c1A - c3A))^$(Q-i)"
        SB3 = SB3*" + "*"$(binomial(Q, i))*((1-β)*2*c3B + β*2*c3A)^$(i)*(1.0 - (1-β)*2*c3B - β*2*c3A)^$(Q-i)"
    end
    SA1 = SA1*")"
    SA2 = SA2*")"
    SA3 = SA3*")"
    SB1 = SB1*")"
    SB2 = SB2*")"
    SB3 = SB3*")"

    def = """
    \tfunction (du, u, params, t)
    \t\tc1A_p = u[1]
    \t\tc3A_p = u[2]
    \t\tc1A_1mp = u[3]
    \t\tc3A_1mp = u[4]

    \t\tc1B_p = u[5]
    \t\tc3B_p = u[6]
    \t\tc1B_1mp = u[7]
    \t\tc3B_1mp = u[8]

    \t\tβ = params[2]
    \t\tp = params[1]

    \t\tc2A_p = p*0.5 - c1A_p - c3A_p
    \t\tc2A_1mp = (1-p)*0.5 - c1A_1mp - c3A_1mp
    \t\tc2B_p = p*0.5 - c1B_p - c3B_p
    \t\tc2B_1mp = (1-p)*0.5 - c1B_1mp - c3B_1mp

    c1A = c1A_p + c1A_1mp
    c2A = c2A_p + c2A_1mp
    c3A = c3A_p + c3A_1mp
    c1B = c1B_p + c1B_1mp
    c2B = c2B_p + c2B_1mp
    c3B = c3B_p + c3B_1mp
   
    \t\tdu[1] = 0.5*c2A_p*$(SA2) - c1A_p*$(SA1)
    \t\tdu[2] = 0.5*c2A_p*$(SA2) - c3A_p*$(SA3)
    \t\tdu[3] = c2A_1mp*$(SA1) - c1A_1mp*$(SA2)
    \t\tdu[4] = c2A_1mp*$(SA3) - c3A_1mp*$(SA2)
    \t\tdu[5] = 0.5*c2B_p*$(SB2) - c1B_p*$(SB1)
    \t\tdu[6] = 0.5*c2B_p*$(SB2) - c3B_p*$(SB3)
    \t\tdu[7] = c2B_1mp*$(SB1) - c1B_1mp*$(SB2)
    \t\tdu[8] = c2B_1mp*$(SB3) - c3B_1mp*$(SB2)
   
    \t\tnothing
    \tend"""
    return eval(Meta.parse(def))
end

function model(q::Int, Q::Int, ::Val{:static2})
    SA1 = "(0"
    SA2 = "(0"
    SA3 = "(0"
    SB1 = "(0"
    SB2 = "(0"
    SB3 = "(0"
    SAleft = "(0"
    SAright = "(0"
    SBleft = "(0"
    SBright = "(0"
    for i in q:Q
        SA1 = SA1*" + "*"$(binomial(Q, i))*((1-β)*2*c1A + β*2*c1B)^$(i)*(1 - (1-β)*2*c1A - β*2*c1B)^$(Q-i)"
        SA2 = SA2*" + "*"$(binomial(Q, i))*((1-β)*2*c2A + β*2*c2B)^$(i)*(1 - (1-β)*2*c2A - β*2*c2B)^$(Q-i)"
        SA3 = SA3*" + "*"$(binomial(Q, i))*((1-β)*2*c3A + β*2*c3B)^$(i)*(1 - (1-β)*2*c3A - β*2*c3B)^$(Q-i)"
        SB1 = SB1*" + "*"$(binomial(Q, i))*((1-β)*2*c1B + β*2*c1A)^$(i)*(1 - (1-β)*2*c1B - β*2*c1A)^$(Q-i)"
        SB2 = SB2*" + "*"$(binomial(Q, i))*((1-β)*2*c2B + β*2*c2A)^$(i)*(1 - (1-β)*2*c2B - β*2*c2A)^$(Q-i)"
        SB3 = SB3*" + "*"$(binomial(Q, i))*((1-β)*2*c3B + β*2*c3A)^$(i)*(1 - (1-β)*2*c3B - β*2*c3A)^$(Q-i)"

        SAleft = SAleft*" + "*"$(binomial(Q, i))*((1-β)*2*(0.5 - c3A) + β*2*(0.5 - c3B))^$(i)*(1 - (1-β)*2*(0.5 - c3A) - β*2*(0.5 - c3B))^$(Q-i)"
        SAright = SAright*" + "*"$(binomial(Q, i))*((1-β)*2*(0.5 - c1A) + β*2*(0.5 - c1B))^$(i)*(1 - (1-β)*2*(0.5 - c1A) - β*2*(0.5 - c1B))^$(Q-i)"
        SBleft = SBleft*" + "*"$(binomial(Q, i))*((1-β)*2*(0.5 - c3B) + β*2*(0.5 - c3A))^$(i)*(1 - (1-β)*2*(0.5 - c3B) - β*2*(0.5 - c3A))^$(Q-i)"
        SBright = SBright*" + "*"$(binomial(Q, i))*((1-β)*2*(0.5 - c1B) + β*2*(0.5 - c1A))^$(i)*(1 - (1-β)*2*(0.5 - c1B) - β*2*(0.5 - c1A))^$(Q-i)"
    end
    SA1 = SA1*")"
    SA2 = SA2*")"
    SA3 = SA3*")"
    SB1 = SB1*")"
    SB2 = SB2*")"
    SB3 = SB3*")"
    SAleft = SAleft*")"
    SAright = SAright*")"
    SBleft = SBleft*")"
    SBright = SBright*")"

    def = """
    \tfunction (du, u, params, t)
    \t\tc1A_p = u[1]
    \t\tc3A_p = u[2]
    \t\tc1A_1mp = u[3]
    \t\tc3A_1mp = u[4]

    \t\tc1B_p = u[5]
    \t\tc3B_p = u[6]
    \t\tc1B_1mp = u[7]
    \t\tc3B_1mp = u[8]

    \t\tβ = params[2]
    \t\tp = params[1]

    \t\tc2A_p = p*0.5 - c1A_p - c3A_p
    \t\tc2A_1mp = (1-p)*0.5 - c1A_1mp - c3A_1mp
    \t\tc2B_p = p*0.5 - c1B_p - c3B_p
    \t\tc2B_1mp = (1-p)*0.5 - c1B_1mp - c3B_1mp
    
    c1A = c1A_p + c1A_1mp
    c2A = c2A_p + c2A_1mp
    c3A = c3A_p + c3A_1mp
    c1B = c1B_p + c1B_1mp
    c2B = c2B_p + c2B_1mp
    c3B = c3B_p + c3B_1mp
    
    \t\tdu[1] = 0.5*c2A_p*$(SA2) - c1A_p*$(SA1)
    \t\tdu[2] = 0.5*c2A_p*$(SA2) - c3A_p*$(SA3)
    \t\tdu[3] = c2A_1mp*$(SA1) - c1A_1mp*$(SAright)
    \t\tdu[4] = c2A_1mp*$(SA3) - c3A_1mp*$(SAleft)
    \t\tdu[5] = 0.5*c2B_p*$(SB2) - c1B_p*$(SB1)
    \t\tdu[6] = 0.5*c2B_p*$(SB2) - c3B_p*$(SB3)
    \t\tdu[7] = c2B_1mp*$(SB1) - c1B_1mp*$(SBright)
    \t\tdu[8] = c2B_1mp*$(SB3) - c3B_1mp*$(SBleft)
   
    \t\tnothing
    \tend"""
    return eval(Meta.parse(def))
end

function model(q::Int, Q::Int, ::Val{:static3})
    SA1 = "(0"
    SA2 = "(0"
    SA3 = "(0"
    SB1 = "(0"
    SB2 = "(0"
    SB3 = "(0"
    SAleft = "(0"
    SAright = "(0"
    SBleft = "(0"
    SBright = "(0"
    for i in q:Q
        SA1 = SA1*" + "*"$(binomial(Q, i))*((1-β)*2*c1A + β*2*c1B)^$(i)*(1 - (1-β)*2*c1A - β*2*c1B)^$(Q-i)"
        SA2 = SA2*" + "*"$(binomial(Q, i))*((1-β)*2*c2A + β*2*c2B)^$(i)*(1 - (1-β)*2*c2A - β*2*c2B)^$(Q-i)"
        SA3 = SA3*" + "*"$(binomial(Q, i))*((1-β)*2*c3A + β*2*c3B)^$(i)*(1 - (1-β)*2*c3A - β*2*c3B)^$(Q-i)"
        SB1 = SB1*" + "*"$(binomial(Q, i))*((1-β)*2*c1B + β*2*c1A)^$(i)*(1 - (1-β)*2*c1B - β*2*c1A)^$(Q-i)"
        SB2 = SB2*" + "*"$(binomial(Q, i))*((1-β)*2*c2B + β*2*c2A)^$(i)*(1 - (1-β)*2*c2B - β*2*c2A)^$(Q-i)"
        SB3 = SB3*" + "*"$(binomial(Q, i))*((1-β)*2*c3B + β*2*c3A)^$(i)*(1 - (1-β)*2*c3B - β*2*c3A)^$(Q-i)"

        SAleft = SAleft*" + "*"$(binomial(Q, i))*((1-β)*2*(0.5 - c3A) + β*2*(0.5 - c3B))^$(i)*(1 - (1-β)*2*(0.5 - c3A) - β*2*(0.5 - c3B))^$(Q-i)"
        SAright = SAright*" + "*"$(binomial(Q, i))*((1-β)*2*(0.5 - c1A) + β*2*(0.5 - c1B))^$(i)*(1 - (1-β)*2*(0.5 - c1A) - β*2*(0.5 - c1B))^$(Q-i)"
        SBleft = SBleft*" + "*"$(binomial(Q, i))*((1-β)*2*(0.5 - c3B) + β*2*(0.5 - c3A))^$(i)*(1 - (1-β)*2*(0.5 - c3B) - β*2*(0.5 - c3A))^$(Q-i)"
        SBright = SBright*" + "*"$(binomial(Q, i))*((1-β)*2*(0.5 - c1B) + β*2*(0.5 - c1A))^$(i)*(1 - (1-β)*2*(0.5 - c1B) - β*2*(0.5 - c1A))^$(Q-i)"
    end
    SA1 = SA1*")"
    SA2 = SA2*")"
    SA3 = SA3*")"
    SB1 = SB1*")"
    SB2 = SB2*")"
    SB3 = SB3*")"
    SAleft = SAleft*")"
    SAright = SAright*")"
    SBleft = SBleft*")"
    SBright = SBright*")"
    
    def = """
    \tfunction (du, u, params, t)
    \t\tc1A_p = u[1]
    \t\tc3A_p = u[2]
    \t\tc1A_1mp = u[3]
    \t\tc3A_1mp = u[4]

    \t\tc1B_p = u[5]
    \t\tc3B_p = u[6]
    \t\tc1B_1mp = u[7]
    \t\tc3B_1mp = u[8]

    \t\tβ = params[2]
    \t\tp = params[1]

    \t\tc2A_p = p*0.5 - c1A_p - c3A_p
    \t\tc2A_1mp = (1-p)*0.5 - c1A_1mp - c3A_1mp
    \t\tc2B_p = p*0.5 - c1B_p - c3B_p
    \t\tc2B_1mp = (1-p)*0.5 - c1B_1mp - c3B_1mp

    c1A = c1A_p + c1A_1mp
    c2A = c2A_p + c2A_1mp
    c3A = c3A_p + c3A_1mp
    c1B = c1B_p + c1B_1mp
    c2B = c2B_p + c2B_1mp
    c3B = c3B_p + c3B_1mp
    
    \t\tdu[1] = c2A_p*$(SA1) - c1A_p*$(SAright)
    \t\tdu[2] = c2A_p*$(SA3) - c3A_p*$(SAleft)
    \t\tdu[3] = c2A_1mp*$(SA1) - c1A_1mp*$(SA2)
    \t\tdu[4] = c2A_1mp*$(SA3) - c3A_1mp*$(SA2)
    \t\tdu[5] = c2B_p*$(SB1) - c1B_p*$(SBright)
    \t\tdu[6] = c2B_p*$(SB3) - c3B_p*$(SBleft)
    \t\tdu[7] = c2B_1mp*$(SB1) - c1B_1mp*$(SB2)
    \t\tdu[8] = c2B_1mp*$(SB3) - c3B_1mp*$(SB2)
   
    \t\tnothing
    \tend"""
    return eval(Meta.parse(def))
end
