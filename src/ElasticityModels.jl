function ElasticModel( Kb, G, Δt, params )
    ηn, ηs, ηb = 1.0, 1.0, 1.0
    if params.el == :Vermeer1990
        ηn = 2/3*G*Δt
        ηs = 1/3*G*Δt
        ηb = 2/3*G*Δt
    elseif params.el == :Standard
        ηn = G*Δt
        ηs = 0.
        ηb = Kb*Δt
    end
    return ηn, ηs, ηb
end
