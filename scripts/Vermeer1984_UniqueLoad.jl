function MainVermeer1984()

    σyyB       = -100e3 # Courbe A - Vermeer
    σxxB       =  -25e3 # Courbe A - Vermeer
    σxxB       = -400e3 # Courbe B - Vermeer
    σyyB       = -100e3 # Courbe B - Vermeer
    σzzB       = σxxB
    PB         = -(σxxB + σyyB + σzzB)/3.0
    τxxB       = PB + σxxB
    τyyB       = PB + σyyB
    τzzB       = PB + σzzB
    τxyB       = 0.0

    ε̇xx        = 0.0
    ε̇yy        = 0.0
    ε̇zz        = 0.0
    ε̇xy        = 1e-14

    Δt         = 1e10
    nt         = 300
    E          = 45e6
    ν          = 0.2
    G          = E/2.0/(1+ν)
    Kb         = E/3.0/(1-2ν)
    ϕ          = 43*π/180
    ψ          = 15*π/180
    c          = 0

    τxx = τxxB
    τyy = τyyB
    τzz = τzzB
    τxy = τxyB
    P   = PB
    γxy = 0

    hist = (τxy=zeros(nt), γxy=zeros(nt))

    for it=1:nt
        τxx0 = τxx
        τyy0 = τyy
        τzz0 = τzz
        τxy0 = τxy
        P0   = P

        ∇v   = ε̇xx + ε̇yy + ε̇zz
        τxx  = τxx0 + 2*G*Δt*(ε̇xx-1/3*∇v)
        τyy  = τyy0 + 2*G*Δt*(ε̇yy-1/3*∇v)
        τzz  = τzz0 + 2*G*Δt*(ε̇zz-1/3*∇v)
        τxy  = τxy0 + 2*G*Δt*(ε̇xy-1/3*∇v)
        P    = P0   -   Kb*Δt*∇v

        τII  = sqrt(0.5*(τxx^2 + τyy^2 + τzz^2) + τxy^2)
        F    = τII - P*sin(ϕ) - c*cos(ϕ)
        @show F
        if F>0
            λ̇    = F / (G*Δt + Kb*Δt*sin(ϕ)*sin(ψ))
            τxx  = τxx0 + 2*G*Δt*(ε̇xx - 1/3*∇v - τxx/2/τII*λ̇)
            τyy  = τyy0 + 2*G*Δt*(ε̇yy - 1/3*∇v - τyy/2/τII*λ̇)
            τzz  = τzz0 + 2*G*Δt*(ε̇zz - 1/3*∇v - τzz/2/τII*λ̇)
            τxy  = τxy0 + 2*G*Δt*(ε̇xy - 1/3*∇v - τxy/2/τII*λ̇)
            P    = P0   - Kb*Δt*(∇v - sin(ψ)*λ̇)
            τII  = sqrt(0.5*(τxx^2 + τyy^2 + τzz^2) + τxy^2)
            F    = τII - P*sin(ϕ) - c*cos(ϕ)
            @show F
        end

        γxy += ε̇xy*Δt
        hist.τxy[it] = τxy
        hist.γxy[it] = γxy
    end
    plot(hist.γxy,  hist.τxy/1e3)
end

MainVermeer1984()