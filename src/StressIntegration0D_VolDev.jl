@doc raw"""
    results = Vermeer1990_StressIntegration_vdev(σi; params)  ;

Function that computes a volumetric-deviatoric return mapping for test 1 in Vermeer (1990). 
The function takes as input:

    σi : a tuple of initial stress, e.g.: `σi = (xx=-25e3, yy=-100e3)`
and returns:

    results : all the data needed for the figure (apparent friction, horizontal stress, volume change and stress orientation)

# Examples
```julia-repl
julia>  results = Vermeer1990_StressIntegration_vdev( (xx=-25e3, yy=-100e3))

```
"""
function Vermeer1990_StressIntegration_vdev(σi; params=(
    K   = 6.6666666667e6, # K = 3/2*Gv in Vermeer (1990)
    G   = 10e6,
    c   = 0.0,
    ϕ   = 40/180*π,
    ψ   = 10/180*π,
    θt  = 25/180*π,
    ηvp = 0.,
    γ̇xy = 0.00001,
    Δt  = 20,
    nt  = 400,
    law = :MC_Vermeer1990, 
    oop = :Vermeer1990,
    pl  = true) # default parameter set
)

    @info "Stress integration: volumetric-deviatoric split v2 with $(params.law)"
    law = params.law
    _type(law) = Val{law}()
    yield(  τ, λ̇, ϕ, c, ηvp, θt, law) = _yield( τ, λ̇, ϕ, c, ηvp, θt, _type(law)) 
    _yield( τ, λ̇, ϕ, c, ηvp, θt, ::Val{:MC_Vermeer1990}) = MohrCoulombVermeer1990_vdev(τ, λ̇, ϕ, c, ηvp, θt)
    _yield( τ, λ̇, ϕ, c, ηvp, θt, ::Val{:DruckerPrager})  =          DruckerPrager_vdev(τ, λ̇, ϕ, c, ηvp, θt)
    _yield( τ, λ̇, ϕ, c, ηvp, θt, ::Val{:MC_AS95})        =       MohrCoulomb_AS95_vdev(τ, λ̇, ϕ, c, ηvp, θt)
    _yield( τ, λ̇, ϕ, c, ηvp, θt, ::Val{:MC_deBorst90})   =  MohrCoulomb_deBorst90_vdev(τ, λ̇, ϕ, c, ηvp, θt)
    
    # Material properties
    Kb  = params.K
    G   = params.G
    c   = params.c
    ϕ   = params.ϕ
    ψ   = params.ψ
    θt  = params.θt
    ηvp = params.ηvp

    # Strain integration
    ε̇xy   = params.γ̇xy/2
    Δt    = params.Δt 
    nt    = params.nt 
    niter = 1000

    # Initial condition
    P    = -1/2*(σi.xx+σi.yy)
    τxxi = σi.xx + P
    τyyi = σi.yy + P
    τxx  = τxxi 
    τyy  = τyyi
    τzz  = 0.
    τxy  = 0. 
    εyy  = 0.
    εxy  = 0.

    # Storage
    τ_vec       = zeros(5)
    app_fric    = zeros(nt)
    γxyv        = zeros(nt)
    θv          = zeros(nt)
    σxxv        = zeros(nt)
    εyyv        = zeros(nt)
    σ1v         = zeros(nt)
    σ2v         = zeros(nt)
    σ3v         = zeros(nt)
    θlv         = zeros(nt)
    app_fric[1] = -τxy/(τyy-P)

    # Time integration loop
    for it=2:nt

        # History
        τxx0 = τxx
        τyy0 = τyy
        τzz0 = τzz
        τxy0 = τxy
        P0   = P
        λ̇    = 0.

        # σ3 angle
        σm  = [τxx-P τxy; τxy τyy-P]
        v    = eigvecs(σm)        
        θ    = atan(v[2,2] / v[1,2])   # σ3 angle
    
        # Trial stress integration
        ε̇xxd_pl, ε̇yyd_pl, ε̇zzd_pl, ε̇xy_pl, ∇v_pl = 0., 0., 0., 0., 0.
        
        # Total strain rates
        ε̇xx = 0.
        ε̇yy = 0. # to be found
        ε̇yy = (4.0 * G .* Δt .* ε̇yyd_pl - Kb .* Δt .* ε̇xx + 2.0 * Kb .* Δt .* ∇v_pl + 2.0 * P0 + 2.0 * σi.yy - 2.0 * τyy0) ./ (Δt .* (4.0 * G + 3.0 * Kb))
        ε̇zz = 1/2*(ε̇xx + ε̇yy)
        ∇v  = ε̇xx + ε̇yy + ε̇zz

        # Deviatoric strain rates
        ε̇xxd = ε̇xx - 1/3*∇v
        ε̇yyd = ε̇yy - 1/3*∇v
        ε̇zzd = ε̇zz - 1/3*∇v
        
        # Deviatoric stress and pressure
        τxx = τxx0 + 2*G*Δt*(ε̇xxd)
        τyy = τyy0 + 2*G*Δt*(ε̇yyd)
        τzz = τzz0 + 2*G*Δt*(ε̇zzd)
        τxy = τxy0 + 2*G*Δt*ε̇xy
        P   = P0   - Kb*Δt*∇v

        # Yield criteria and plastic flow direction
        τ_vec .= [τxx; τyy; τzz; τxy; P]
        f      = yield( τ_vec, λ̇, ϕ, c, ηvp, θt, law)
        q      = τ -> yield(τ, λ̇, ψ, c, ηvp, θt, law)
        # dqdτ   = gradient(Forward, q, τ_vec )
        dqdτ   = ForwardDiff.gradient(q, τ_vec )
        fc     = f

        # Plastic correction: return mapping
        if f>0 && params.pl
            dQdτxx = dqdτ[1]
            dQdτyy = dqdτ[2]
            dQdτzz = dqdτ[3]
            dQdτxy = dqdτ[4]
            dQdP   = dqdτ[5]
            λ̇      = 0.0
             
            for iter=1:niter

                # Plastic total strain rates
                ε̇xy_pl  =  λ̇*dQdτxy/2
                ε̇xxd_pl =  λ̇*dQdτxx 
                ε̇yyd_pl =  λ̇*dQdτyy
                # with this one: ezzp'!=0, ezzp=0, tzz!=0, ezz'=0
                # ∇vp = ε̇xxp + ε̇yyp + ε̇zzp so ideally ε̇zzp = 0 such that ∇vp = -λ*dQdP
                ε̇zzd_pl =  λ̇*dQdτzz
                ∇v_pl   = -λ̇*dQdP
            
                if params.oop == :Vermeer1990
                    # with this one: ezzp'=0, ezzp!=0, tzz=0, ezz'=0
                    # this allows every return mapping to match on case A and B of Vermeer (1990)
                    ε̇zzd_pl =  1/2*λ̇*(dQdτxx + dQdτyy)
                    ∇v_pl   = -3/2*λ̇*dQdP
                end

                # Total strain rates
                ε̇xx = 0.      
                ε̇yy = (4.0 * G .* Δt .* ε̇yyd_pl - Kb .* Δt .* ε̇xx + 2.0 * Kb .* Δt .* ∇v_pl + 2.0 * P0 + 2.0 * σi.yy - 2.0 * τyy0) ./ (Δt .* (4.0 * G + 3.0 * Kb))
                ε̇zz = 1/2*(ε̇xx + ε̇yy)
                ∇v  = ε̇xx + ε̇yy + ε̇zz

                # Deviatoric strain rates (recompute because ∇v changed just above)
                ε̇xxd = ε̇xx - 1/3*∇v
                ε̇yyd = ε̇yy - 1/3*∇v
                ε̇zzd = ε̇zz - 1/3*∇v

                # Make sure to make to make plastic strain deviatoric for the case of Vermeer's yield surface
                if params.law == :MC_Vermeer1990
                    ε̇xxd_pl = ε̇xxd_pl - 1/3*∇v_pl
                    ε̇yyd_pl = ε̇yyd_pl - 1/3*∇v_pl
                    ε̇zzd_pl = ε̇zzd_pl - 1/3*∇v_pl
                end
               
                # Deviatoric stress and pressure
                τxx = τxx0 + 2*G*Δt*(ε̇xxd - ε̇xxd_pl)
                τyy = τyy0 + 2*G*Δt*(ε̇yyd - ε̇yyd_pl)
                τzz = τzz0 + 2*G*Δt*(ε̇zzd - ε̇zzd_pl)
                τxy = τxy0 + 2*G*Δt*(ε̇xy  - ε̇xy_pl )
                P   = P0   -  Kb*Δt*(∇v   - ∇v_pl  )

                # if iter==2 && it==30
                #     @show dQdτzz, dQdP
                #     @show (ε̇xxd_pl + ε̇yyd_pl + ε̇zzd_pl)
                #     @show ∇v_pl
                #     @show ε̇zzd_pl
                #     @show (ε̇xxd_pl+ε̇yyd_pl)
                # end
                
                # Yield criteria
                τ_vec .= [τxx; τyy; τzz; τxy; P]
                fc     = yield(τ_vec, λ̇, ϕ, c, ηvp, θt, law)
                λ̇      += fc / ((Kb + G)*Δt) /5
                if abs(fc)<1e-8 break end
                # @show f, fc
                if iter==niter error("Failed return mapping: f = $(f) --- fc = $(fc)") end
            end
        end

        # Increment total strains
        εyy += ε̇yy*Δt
        εxy += ε̇xy*Δt

        # Principal stresses
        σm  = [τxx-P τxy 0; τxy τyy-P 0; 0 0 τzz-P]
        v   = eigvals(σm) 

        # Lode
        τ   = [τxx; τyy; τzz; τxy; P]
        J2  = 0.5*(τ[1]^2 + τ[2]^2 + τ[3]^2) + τ[4]^2
        J3  = τ[1]*τ[2]*τ[3] + τ[3]*τ[4]^2
        θl  = -3/2*sqrt(3)*J3/J2/sqrt(J2)

        # Storage
        app_fric[it] = -τxy/(τyy - P)
        γxyv[it]     = 2.0.*εxy
        εyyv[it]     = εyy
        θv[it]       = θ
        σxxv[it]     = τxx - P
        σ1v[it]      = v[1]
        σ2v[it]      = v[2]
        σ3v[it]      = v[3]
        θlv[it]      = θl
    end
    return (γxy=γxyv[2:end].*100, εyy=εyyv[2:end].*100, app_fric=app_fric[2:end], σxx=.-σxxv[2:end]./1e3, θ=θv[2:end].*(180/π), σ1=σ1v, σ2=σ2v, σ3=σ3v, θl=θlv)
end