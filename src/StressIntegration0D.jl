using Printf, LinearAlgebra, Enzyme
import ForwardDiff

@doc raw"""
    results = Vermeer90_StressIntegration_tot(σi; params)  ;

Function that computes a total stress return mapping for test 1 in Vermeer (1990). 
The function takes as input:
    σi : a tuple of initial stress, e.g.: `σi = (xx=-25e3, yy=-100e3)`
and returns:
    results : all the data needed for the figure (apparent friction, horizontal stress, volume change and stress orientation)

# Examples
```julia-repl
julia>  results = Vermeer90_StressIntegration_tot( (xx=-25e3, yy=-100e3))

```
"""
function Vermeer90_StressIntegration_tot(σi; params=(
        G   = 10e6,
        c   = 0.,
        ϕ   = 40/180*π,
        ψ   = 10/180*π,
        θt  = 25/180*π,
        ηvp = 0.,
        γ̇xy = 0.00001,
        Δt  = 20,
        nt  = 400,
        law = :MC_Vermeer1990,
        pl  = true) # default parameter set
    )

    @info "Stress integration: total stress"
    law = params.law
    _type(law) = Val{law}()
    yield(  τ, λ̇, ϕ, c, ηvp, θt, law) = _yield( τ, λ̇, ϕ, c, ηvp, θt, _type(law)) 
    _yield( τ, λ̇, ϕ, c, ηvp, θt, ::Val{:MC_Vermeer1990}) = MohrCoulombVermeer1990_tot(τ, λ̇, ϕ, c, ηvp, θt)
    _yield( τ, λ̇, ϕ, c, ηvp, θt, ::Val{:DruckerPrager})  =          DruckerPrager_tot(τ, λ̇, ϕ, c, ηvp, θt)

    # Material properties
    G   = params.G
    c   = params.c
    ϕ   = params.ϕ
    ψ   = params.ψ
    θt  = params.θt
    ηvp = params.ηvp
    
    # Strain integration
    γ̇xy = params.γ̇xy
    Δt  = params.Δt 
    nt  = params.nt 

    # Initial condition
    σxx = σi.xx 
    σyy = σi.yy
    σxy = 0.  
    σzz = 0.5*(σi.xx + σi.yy)
    σm  = [σxx σxy; σxy σyy]
    v   = eigvecs(σm)        
    θ   = atan(v[2,2] / v[1,2])   # σ3 angle
    ε̇xx = 0.
    ε̇yy = 0. # to be found
    εyy = 0.
    γxy = 0.
    
    # Storage
    σ_vec       = zeros(4)
    app_fric    = zeros(nt)
    σxxv        = zeros(nt) 
    εyyv        = zeros(nt) 
    γxyv        = zeros(nt)
    θv          = zeros(nt)
    app_fric[1] = -σxy/σyy

    for it=2:nt

        # History
        σxx0 = σxx
        σyy0 = σyy
        σzz0 = σzz
        σxy0 = σxy
        λ̇    = 0.

        # σ3 angle
        σm   = [σxx σxy; σxy σyy]
        v    = eigvecs(σm)        
        θ    = atan(v[2,2] / v[1,2])   # σ3 angle

        # Trial stress integration
        ε̇yy  = (σyy0 - σi.yy) / (2*G*Δt)
        σxx  = σxx0 + 2*G*Δt*ε̇xx
        σyy  = σyy0 + 2*G*Δt*ε̇yy # this will not change
        σxy  = σxy0 +   G*Δt*γ̇xy
        ε̇zz  = 0. + (σzz0 - σzz) / (2*G*Δt)
        σzz  = σzz0 + 2*G*Δt*ε̇zz # σzz  = (σxx + σyy)/2

        # Yield criteria and plastic flow direction
        σ_vec .= [σxx; σyy; σzz; σxy;]
        f      = yield(σ_vec, λ̇, ϕ, c, ηvp, θt, law)
        q      = σ -> yield(σ, λ̇, ψ, c, ηvp, θt, law)
        dqdσ   = ForwardDiff.gradient(q, σ_vec)
        # dqdσ   = gradient(Forward, q, σ_vec)
        fc     = f

        # Plastic correction: return mapping
        if f>0 && params.pl
            λ = 0.0
            for iter=1:1000

                # Plastic total strains
                ε̇xxp   = λ̇*dqdσ[1]
                ε̇yyp   = λ̇*dqdσ[2]
                ε̇zzp   = λ̇*dqdσ[3]
                γ̇xyp   = λ̇*dqdσ[4]

                # Total strains
                ε̇yy    = ε̇yyp + (σyy0 - σi.yy) / (2*G*Δt)
                ε̇zz    = 1/2*(ε̇xx + ε̇yy)

                # Total stresses
                σxx    = σxx0 + 2*G*Δt*(ε̇xx - ε̇xxp)
                σyy    = σyy0 + 2*G*Δt*(ε̇yy - ε̇yyp) # this will not change
                σxy    = σxy0 +   G*Δt*(γ̇xy - γ̇xyp)
                σzz    = σzz0 + 2*G*Δt*(ε̇zz - ε̇zzp)
                
                # Yield criteria
                σ_vec .= [σxx; σyy; σzz; σxy;]
                fc     = yield(σ_vec, λ̇, ϕ, c, ηvp, θt, law)

                λ̇     += fc / G / 20
                if abs(fc)<1e-8 break end
            end
        end

        # Increment total strains
        εyy += ε̇yy*Δt
        γxy += γ̇xy*Δt

        # Storage
        app_fric[it] = -σxy/σyy
        γxyv[it]     = γxy
        εyyv[it]     = εyy
        θv[it]       = θ
        σxxv[it]     = σxx
    end
    return (γxy=γxyv[2:end].*100, εyy=εyyv[2:end].*100, app_fric=app_fric[2:end], σxx=.-σxxv[2:end]./1e3, θ=θv[2:end].*(180/π))
end

@doc raw"""
    results = Vermeer90_StressIntegration_vdev(σi; params)  ;

Function that computes a volumetric-deviatoric return mapping for test 1 in Vermeer (1990). 
The function takes as input:

    σi : a tuple of initial stress, e.g.: `σi = (xx=-25e3, yy=-100e3)`
and returns:

    results : all the data needed for the figure (apparent friction, horizontal stress, volume change and stress orientation)

# Examples
```julia-repl
julia>  results = Vermeer90_StressIntegration_vdev( (xx=-25e3, yy=-100e3))

```
"""
function Vermeer90_StressIntegration_vdev(σi; params=(
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
    pl  = true) # default parameter set
)

    @info "Stress integration: volumetric-deviatoric split"
    law = params.law
    _type(law) = Val{law}()
    yield(  τ, λ̇, ϕ, c, ηvp, θt, law) = _yield( τ, λ̇, ϕ, c, ηvp, θt, _type(law)) 
    _yield( τ, λ̇, ϕ, c, ηvp, θt, ::Val{:MC_Vermeer1990}) = MohrCoulombVermeer1990_vdev(τ, λ̇, ϕ, c, ηvp, θt)
    _yield( τ, λ̇, ϕ, c, ηvp, θt, ::Val{:DruckerPrager})  =          DruckerPrager_vdev(τ, λ̇, ϕ, c, ηvp, θt)
    _yield( τ, λ̇, ϕ, c, ηvp, θt, ::Val{:MC_AS95})        =       MohrCoulomb_AS95_vdev(τ, λ̇, ϕ, c, ηvp, θt)

    # Material properties
    G   = params.G
    c   = params.c
    ϕ   = params.ϕ
    ψ   = params.ψ
    θt  = params.θt
    ηvp = params.ηvp

    # Strain integration
    ε̇xy = params.γ̇xy/2
    Δt  = params.Δt 
    nt  = params.nt 

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
        ε̇xxp, ε̇yyp, ε̇zzp, ε̇xyp, ∇vp = 0., 0., 0., 0., 0.
        
        # Total strain rates
        ε̇xx = 0.
        ε̇yy = 0. # to be found
        ε̇yy = (0.333333333333333 * G .* Δt .* (-ε̇xxp + 2.0 * ε̇yyp - ε̇zzp + ∇vp) + 0.5 * P0 + 0.5 * σi.yy - 0.5 * τyy0) ./ (G .* Δt)
        ε̇zz = 1/2*(ε̇xx + ε̇yy)
        ∇v  = ε̇xx + ε̇yy + ε̇zz

        # Deviatoric strain rates
        ε̇xxd = ε̇xx - 1/3*(ε̇xx + ε̇yy + ε̇zz)
        ε̇yyd = ε̇yy - 1/3*(ε̇xx + ε̇yy + ε̇zz)
        ε̇zzd = ε̇zz - 1/3*(ε̇xx + ε̇yy + ε̇zz)
        
        # Deviatoric stress and pressure
        τxx = τxx0 + 4/3*G*Δt*(ε̇xxd + 1/3*∇v) - 2/3*G*Δt*(ε̇yyd + 1/3*∇v) - 2/3*G*Δt*(ε̇zzd + 1/3*∇v)
        τyy = τyy0 + 4/3*G*Δt*(ε̇yyd + 1/3*∇v) - 2/3*G*Δt*(ε̇xxd + 1/3*∇v) - 2/3*G*Δt*(ε̇zzd + 1/3*∇v)
        τzz = τzz0 + 4/3*G*Δt*(ε̇zzd + 1/3*∇v) - 2/3*G*Δt*(ε̇xxd + 1/3*∇v) - 2/3*G*Δt*(ε̇yyd + 1/3*∇v)
        τxy = τxy0 +   2*G*Δt*ε̇xy
        P   = P0   - 2/3*G*Δt*∇v

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
            for iter=1:1000

                # Plastic total strain rates
                ε̇xyp =  λ̇*dQdτxy/2
                ε̇xxp =  λ̇*dQdτxx 
                ε̇yyp =  λ̇*dQdτyy
                # with this one: ezzp'=0, ezzp!=0, tzz=0, ezz'=0
                # ε̇zzp =  1/2*(ε̇xxp + ε̇yyp)
                # ∇vp  = -3/2*λ̇*dQdP
                # with this one: ezzp'!=0, ezzp=0, tzz!=0, ezz'=0
                ε̇zzp =  λ̇*dQdτzz
                ∇vp  = -λ̇*dQdP
                # ∇vp = ε̇xxp + ε̇yyp + ε̇zzp so ideally ε̇zzp = 0 such that ∇vp = -λ*dQdP

                # Total strain rates
                ε̇xx  = 0.
                ε̇yy  = (0.333333333333333 * G .* Δt .* (-ε̇xxp + 2.0 * ε̇yyp - ε̇zzp + ∇vp) + 0.5 * P0 + 0.5 * σi.yy - 0.5 * τyy0) ./ (G .* Δt)
                ε̇zz  = 1/2*(ε̇xx + ε̇yy)
                ∇v   = ε̇xx + ε̇yy + ε̇zz

                # Deviatoric strain rates
                ε̇xxd = ε̇xx - 1/3*∇v
                ε̇yyd = ε̇yy - 1/3*∇v
                ε̇zzd = ε̇zz - 1/3*∇v

                # Deviatoric stress and pressure
                τxx = τxx0 + 4/3*G*Δt*(ε̇xxd - ε̇xxp + 1/3*∇vp) - 2/3*G*Δt*(ε̇yyd - ε̇yyp + 1/3*∇vp) - 2/3*G*Δt*(ε̇zzd - ε̇zzp + 1/3*∇vp)
                τyy = τyy0 - 2/3*G*Δt*(ε̇xxd - ε̇xxp + 1/3*∇vp) + 4/3*G*Δt*(ε̇yyd - ε̇yyp + 1/3*∇vp) - 2/3*G*Δt*(ε̇zzd - ε̇zzp + 1/3*∇vp)
                τzz = τzz0 - 2/3*G*Δt*(ε̇xxd - ε̇xxp + 1/3*∇vp) - 2/3*G*Δt*(ε̇yyd - ε̇yyp + 1/3*∇vp) + 4/3*G*Δt*(ε̇zzd - ε̇zzp + 1/3*∇vp)
                τxy = τxy0 +   2*G*Δt*(ε̇xy  - ε̇xyp)
                P   = P0   - 2/3*G*Δt*(∇v- ∇vp) 
                
                # Yield criteria
                τ_vec .= [τxx; τyy; τzz; τxy; P]
                fc     = yield(τ_vec, λ̇, ϕ, c, ηvp, θt, law)
                λ̇     += fc / G / 20
                if abs(fc)<1e-8 break end
            end
        end

        # Increment total strains
        εyy += ε̇yy*Δt
        εxy += ε̇xy*Δt

        # Storage
        app_fric[it] = -τxy/(τyy - P)
        γxyv[it]     = 2.0.*εxy
        εyyv[it]     = εyy
        θv[it]       = θ
        σxxv[it]     = τxx - P
    end
    return (γxy=γxyv[2:end].*100, εyy=εyyv[2:end].*100, app_fric=app_fric[2:end], σxx=.-σxxv[2:end]./1e3, θ=θv[2:end].*(180/π))
end