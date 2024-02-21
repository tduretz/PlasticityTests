@doc raw"""
    results = Vermeer1990_StressIntegration_tot(σi; params)  ;

Function that computes a total stress return mapping for test 1 in Vermeer (1990). 
The function takes as input:
    σi : a tuple of initial stress, e.g.: `σi = (xx=-25e3, yy=-100e3)`
and returns:
    results : all the data needed for the figure (apparent friction, horizontal stress, volume change and stress orientation)

# Examples
```julia-repl
julia>  results = Vermeer1990_StressIntegration_tot( (xx=-25e3, yy=-100e3))

```
"""
function Vermeer1990_StressIntegration_tot(σi; params=(
        K   = 6.6666666667e6, # K = 3/2*Gv in Vermeer (1990)
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
        oop = :Vermeer1990,
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
    γ̇xy   = params.γ̇xy
    Δt    = params.Δt 
    nt    = params.nt
    niter = 1000 

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
            for iter=1:niter

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
                if iter==niter error("Failed return mapping") end
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