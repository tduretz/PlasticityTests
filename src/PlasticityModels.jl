@doc raw"""
    f = MohrCoulombVermeer1990_vdev(τ, λ̇, ϕ, c, ηvp)  

Computes Mohr-Coulomb yield function for test 1 in Vermeer (1990) using deviatoric stresses and pressure. 
The function takes as input:

    τ   : vector of deviatoric stress components and pressure 
    λ̇   : plastic multiplier rate
    ϕ   : friction angle
    c   : cohesion
    ηvp : viscoplastic viscosity
and returns:

    f   : yield function value

# Examples
```julia-repl
julia>  MohrCoulombVermeer1990_vdev(τ, λ̇, ϕ, c, ηvp)

```
"""
function MohrCoulombVermeer1990_vdev(τ, λ̇, ϕ, c, ηvp)
    # τII    = sqrt(0.25*(τxx-τyy)^2 + τxy^2)
    # dQdτxx =  0.25*(τxx-τyy)/τII + 0.5*sin(ψ)
    # dQdτyy = -0.25*(τxx-τyy)/τII + 0.5*sin(ψ)
    # dQdτzz = 0
    # dQdτxy =  τxy/τII/2
    # dQdP   = -sin(ψ)
    τII  = sqrt(0.25*(τ[1]-τ[2])^2 + τ[4]^2)
    f    = τII + 0.5*(τ[1]+τ[2]-2*τ[5])*sin(ϕ) - c*cos(ϕ) - λ̇*ηvp
    return f
end

@doc raw"""
    f = MohrCoulombVermeer1990_tot(σ, λ̇, ϕ, c, ηvp)  

Computes Mohr-Coulomb yield function for test 1 in Vermeer (1990) using total stresses. 
The function takes as input:

    σ   : vector of total stress components
    λ̇   : plastic multiplier rate
    ϕ   : friction angle
    c   : cohesion
    ηvp : viscoplastic viscosity
and returns:

    f   : yield function value

# Examples
```julia-repl
julia>  MohrCoulombVermeer1990_tot(σ, λ̇, ϕ, c, ηvp)

```
"""
function MohrCoulombVermeer1990_tot(σ, λ̇, ϕ, c, ηvp)
    # τII    = sqrt(0.25*(σxx-σyy)^2 + σxy^2)
    # dQdσxx = 0.25*(σxx-σyy)/τII + 0.5*sin(ψ)
    # dQdσyy =-0.25*(σxx-σyy)/τII + 0.5*sin(ψ)
    # dQdσzz = 0.0
    # dQdσxy = σxy/τII
    τII  = sqrt(0.25*(σ[1]-σ[2])^2 + σ[4]^2)
    f    = τII + 0.5*(σ[1]+σ[2])*sin(ϕ) - c*cos(ϕ) - λ̇*ηvp
    return f
end

@doc raw"""
    f = DruckerPrager_vdev(τ, λ̇, ϕ, c, ηvp)  

Computes Drucker-Prager yield function for test 1 in Vermeer (1990) using deviatoric stresses and pressure. 
The function takes as input:

    τ   : vector of deviatoric stress components and pressure 
    λ̇   : plastic multiplier rate
    ϕ   : friction angle
    c   : cohesion
    ηvp : viscoplastic viscosity
and returns:

    f   : yield function value

# Examples
```julia-repl
julia>  DruckerPrager_vdev(τ, λ̇, ϕ, c, ηvp)

```
"""
function DruckerPrager_vdev(τ, λ̇, ϕ, c, ηvp)
    # τII    = sqrt(0.25*(τxx-τyy)^2 + τxy^2)
    # dQdτxx =  0.5*τxx/τII 
    # dQdτyy =  0.5*τyy/τII
    # dQdτzz =  0.5*τzz/τII
    # dQdτxy =  τxy/τII/2
    # dQdP   = -sin(ψ)
    τII  = sqrt(0.5*(τ[1]^2 + τ[2]^2 + τ[3]^2) + τ[4]^2)
    f    = τII - τ[5]*sin(ϕ) - c*cos(ϕ) - λ̇*ηvp
    return f
end

@doc raw"""
    f = DruckerPrager_tot(τ, λ̇, ϕ, c, ηvp)  

Computes Drucker-Prager yield function for test 1 in Vermeer (1990) using total stresses. 
The function takes as input:

    τ   : vector of deviatoric stress components and pressure 
    λ̇   : plastic multiplier rate
    ϕ   : friction angle
    c   : cohesion
    ηvp : viscoplastic viscosity
and returns:

    f   : yield function value

# Examples
```julia-repl
julia>  DruckerPrager_tot(τ, λ̇, ϕ, c, ηvp)

```
"""
function DruckerPrager_tot(σ, λ̇, ϕ, c, ηvp)
    # τII    = sqrt(0.25*(τxx-τyy)^2 + τxy^2)
    # dQdτxx =  0.5*τxx/τII 
    # dQdτyy =  0.5*τyy/τII
    # dQdτzz =  0.5*τzz/τII
    # dQdτxy =  τxy/τII/2
    # dQdP   = -sin(ψ)
    P    = -1/3*(σ[1] + σ[2] + σ[3])
    τII  = sqrt(0.5*((σ[1]+P)^2 + (σ[2]+P)^2 + (σ[3]+P)^2) + σ[4]^2)
    f    = τII - P*sin(ϕ) - c*cos(ϕ) - λ̇*ηvp
    return f
end



