@doc raw"""
    f = MohrCoulombVermeer1990_vdev(τ, λ̇, ϕ, c, ηvp, θt)  

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
julia>  MohrCoulombVermeer1990_vdev(τ, λ̇, ϕ, c, ηvp, θt)

```
"""
function MohrCoulombVermeer1990_vdev(τ, λ̇, ϕ, c, ηvp, θt)
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
    f = MohrCoulombVermeer1990_tot(σ, λ̇, ϕ, c, ηvp, θt)  

Computes Mohr-Coulomb yield function for test 1 in Vermeer (1990) using total stresses. 
The function takes as input:

    σ   : vector of total stress components
    λ̇   : plastic multiplier rate
    ϕ   : friction angle
    c   : cohesion
    ηvp : viscoplastic viscosity
    θt  : unused
and returns:

    f   : yield function value

# Examples
```julia-repl
julia>  MohrCoulombVermeer1990_tot(σ, λ̇, ϕ, c, ηvp, θt)

```
"""
function MohrCoulombVermeer1990_tot(σ, λ̇, ϕ, c, ηvp, θt)
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
    f = DruckerPrager_vdev(τ, λ̇, ϕ, c, ηvp, θt)  

Computes Drucker-Prager yield function for test 1 in Vermeer (1990) using deviatoric stresses and pressure. 
The function takes as input:

    τ   : vector of deviatoric stress components and pressure 
    λ̇   : plastic multiplier rate
    ϕ   : friction angle
    c   : cohesion
    ηvp : viscoplastic viscosity
    θt  : unused
and returns:

    f   : yield function value

# Examples
```julia-repl
julia>  DruckerPrager_vdev(τ, λ̇, ϕ, c, ηvp, θt)

```
"""
function DruckerPrager_vdev(τ, λ̇, ϕ, c, ηvp, θt)
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
    f = DruckerPrager_tot(τ, λ̇, ϕ, c, ηvp, θt)  

Computes Drucker-Prager yield function for test 1 in Vermeer (1990) using total stresses. 
The function takes as input:

    τ   : vector of deviatoric stress components and pressure 
    λ̇   : plastic multiplier rate
    ϕ   : friction angle
    c   : cohesion
    ηvp : viscoplastic viscosity
    θt  : unused
and returns:

    f   : yield function value

# Examples
```julia-repl
julia>  DruckerPrager_tot(τ, λ̇, ϕ, c, ηvp, θt)

```
"""
function DruckerPrager_tot(σ, λ̇, ϕ, c, ηvp, θt)
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

Lode(τII, J3) = -3.0*sqrt(3.0)/2.0*J3/τII^3

function MohrCoulomb_AS95_vdev(τ, λ̇, ϕ, c, ηvp, θt)
    τII = sqrt(0.5*(τ[1]^2 + τ[2]^2 + τ[3]^2) + τ[4]^2)
    J3  = τ[1]*τ[2]*τ[3] + τ[3]*τ[4]^2 # + 2*τ[4]*τ[5]*τ[6] + τ[1]*τ[6]^2 + τ[2]*τ[5]^2
    L   = Lode(τII,J3)
    L> 1.0 ? L= 1.0 : nothing
    L<-1.0 ? L=-1.0 : nothing
    θ   =  1.0/3.0*asin(-L) 
    if abs(θ)>θt
        sgnθ = sign(θ)
        A    = 1/3*cos(θt)*(3+tan(θt)*tan(3*θt) + 1/sqrt(3)*sgnθ*(tan(3*θt)-3*tan(θt))*sin(ϕ))
        B    = 1/(3*cos(3*θt))*(sgnθ*sin(θt) + 1/sqrt(3)*sin(ϕ)*cos(θt))
        k    = A - B*sin(3*θ)
    else
        k    = cos(θ) - 1/sqrt(3)*sin(ϕ)*sin(θ)
    end
    F   = k*τII - τ[5]*sin(ϕ) - c*cos(ϕ) - ηvp*λ̇
    return F
end

function MohrCoulomb_deBorst90_vdev(τ, λ̇, ϕ, c, ηvp, θt)
    P   = τ[5]
    J2  = 0.5*(τ[1]^2 + τ[2]^2 + τ[3]^2) + τ[4]^2
    J3  = τ[1]*τ[2]*τ[3] + τ[3]*τ[4]^2
    L   = -3/2*sqrt(3)*J3/J2/sqrt(J2)
    L> 1.0 ? L= 1.0 : nothing
    L<-1.0 ? L=-1.0 : nothing
    α   = 1/3*asin(L)
    σ3  = -(P + 2*sqrt(1/3*J2) * sin(α - 2/3*π))
    σ2  =   P + 2*sqrt(1/3*J2) * sin(α)
    σ1  = -(P + 2*sqrt(1/3*J2) * sin(α + 2/3*π))
    F   = 1/2*(σ3 - σ1) + 1/2*(σ3 + σ1)*sin(ϕ) - ηvp*λ̇
    return F
end