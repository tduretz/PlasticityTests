function Vermeer1990_Test1_Components(σi; params=(
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

    # Material properties
    G   = params.G
    ϕ   = params.ϕ
    ψ   = params.ψ
    
    # Strain integration
    γ̇xy = params.γ̇xy
    Δt  = params.Δt 
    nt  = params.nt 

    # Initial condition
    σxx = σi.xx 
    σyy = σi.yy
    σxy = 0. 
    σm  = [σxx σxy; σxy σyy]
    v   = eigvecs(σm)        
    θ   = atan(v[2,2] / v[1,2])   # σ3 angle
    
    # Storage
    app_fric    = zeros(nt)
    γxyv        = zeros(nt)
    εyyv        = zeros(nt)
    σxxv        = zeros(nt)
    θv          = zeros(nt)
    app_fric[1] = -σxy/σyy
    γxy,  εyy   = 0., 0. 

    for it=2:nt

        # σ3 angle
        σm   = [σxx σxy; σxy σyy]
        v    = eigvecs(σm)        
        θ    = atan(v[2,2] / v[1,2])   # σ3 angle
        τ    = sqrt(0.25*(σxx-σyy)^2 + σxy^2)
        # θ    =  0.5*asin(σxy/τ)
        θ    = 0.5*acos((σxx-σyy)/2/τ)

        # Stress integration
        τII  = sqrt(0.25*(σxx-σyy)^2 + σxy^2)
        f    = τII + 0.5*(σxx+σyy)*sin(ϕ) 
        if f>0
            Gs   = G*(2*sin(2θ)^2 + (cos(2θ) + sin(ψ))*(cos(2θ) + sin(ϕ)))^(-1)
            σ̇xx  = -2Gs*sin(2θ)*(cos(2θ) + sin(ψ))*γ̇xy
            σ̇yy  = 0.0
            σ̇xy  = Gs*(cos(2θ) + sin(ψ))*(cos(2θ) + sin(ϕ))*γ̇xy
        else
            σ̇xx  = 0.0
            σ̇yy  = 0.0
            σ̇xy  = G*γ̇xy
        end
        σxx += σ̇xx*Δt 
        σyy += σ̇yy*Δt
        σxy += σ̇xy*Δt

        γxy += γ̇xy*Δt
        # εyy += ε̇yy*Δt

        # Storage
        app_fric[it] = -σxy/σyy
        γxyv[it]     = γxy  
        # εyyv[it]     = εyy
        σxxv[it]     = σxx
        θv[it]       = θ
    end
    return (γxy=γxyv[2:end].*100, εyy=εyyv[2:end].*100, app_fric=app_fric[2:end], σxx=.-σxxv[2:end]./1e3, θ=θv[2:end].*(180/π))
end

function Vermeer1990_Test1_MatVec(σi; params=(
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

    # Material properties
    G   = params.G
    ϕ   = params.ϕ
    ψ   = params.ψ
    
    # Strain integration
    γ̇xy = params.γ̇xy
    Δt  = params.Δt 
    nt  = params.nt 

    𝐃ᵉ  = [2G 0 0; 0 2G 0; 0 0 G]
    
    # Initial condition
    σ   = [σi.xx; σi.yy; 0.]
    σm  = [σ[1] σ[3]; σ[3] σ[2]]
    v   = eigvecs(σm)        
    θ   = atan(v[2,2] / v[1,2])   # σ3 angle
 
    # Storage
    app_fric    = zeros(nt)
    γxyv        = zeros(nt)
    εyyv        = zeros(nt)
    σxxv        = zeros(nt)
    θv          = zeros(nt)
    app_fric[1] = -σ[3]/σ[2]
    γxy,  εyy   = 0., 0. 

    for it=2:nt

        # σ3 angle
        σm   = [σ[1] σ[3]; σ[3] σ[2]]
        v    = eigvecs(σm)        
        θ    = atan(v[2,2] / v[1,2])   # σ3 angle
        τ    = sqrt(0.25*(σ[1]-σ[2])^2 + σ[3]^2)
        θ    = 0.5*acos((σ[1]-σ[2])/2/τ)
        
        f      = τ + 0.5*(σ[1]+σ[2])*sin(ϕ) 

        # Stress integration
        𝐚    = G.*[ cos(2θ) + sin(ψ)
                   -cos(2θ) + sin(ψ) 
                    sin(2θ)         ] 
        𝐛    = G.*[ cos(2θ) + sin(ϕ)
                -cos(2θ) + sin(ϕ) 
                 sin(2θ)         ] 
        d    = G*(1 + sin(ψ)*sin(ϕ))
        𝐌    = 𝐃ᵉ - (f>=0)*1/d*𝐚*𝐛'
        ε̇    = [0.;   -(𝐌[2,3]/𝐌[2,2])*γ̇xy;   γ̇xy]
        σ  .+= 𝐌*ε̇*Δt

        τ    = sqrt(0.25*(σ[1]-σ[2])^2 + σ[3]^2)
        θ    = 0.5*acos((σ[1]-σ[2])/2/τ)
        
        γxy += γ̇xy*Δt
        # εyy += ε̇yy*Δt

        # Storage
        app_fric[it] = -σ[3]/σ[2]
        γxyv[it]     = γxy  
        # εyyv[it]     = εyy
        σxxv[it]     = σ[1]
        θv[it]       = θ
    end
    return (γxy=γxyv[2:end].*100, εyy=εyyv[2:end].*100, app_fric=app_fric[2:end], σxx=.-σxxv[2:end]./1e3, θ=θv[2:end].*(180/π))
end