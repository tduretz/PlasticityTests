using Plots, LinearAlgebra

# Vermeer (1990)
# cohesion is 0!

function Vermeer1990_FirstTest_Components(σxxi, σyyi)

    # Material properties
    G   = 10e6
    ϕ   = 40/180*π
    ψ   = 10/180*π
    
    # Strain integration
    γ̇xy = 0.00001
    Δt  = 20
    nt  = 400

    # Initial condition
    σxx = σxxi 
    σyy = σyyi
    σxy = 0. 
    σm  = [σxx σxy; σxy σyy]
    v   = eigvecs(σm)        
    θ   = atan(v[2,2] / v[1,2])   # σ3 angle
    
    # Storage
    app_fric    = zeros(nt)
    γxy         = zeros(nt)
    θv          = zeros(nt)
    app_fric[1] = -σxy/σyy

    for it=2:nt

        # σ3 angle
        σm   = [σxx σxy; σxy σyy]
        v    = eigvecs(σm)        
        θ    = atan(v[2,2] / v[1,2])   # σ3 angle
        τ    = sqrt(0.25*(σxx-σyy)^2 + σxy^2)
        # θ    =  0.5*asin(σxy/τ)
        θ    = 0.5*acos((σxx-σyy)/2/τ)

        # Stress integration
        Gs   = G*(2*sin(2θ)^2 + (cos(2θ) + sin(ψ))*(cos(2θ) + sin(ϕ)))^(-1)
        σ̇xx  = -2Gs*sin(2θ)*(cos(2θ) + sin(ψ))*γ̇xy
        σ̇yy  = 0.0
        σ̇xy  = Gs*(cos(2θ) + sin(ψ))*(cos(2θ) + sin(ϕ))*γ̇xy
        σxx += σ̇xx*Δt 
        σyy += σ̇yy*Δt
        σxy += σ̇xy*Δt

        # Storage
        app_fric[it] = -σxy/σyy
        γxy[it]      = γxy[it-1] + γ̇xy*Δt
        θv[it]       = θ
    end
    return γxy, θv, app_fric
end


function Vermeer1990_FirstTest_Matrix(σxxi, σyyi)

    # Material properties
    G   = 10e6
    ϕ   = 40/180*π
    ψ   = 10/180*π
    𝐃ᵉ  = [2G 0 0; 0 2G 0; 0 0 G]
    
    # Strain integration
    γ̇xy = 0.00001
    Δt  = 20
    nt  = 400

    # Initial condition
    σ   = [σxxi; σyyi; 0.]
    σm  = [σ[1] σ[3]; σ[3] σ[2]]
    v   = eigvecs(σm)        
    θ   = atan(v[2,2] / v[1,2])   # σ3 angle
 
    # Storage
    app_fric    = zeros(nt)
    γxy         = zeros(nt)
    θv          = zeros(nt)
    app_fric[1] = -σ[3]/σ[2]

    for it=2:nt

        # σ3 angle
        σm   = [σ[1] σ[3]; σ[3] σ[2]]
        v    = eigvecs(σm)        
        θ    = atan(v[2,2] / v[1,2])   # σ3 angle
        τ    = sqrt(0.25*(σ[1]-σ[2])^2 + σ[3]^2)
        # θ    =  0.5*asin(σ[3]/τ)
        θ    = 0.5*acos((σ[1]-σ[2])/2/τ)

        # Stress integration
        𝐚    = G.*[ cos(2θ) + sin(ψ)
                   -cos(2θ) + sin(ψ) 
                    sin(2θ)         ] 
        𝐛    = G.*[ cos(2θ) + sin(ϕ)
                -cos(2θ) + sin(ϕ) 
                 sin(2θ)         ] 
        d    = G*(1 + sin(ψ)*sin(ϕ))
        𝐌    = 𝐃ᵉ - 1/d*𝐚*𝐛'
        ε̇    = [0.;   -(𝐌[2,3]/𝐌[2,2])*γ̇xy;   γ̇xy]
        σ  .+= 𝐌*ε̇*Δt

        # Storage
        app_fric[it] = -σ[3]/σ[2]
        γxy[it]      = γxy[it-1] + γ̇xy*Δt
        θv[it]       = θ
    end
    return γxy, θv, app_fric
end

function main()

    # Material properties
    G   = 10e6
    ϕ   = 40/180*π
    ψ   = 10/180*π
    #----------------------
    σh       = -100e3  
    σv       = (1+sin(ϕ))/(1-sin(ϕ))*σh
    #----------------------


    𝐃ᵉ  = [2G 0 0; 0 2G 0; 0 0 G]
    
    # Strain integration
    γ̇xy = 0.00001
    Δt  = 20
    nt  = 400

    # # Initial condition
    # σ   = [σxxi; σyyi; 0.]
    # σm  = [σ[1] σ[3]; σ[3] σ[2]]
    # v   = eigvecs(σm)        
    # θ   = atan(v[2,2] / v[1,2])   # σ3 angle
 
    # # Storage
    # app_fric    = zeros(nt)
    # γxy         = zeros(nt)
    # θv          = zeros(nt)
    # app_fric[1] = -σ[3]/σ[2]

    # for it=2:nt

    #     # σ3 angle
    #     σm   = [σ[1] σ[3]; σ[3] σ[2]]
    #     v    = eigvecs(σm)        
    #     θ    = atan(v[2,2] / v[1,2])   # σ3 angle
    #     τ    = sqrt(0.25*(σ[1]-σ[2])^2 + σ[3]^2)
    #     # θ    =  0.5*asin(σ[3]/τ)
    #     θ    = 0.5*acos((σ[1]-σ[2])/2/τ)

    #     # Stress integration
    #     𝐚    = G.*[ cos(2θ) + sin(ψ)
    #                -cos(2θ) + sin(ψ) 
    #                 sin(2θ)         ] 
    #     𝐛    = G.*[ cos(2θ) + sin(ϕ)
    #             -cos(2θ) + sin(ϕ) 
    #              sin(2θ)         ] 
    #     d    = G*(1 + sin(ψ)*sin(ϕ))
    #     𝐌    = 𝐃ᵉ - 1/d*𝐚*𝐛'
    #     ε̇    = [0.;   -(𝐌[2,3]/𝐌[2,2])*γ̇xy;   γ̇xy]
    #     σ  .+= 𝐌*ε̇*Δt

    #     # Storage
    #     app_fric[it] = -σ[3]/σ[2]
    #     γxy[it]      = γxy[it-1] + γ̇xy*Δt
    #     θv[it]       = θ
    # end

    # if style == :componentwise
    #     @info "Run with componentwise approach"
    #     γxy_A, θv_A, app_fric_A = Vermeer1990_FirstTest_Components(σxxA, σyyA)
    #     γxy_B, θv_B, app_fric_B = Vermeer1990_FirstTest_Components(σxxB, σyyB)
    # else
    #     @info "Run with matrix-vector approach"
    #     γxy_A, θv_A, app_fric_A = Vermeer1990_FirstTest_Matrix(σxxA, σyyA)
    #     γxy_B, θv_B, app_fric_B = Vermeer1990_FirstTest_Matrix(σxxB, σyyB)
    # end

    # p1 = plot( xlabel = "γxy",  ylabel = "-σxy/σyy" )
    # p1 = plot!(γxy_A[2:end], app_fric_A[2:end], label="case A")
    # p1 = plot!(γxy_B[2:end], app_fric_B[2:end], label="case B")
    # p2 = plot( xlabel = "γxy",  ylabel = "θ (σ₃ angle)" )
    # p2 = plot!(γxy_A[2:end], θv_A[2:end]*180/π, label="case A")
    # p2 = plot!(γxy_B[2:end], θv_B[2:end]*180/π, label="case B")

    # display(plot(p1, p2, layout=(2,1)))

end

main()

