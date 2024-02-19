using PlasticityTests, GeoParams, Plots, Printf, MathTeXEngine, LinearAlgebra, StaticArrays, Statistics
import LinearAlgebra:norm
# Makie.update_theme!(fonts = (regular = texfont(), bold = texfont(:bold), italic = texfont(:italic)))

function PrincipalStress!(σ1, σ3, τxx, τyy, τzz, τxy, P)
    for i in eachindex(τxy)
        σ  = @SMatrix[-P[i]+τxx[i] τxy[i] 0.; τxy[i] -P[i]+τyy[i] 0.; 0. 0. -P[i]+τzz[i]]
        v  = eigvecs(σ)
        σp = eigvals(σ)
        σ1.x[i] = v[1,1]
        σ1.z[i] = v[2,1]
        σ3.x[i] = v[1,3]
        σ3.z[i] = v[2,3]
        σ1.v[i] = σp[1]
        σ3.v[i] = σp[3]
    end
end

function Main_VEP_1D_standalone(σi; params=(
    #---------------#

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
    pl  = true), 
    #---------------#
    visu     = true, 
    make_gif = false,  
    Ncy      = 10# default parameter set
    #---------------#
    )

    # Visualisation is important!
    if visu==false 
        make_gif = false 
    end
    
    # Unit system
    CharDim    = SI_units(length=1000m, temperature=1000C, stress=1e7Pa, viscosity=1e20Pas)
    # CharDim    = SI_units(length=1m, temperature=1C, stress=1Pa, viscosity=1Pas)

    # External load
    ϕ          = params.ϕ
    ψ          = params.ψ  

    σh = nondimensionalize( (-100e3)Pa, CharDim) 
    σv   = (1 + sin(ϕ))/(1 - sin(ϕ))*σh
    θ_A  = π/4 + 0.25*(ϕ + ψ)
    θ_C  = π/4 + 0.5*(ϕ)
    θ_R  = π/4 + 0.5*(ψ)
    θ_SB = θ_A

    # to Cartesian
    σxxi = 1/2*(σh + σv) +  1/2*(σh - σv)*cos(2*θ_SB)
    σyyi = 1/2*(σh + σv) -  1/2*(σh - σv)*cos(2*θ_SB)
    σxyi = 1/2*(σh - σv)*sin(2*θ_SB)

    if σi.xx>σi.yy
        Vemeer90 = ExtractDataCase("CaseA")
    else
        Vemeer90 = ExtractDataCase("CaseB")
    end

    # Physical parameters
    σzzi       = 0.5*(σxxi + σyyi)
    Pi         = -(σxxi + σyyi)/2.0
    τxxi       = Pi + σxxi
    τyyi       = Pi + σyyi
    τzzi       = Pi + σzzi
    τxyi       = σxyi

    Ly         = nondimensionalize(2e4m, CharDim)
    Ẇ0         = nondimensionalize(5e-5Pa/s, CharDim)
    ε0         = nondimensionalize((params.γ̇xy)s^-1, CharDim)
    G          = nondimensionalize((params.G)Pa, CharDim)
    Coh0       = nondimensionalize((params.c)Pa, CharDim)
    μs         = nondimensionalize(1e52Pa*s, CharDim) 
    ηvp        = nondimensionalize((params.ηvp)Pa*s, CharDim)     # 7.5e7

    # Numerical parameters
    Nt         = params.nt-1
    Δy         = Ly/Ncy
    yc         = LinRange(-Ly/2-Δy/2, Ly/2+Δy/2, Ncy+2)
    yv         = LinRange(-Ly/2,      Ly/2,      Ncy+1)
    Δt         = nondimensionalize((params.Δt)s, CharDim) /10
    ηe         = G*Δt
    nout_viz   = 10
    if make_gif 
        nout_viz = 1 
    end

    # Allocate arrays
    Kb         = 2/3*G*ones(Ncy+1);# Kb[1] = Kb[1]*2
    Pt         =  Pi*ones(Ncy+1) 
    Ptc        =  Pi*ones(Ncy+1) 
    Pt0        =  Pi*ones(Ncy+1)
    τxx        =  τxxi*ones(Ncy+1); 
    τxy        =  τxyi*ones(Ncy+1)
    τyy        =  τyyi*ones(Ncy+1);  
    τzz        =  τzzi*ones(Ncy+1)  
    τxx0       =  τxxi*ones(Ncy+1)
    τxy0       =  τxyi*ones(Ncy+1)
    τyy0       =  τyyi*ones(Ncy+1)
    τzz0       =  τzzi*ones(Ncy+1)
    Coh        =  zeros((Ncy+1));    
    F          =  zeros((Ncy+1))
    Fc         =  zeros((Ncy+1))
    λ̇          =  zeros((Ncy+1))
    λ̇rel       =  zeros((Ncy+1))
    ispl       =  zeros(Int, (Ncy+1))
    ηve        =  zeros((Ncy+1)); ηve        .= 1.0./(1.0/μs .+ 1.0/ηe)
    ηvep       =  zeros((Ncy+1))
    εxy        =    zeros(Ncy+1)
    εyy        =    zeros(Ncy+1)
    ε̇xy        =  ε0*ones(Ncy+1)
    ε̇xx        =    zeros(Ncy+1)
    ε̇yy        =    zeros(Ncy+1)
    ε̇zz        =    zeros(Ncy+1)
    ε̇xxd       =    zeros(Ncy+1)
    ε̇yyd       =    zeros(Ncy+1)
    ε̇zzd       =    zeros(Ncy+1)
    ε̇iiᵉᶠᶠ     =    zeros(Ncy+1)
    τii        =    zeros(Ncy+1)
    η          =    zeros(Ncy+1)
    ΔτV        =    zeros(Ncy)
    ΔτPt       =    zeros(Ncy+1)
    η_mm       =    zeros(Ncy)
    Vx         =    zeros(Ncy+2);  Vx .= ε0.*yc
    Vy         =    zeros(Ncy+2);
    RPt        =    zeros(Ncy+1)
    RVx        =    zeros(Ncy+2)
    RVy        =    zeros(Ncy+2)
    ∂Pt∂τ      =    zeros(Ncy+1)
    ∂Vx∂τ      =    zeros(Ncy+2)
    ∂Vy∂τ      =    zeros(Ncy+2)
    σ1         = (x=zeros(size(τxx)), z=zeros(size(τxx)), v=zeros(size(τxx)) )
    σ3         = (x=zeros(size(τxx)), z=zeros(size(τxx)), v=zeros(size(τxx)) )
 
    ε̇xy_pl     =   zeros(Ncy+1)
    ε̇xx_pl     =   zeros(Ncy+1)
    ε̇yy_pl     =   zeros(Ncy+1)
    ε̇zz_pl     =   zeros(Ncy+1)
    ∇v         =   zeros(Ncy+1)
    ∇v_pl      =   zeros(Ncy+1)

    # Monitoring
    probes    = (Ẇ0 = zeros(Nt), τxy0 = zeros(Nt), σyy0 = zeros(Nt), Vx0 = zeros(Nt), τii, θs3 = zeros(Nt), θs3_out = zeros(Nt), θs3_in = zeros(Nt), fric_in = zeros(Nt), fric_out = zeros(Nt), εyy = zeros(Nt), σxx=zeros(Nt), fric=zeros(Nt), σv_σh=zeros(Nt), γxySB=zeros(Nt))
    η        .= μs
   
    # BC
    BC_Vy = :Neumann
    # BC_Vy = :Dirichlet
    VxS   =  ε0*yv[1]
    VxN   =  ε0*yv[end]
    VyS   =  0.0
    VyN   =  0.0

    # PT solver
    niter = 10000
    θVx   = 0.6
    θVy   = 0.6
    θPt   = 1.0
    nout  = 1000
    ϵ     = 1e-10
    rel   = 1e-2
    errPt, errVx, errVy = 0., 0., 0.

    PrincipalStress!(σ1, σ3, τxx, τyy, τzz, τxy, Pt)
    println(ustrip.(dimensionalize( σ1[1], Pa, CharDim)./1e3))
    println(ustrip.(dimensionalize( σ3[1], Pa, CharDim)./1e3))

    # fsxx = σyyi - 1/2*(σh + σv) -  1/2*(σh - σv)*cos(2*θ_SB)
    # σxxi = (τxx[end]-Pt[end])*cos(θ_SB)^2 + (τyy[end]-Pt[end])*sin(θ_SB)^2 + (τxy[end])*sin(2*θ_SB)
    # σyyi = (τxx[end]-Pt[end])*sin(θ_SB)^2 + (τyy[end]-Pt[end])*cos(θ_SB)^2 - (τxy[end])*sin(2*θ_SB)

    # println(ustrip.(dimensionalize( σxxi, Pa, CharDim)./1e3))
    # println(ustrip.(dimensionalize( σyyi, Pa, CharDim)./1e3))
    # println(ustrip.(dimensionalize( σv, Pa, CharDim)./1e3))
    # println(ustrip.(dimensionalize( σh, Pa, CharDim)./1e3))

    σrel = 0.99999

    dt_red = 1#e4

    anim = @animate for it=1:Nt
        # History
        @. τxy0  = τxy
        @. τxx0  = τxx
        @. τyy0  = τyy
        @. τzz0  = τzz
        @. Pt0   = Pt
        @. λ̇     = 0.0
        @. λ̇rel  = 0.0
        
        τyy_blah = τyy[end]
        τxx_blah = τxx[end]
        τxy_blah = τxy[end]
        Pt_blah  = Pt[end]
        σxx_trial, σyy_trial = 0., 0.

        @views for iter=1:niter

            σxx_trial = (τxx[end] - Pt[end])*cos(θ_SB)^2 + (τyy[end] - Pt[end])*sin(θ_SB)^2 + (τxy[end])*sin(2*θ_SB)
            σyy_trial = (τxx[end] - Pt[end])*sin(θ_SB)^2 + (τyy[end] - Pt[end])*cos(θ_SB)^2 - (τxy[end])*sin(2*θ_SB)

            scale     = abs(σh/σxx_trial)

            τyy_blah = τyy_blah*σrel + (1-σrel)*τyy[end]*scale 
            τxx_blah = τxx_blah*σrel + (1-σrel)*τxx[end]*scale
            τxy_blah = τxy_blah*σrel + (1-σrel)*τxy[end]*scale
            Pt_blah  = Pt_blah *σrel + (1-σrel)* Pt[end]*scale

            # Kinematics
            # Vx BC
            Vx[1]   = -Vx[2]     + 2VxS
            Vx[end] = -Vx[end-1] + 2VxN

            # Vx[1]   = -Vx[2]
            # Vx[end] = Vx[end-1] + 2 * Δy .* ε̇xy_pl[end] + Δy .* σxyi ./ ηve[end] - Δy .* τxy0[end] ./ ηe[end]            

            # Vy BC
            Vy[1]   = - Vy[2]     + 2VyS
            if BC_Vy == :Dirichlet
                Vy[end] = - Vy[end-1] + 2VyN
            elseif BC_Vy == :Neumann  
                # σyyi    such that σh is -100 kPa
                σyyi =  τyy_blah - Pt_blah
                # @printf("σyy_desire = %2.4e\n", ustrip(dimensionalize(τyy_blah - Pt_blah, Pa, CharDim)/1000))
                # Vy[end] = Vy[end-1]
                # Vy[end] = (Pt0[end] .* Δy .* ηe[end] + Vy[end-1] .* ηe[end] .^ 2 + Vy[end-1] .* ηe[end] .* ηve[end] + Δy .* ηe[end] .* σyyi - Δy .* ηve[end] .* τyy0[end]) ./ (ηe[end] .* (ηe[end] + ηve[end]))
                # Vy[end]  = (Pt0[end] .* Δy .* ηe[end] + Vy[end-1] .* ηe[end] .^ 2 + Vy[end-1] .* ηe[end] .* ηve[end] - 0.666666666666667 * Δy .* ε̇xx_pl[end] .* ηe[end] .* ηve[end] + 1.33333333333333 * Δy .* ε̇yy_pl[end] .* ηe[end] .* ηve[end] - 0.666666666666667 * Δy .* ε̇zz_pl[end] .* ηe[end] .* ηve[end] + Δy .* ηe[end] .* σyyi - Δy .* ηve[end] .* τyy0[end]) ./ (ηe[end] .* (ηe[end] + ηve[end]))
                Vy[end] = (3.0 * Kb[end] .* Vy[end-1] .* Δt .* ηe[end] + 3.0 * Kb[end] .* Δt .* Δy .* ηe[end] .* ∇v_pl[end] + 3.0 * Pt0[end] .* Δy .* ηe[end] + 4.0 * Vy[end-1] .* ηe[end] .* ηve[end] + 6.0 * Δy .* ε̇yy_pl[end] .* ηe[end] .* ηve[end] + 3.0 * Δy .* ηe[end] .* σyyi - 3.0 * Δy .* ηve[end] .* τyy0[end]) ./ (ηe[end] .* (3.0 * Kb[end] .* Δt + 4.0 * ηve[end]))
            end
           
            @. ε̇xy  =  0.5*(Vx[2:end] - Vx[1:end-1])/Δy
            @. ε̇yy  =      (Vy[2:end] - Vy[1:end-1])/Δy # total
            @. ε̇zz  = 1/2*(ε̇xx + ε̇yy)
            @. ∇v   = ε̇xx + ε̇yy + ε̇zz

            # Deviatoric strain rates
            @. ε̇xxd = ε̇xx - 1/3*∇v
            @. ε̇yyd = ε̇yy - 1/3*∇v
            @. ε̇zzd = ε̇zz - 1/3*∇v

            # Stress
            @. τxx     =  4/3 * ηve * (ε̇xxd + 1/3*∇v + τxx0/(4/3*ηe)) - 2/3 * ηve * (ε̇yyd + 1/3*∇v)                 - 2/3 * ηve * (ε̇zzd + 1/3*∇v)
            @. τyy     = -2/3 * ηve * (ε̇xxd + 1/3*∇v)                 + 4/3 * ηve * (ε̇yyd + 1/3*∇v + τyy0/(4/3*ηe)) - 2/3 * ηve * (ε̇zzd + 1/3*∇v)
            @. τzz     = -2/3 * ηve * (ε̇xxd + 1/3*∇v)                 - 2/3 * ηve * (ε̇yyd + 1/3*∇v)                 + 4/3 * ηve * (ε̇zzd + 1/3*∇v + τzz0/(4/3*ηe))
            @. τxy     =    2 * ηve * (ε̇xy  + τxy0/(2*ηe)) 
            @. τii     = sqrt(τxy^2 + 0.5*(τyy.^2 + τxx.^2 + τzz.^2))
            # @. Pt     = Pt0 - 2/3*G*Δt*∇v

            # Plasticity
            @. F    = τii - Coh*cos(ϕ) - Pt*sin(ϕ)
            @. Ptc  = Pt
            @. ηvep = ηve
            @. ispl = 0
            @. ispl[F>=0] = 1
            @. ε̇iiᵉᶠᶠ   = sqrt( (ε̇xy + τxy0/2/ηe)^2 + 0.5*( (ε̇xxd + τxx0/(4/3*ηe))^2 + ((ε̇yyd + τyy0/(4/3*ηe))).^2 + ((ε̇zzd + τzz0/(4/3*ηe))).^2 ) ) 
            for it=1:50
                # @. λ̇      = F / (ηvp + ηve + Kb*Δt*sin(ϕ)*sin(ψ))
                # @. λ̇rel   = (1.0-rel)*λ̇rel + rel*λ̇   
                # @. λ̇rel   = λ̇   
                @. ηvep   = (Coh*cos(ϕ) + Ptc*sin(ϕ) + ηvp*λ̇rel) / 2.0 / ε̇iiᵉᶠᶠ
                @. ε̇xx_pl = λ̇rel*(τxx/2/τii)
                @. ε̇yy_pl = λ̇rel*(τyy/2/τii)
                @. ε̇zz_pl = λ̇rel*(τxx/2/τii + τyy/2/τii)/2   # dqdτzz*λ̇rel
                @. ε̇xy_pl = λ̇rel*(τxy/2/τii)
                @. ∇v_pl  = 3/2*sin(ψ)*λ̇rel
                @. Ptc    = Pt0 - Kb*Δt*(∇v - ∇v_pl)
                @. τxx    =  4/3 * ηve * (ε̇xxd + τxx0/(4/3*ηe) -  ε̇xx_pl) - 2/3 * ηve * (ε̇yyd                 -  ε̇yy_pl) - 2/3 * ηve * (ε̇zzd                 -  ε̇zz_pl)
                @. τyy    = -2/3 * ηve * (ε̇xxd                 -  ε̇xx_pl) + 4/3 * ηve * (ε̇yyd + τyy0/(4/3*ηe) -  ε̇yy_pl) - 2/3 * ηve * (ε̇zzd                 -  ε̇zz_pl)
                @. τzz    = -2/3 * ηve * (ε̇xxd                 -  ε̇xx_pl) - 2/3 * ηve * (ε̇yyd                 -  ε̇yy_pl) + 4/3 * ηve * (ε̇zzd + τzz0/(4/3*ηe) -  ε̇zz_pl)
                @. τxy    =    2 * ηve * (ε̇xy  + τxy0/(2*ηe)   -  ε̇xy_pl) 
                @. τii    = sqrt(τxy^2 + 0.5*(τyy^2 + τxx^2 + τzz^2))
                @. Fc     = τii - Coh*cos(ϕ) - Ptc*sin(ϕ) - ηvp*λ̇rel
                @. λ̇rel  += (F>=0) .* Fc / (ηvp + ηve + Kb*Δt*sin(ϕ)*sin(ψ))
                @. λ̇rel[F<0] = 0.0
                if maximum(Fc) < ϵ break end 
            end

            # PT time steps
            @. η_mm = min.(ηvep[1:end-1], ηvep[2:end]); 
            @. ΔτV   = Δy^2/(η_mm)/2.1 /4 /dt_red
            @. ΔτPt  = ηvep/Δy/G/Δt/3 / dt_red
            
            # Residuals
            @. RPt          =  (- Kb*Δt*∇v - (Pt - Pt0))
            @. RVx[2:end-1] =  ((τxy[2:end] - τxy[1:end-1])/Δy )
            @. RVy[2:end-1] =  ((τyy[2:end] - τyy[1:end-1])/Δy - (Ptc[2:end] - Ptc[1:end-1])/Δy)
            
            # Damp residuals
            @. ∂Vx∂τ = RVx + (1.0 - θVx)*∂Vx∂τ
            @. ∂Vy∂τ = RVy + (1.0 - θVy)*∂Vy∂τ
            @. ∂Pt∂τ = RPt + (1.0 - θPt)*∂Pt∂τ

            # Update solutions
            @. Vx[2:end-1] += ΔτV * ∂Vx∂τ[2:end-1] 
            @. Vy[2:end-1] += ΔτV * ∂Vy∂τ[2:end-1] 
            @. Pt += ΔτPt * ∂Pt∂τ 

            if mod(iter, nout) == 0 || iter==1
                errPt = norm(RPt)/sqrt(length(RPt))
                errVx = norm(RVx)/sqrt(length(RVx))
                errVy = norm(RVy)/sqrt(length(RVy))
                σyyBC = τyy[end] - Ptc[end]
                @printf("Iteration %05d --- Time step %4d --- Δt = %2.2e --- ΔtC = %2.2e --- εxy = %2.2e --- σyyBC = %2.7e --- max(F) = %2.2e --- max(Fc) = %2.2e \n", iter, it, ustrip(dimensionalize(Δt, s, CharDim)), ustrip(dimensionalize(Δy/2/maximum(Vx), s, CharDim)), ε0*it*Δt, ustrip(dimensionalize(σyyBC, Pa, CharDim)/1e3), maximum(ustrip.(dimensionalize(F, Pa, CharDim))), maximum(ustrip.(dimensionalize(Fc, Pa, CharDim))) )
                # @printf("Exy_net = %2.2e --- Eyy_net = %2.2e --- Div net = %2.2e\n", mean(abs.(ε̇xy_net)), mean(abs.(ε̇yy_net)), mean(abs.(∇v_net)) )
                # @printf("Exy_el  = %2.2e --- Exy_pl  = %2.2e --- Exy net = %2.2e\n", mean(abs.(ε̇xy_el)), mean(abs.(ε̇xy_pl)), mean(abs.(ε̇xy_net)) )
                @printf("fPt = %2.4e\n", errPt)
                @printf("fVx = %2.4e\n", errVx)
                @printf("fVy = %2.4e\n", errVy)
                σxx_trial = (τxx[end]-Pt[end])*cos(θ_SB)^2 + (τyy[end]-Pt[end])*sin(θ_SB)^2 + (τxy[end])*sin(2*θ_SB)
                @printf("σxx_trial = %2.4e\n", ustrip(dimensionalize(σxx_trial, Pa, CharDim)/1000))
                (errVx < ϵ && errVy < ϵ) && break 
                (isnan(errPt) || isnan(errVx) || isnan(errVx)) && error("NaNs")        
            end
        end        

        @. Pt   = Ptc
        @. εxy += ε̇xy*Δt 
        @. εyy += ε̇yy*Δt 

        # Show array infos
        # @minmax(Pt)
        # @minmax(τxx)
        # @minmax(τyy)
        # @minmax(τzz)
        # @minmax(τxy)

        # if (errPt > ϵ || errVx > ϵ || errVx > ϵ) error("non converged") end
        PrincipalStress!(σ1, σ3, τxx, τyy, τzz, τxy, Pt)

        _, iA = findmin(σ3.v) # out
        _, iB = findmax(σ1.v) # in
        # iA = Ncy+1


        probes.Ẇ0[it]       = τxy[end]*ε̇xy[end]
        probes.τxy0[it]     = τxy[end]
        probes.Vx0[it]      = 0.5*(Vx[end] + Vx[end-1])
        probes.σyy0[it]     = τyy[end] - Pt[end]
        probes.θs3_out[it]  = atand.(σ3.z[iA] ./ σ3.x[iA])
        probes.θs3_in[it]   = atand.(σ3.z[iB] ./ σ3.x[iB])
        probes.fric_out[it] = -τxy[iA]./(τyy[iA] .- Pt[iA])
        probes.fric_in[it]  = -τxy[iB]./(τyy[iB] .- Pt[iB])
        probes.εyy[it]      = εyy[iB]
        probes.σxx[it]      = τxx[iB] - Pt[iB]
        probes.fric[it]     = -τxy[iB]./(τyy[iB] .- Pt[iB])
        probes.θs3[it]      = atand.(σ3.z[iB] ./ σ3.x[iB])
        probes.σv_σh[it]    = σyy_trial/σxx_trial
        probes.γxySB[it]    = maximum(εxy)

        # Visualisation
        # if visu==true && (mod(it, nout_viz)==0 || it==1 || it==Nt)
    
            θ   = LinRange(-π, 0, 100)
            σMC = LinRange(-500, 0, 100 ) .*1e3
            τMC = -σMC.*tan(ϕ) 

            PA    = dimensionalize( (σ1.v[iA] + σ3.v[iA])/2, Pa, CharDim)
            τA    = dimensionalize( (σ1.v[iA] - σ3.v[iA])/2, Pa, CharDim)
            τB    = dimensionalize( (σ1.v[iB] - σ3.v[iB])/2, Pa, CharDim)
            PB    = dimensionalize( (σ1.v[iB] + σ3.v[iB])/2, Pa, CharDim)
            
            yield = (x = σMC./1e3, y = τMC./1e3)
            MC_A  = (x = (τA.*cos.(θ) .+ PA)./1e3, y = (τA.*sin.(θ))./1e3) 
            MC_B  = (x = (τB.*cos.(θ) .+ PB)./1e3, y = (τB.*sin.(θ))./1e3)
    
            # p1=plot( title = "Total pressure", xlabel = L"$P$ [kPa]", ylabel = L"$y$ [km]" )
            # p1=plot!(ustrip.(dimensionalize(τyy.-Pt, Pa, CharDim))/1e3, ustrip.(dimensionalize(yv, m, CharDim)./1e3), label="σyy" )
            # p1=plot!(ustrip.(dimensionalize(τyy[ispl.==1].-Pt[ispl.==1], Pa, CharDim))/1e3, ustrip.(dimensionalize(yv[ispl.==1], m, CharDim)./1e3), linewidth=5 )

            # p1 = plot( title = "γxy [%]", xlabel = "γxy [%]", ylabel = "y [-]", legend=:none )
            # p1 = plot!(ustrip.(2 .* εxy*100),           ustrip.(dimensionalize(yv, m, CharDim)./1e3), label="γxy" )
            p1 = plot(atand.(σ3.z./σ3.x), ustrip.(dimensionalize(yv, m, CharDim)./1e3), linewidth=5 )

            # p1 = plot( xlabel = "γxy shear band [%]", ylabel = "σxx/σxx0", foreground_color_legend = nothing, background_color_legend = nothing )
            # p1 = plot!(probes.γxySB[1:it]*100, probes.σxx[1:it]./σxxi, label="σxx/σxx0", color=:blue,  xlims=(0,8) )

            # p2=plot(title = "Velocity", xlabel = L"$Vx$ [cm/y]", ylabel = L"$y$ [-]" )
            # p2=plot!(ustrip.(dimensionalize(Vx, m/s, CharDim)), ustrip.(dimensionalize(yc, m, CharDim)./1e3) )
            p2 = plot(title="Mohr circles", ylabel="τ [kPa]", xlabel="σₙ [kPa]", size=(300,300), aspect_ratio=1, xlim=(-500,0), ylim=(0,400))
            p2 = plot!( MC_A... , color=:blue, label="out" )
            p2 = plot!( MC_B...,  color=:green, label="in"  )
            p2 = plot!( yield..., color=:red, label="Yield"  )
            p2 = plot!(ylabel="τ [kPa]", xlabel="σₙ [kPa]", foreground_color_legend = nothing, background_color_legend = nothing, legend=:topright)

            p3 = plot( ylabel = "θ σ3 [ᵒ]", xlabel = "γxy shear band [%]", xlims=(0,8), foreground_color_legend = nothing, background_color_legend = nothing  )
            p3 = plot!(probes.γxySB[1:it]*100, ustrip.(probes.θs3_out[1:it]), label="out", color=:blue  )
            p3 = plot!(probes.γxySB[1:it]*100, ustrip.(probes.θs3_in[1:it]),  label="in" , color=:green )

            p4 = plot( xlabel = "γxy shear band [%]", ylabel = "σv/σh", foreground_color_legend = nothing, background_color_legend = nothing )
            p4 = plot!(probes.γxySB[1:it]*100, probes.σv_σh[1:it], label="σv/σh", color=:blue,  xlims=(0,8) )
            # p4 = plot!((1:it)*ε0*Δt*100, probes.fric_out[1:it], label="out", color=:blue,  xlims=(0,8) )
            # p4 = plot!((1:it)*ε0*Δt*100, probes.fric_in[1:it],  label="in" , color=:green, xlims=(0,8) )
            # p4 = plot!(title=@sprintf("-σxy/σyy = %1.4f", probes.fric_in[it]))
            display(plot(p1,p2,p3,p4))
            
        # end
    end

    if make_gif gif(anim, "figures/Test1_MohrCircles_1D.gif", fps = 15) end
    return (γxy=(1:Nt)*ε0*Δt*100, εyy=probes.εyy.*100, app_fric=probes.fric, σxx=.-probes.σxx/1e3, θ=probes.θs3)
end

Main_VEP_1D_standalone((xx = -400e3, yy=-100e3));