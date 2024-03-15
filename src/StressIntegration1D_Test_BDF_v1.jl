
using Plots, Printf, LinearAlgebra

# This works with bdf2 and 3

function bdf2(dt0, dt)
    a = 1 ./ (dt + dt0) + 1 ./ dt
    b = dt .* (-1 ./ dt0 - 1 ./ dt) ./ (dt + dt0) - 1 ./ dt
    c = dt ./ (dt0 .* (dt + dt0))
    return (a, b, c)
end

function bdf3(dt00, dt0, dt)
    a = 1 ./ (dt + dt0 + dt00) + 1 ./ (dt + dt0) + 1 ./ dt
    b = dt .* (dt + dt0) .* ((-1 ./ dt0 - 1 ./ dt) ./ (dt + dt0) - 1 ./ (dt0 .* (dt0 + dt00))) ./ (dt + dt0 + dt00) + dt .* (-1 ./ dt0 - 1 ./ dt) ./ (dt + dt0) - 1 ./ dt
    c = dt .* (dt + dt0) .* (-(-1 ./ dt00 - 1 ./ dt0) ./ (dt0 + dt00) + 1 ./ (dt0 .* (dt + dt0))) ./ (dt + dt0 + dt00) + dt ./ (dt0 .* (dt + dt0))
    d = -dt .* (dt + dt0) ./ (dt00 .* (dt0 + dt00) .* (dt + dt0 + dt00))
    return (a, b, c, d)
end

function Main_VEP_1D_test_BDF(σi; params=(
    #---------------#
    K    = 6.6666666667e6, # K = 3/2*Gv in Vermeer (1990)
    G    = 10e6,
    c    = 0.,
    ϕ    = 40/180*π,
    ψ    = 10/180*π,
    θt   = 25/180*π,
    ηvp  = 0.,
    lc   = 5e2,
    γ̇xy  = 0.00001,
    Δt   = 20/4,
    nt   = 50*4,
    law  = :MC_Vermeer1990,
    coss = true,
    oop  = :Vermeer1990,
    pl   = true,
    noisy = false,
    bdf   = 3 ), 
    #---------------#
    visu     = true, 
    make_gif = false,  
    Ncy      = 40   # default parameter set
    #---------------#
    )
    sc = (σ = params.G, L = 1.0, t = 1.0/params.γ̇xy)

    coss  = params.coss
    noisy = params.noisy

    # Visualisation is important!
    if visu==false 
        make_gif = false 
    end
    
    # Load digitised data from test 1 of Vermeer (1990)
    σxxi     = σi.xx/sc.σ 
    σyyi     = σi.yy/sc.σ 
 
    # Physical parameters
    σzzi       = 0.5*(σxxi + σyyi)
    Pi         = -(σxxi + σyyi)/2.0
    τxxi       = Pi + σxxi
    τyyi       = Pi + σyyi
    τzzi       = Pi + σzzi
    τxyi       = 0.0*2e5/sc.σ

    # Dimensionally dependent
    Ly         = 2e4/sc.L
    Ẇ0         = 5e-5/(sc.σ/sc.t)
    ε0         = params.γ̇xy/(1.0/sc.t)
    G          = params.G/sc.σ
    Kb         = params.K/sc.σ
    Coh0       = params.c/sc.σ
    Coh1       =  1.0/100/sc.σ
    μs         = 1e10/(sc.σ*sc.t)
    ηvp        = params.ηvp/(sc.σ*sc.t)
    ϕ          = params.ϕ
    ψ          = params.ψ   
    Gc         = G
    lc         = params.lc/sc.L

    # Numerical parameters
    Nt         = params.nt-1
    Δy         = Ly/Ncy
    yc         = LinRange(-Ly/2-Δy/2, Ly/2+Δy/2, Ncy+2)
    yv         = LinRange(-Ly/2,      Ly/2,      Ncy+1)
    Δt         = params.Δt/sc.t
    a, b, c, d, e = 1/Δt, -1/Δt, 0/Δt, 0/Δt, 0
    
    nout_viz   = params.nt
    if make_gif 
        nout_viz = 1 
    end

    # Allocate arrays
    ηe         = G/a*ones(Ncy+1); #ηe[Int64(ceil(Ncy/2))] = G*Δt/5
    Pt         =  Pi*ones(Ncy+1) 
    Ptc        =  Pi*ones(Ncy+1) 
    Pt0        =  Pi*ones(Ncy+1)
    Pt00       =  Pi*ones(Ncy+1)
    Pt000      =  Pi*ones(Ncy+1)
    Pt0000     =  Pi*ones(Ncy+1)
    τxx        =  τxxi*ones(Ncy+1)
    τyy        =  τyyi*ones(Ncy+1)    
    τzz        =  τzzi*ones(Ncy+1) 
    τxy        =  τxyi*ones(Ncy+1)
    τxx0       =  zeros(Ncy+1)
    τyy0       =  zeros(Ncy+1)
    τzz0       =  zeros(Ncy+1)
    τxy0       =  zeros(Ncy+1)
    τxx00      =  zeros(Ncy+1)
    τyy00      =  zeros(Ncy+1)
    τzz00      =  zeros(Ncy+1)
    τxy00      =  zeros(Ncy+1)
    τxx000     =  zeros(Ncy+1)
    τyy000     =  zeros(Ncy+1)
    τzz000     =  zeros(Ncy+1)
    τxy000     =  zeros(Ncy+1)
    τxx0000    =  zeros(Ncy+1)
    τyy0000    =  zeros(Ncy+1)
    τzz0000    =  zeros(Ncy+1)
    τxy0000    =  zeros(Ncy+1)
    Coh        =  ones((Ncy+1)).*Coh1; 
    Coh[Int64(ceil(Ncy/2))] = Coh0  
    F          =  zeros((Ncy+1))
    Fc         =  zeros((Ncy+1))
    λ̇          =  zeros((Ncy+1))
    λ̇rel       =  zeros((Ncy+1))
    ispl       =  zeros(Int, (Ncy+1))
    ηve        =  zeros((Ncy+1)); 
    ηvep       =  zeros((Ncy+1))
    εxy        =    zeros(Ncy+1)
    Ẇz         =    zeros(Ncy+1)
    κ̇yz        =    zeros(Ncy+1)
    εyy        =    zeros(Ncy+1)
    ε̇xy        =  ε0*ones(Ncy+1)
    ε̇xx        =    zeros(Ncy+1)
    ε̇yy        =    zeros(Ncy+1)
    ε̇zz        =    zeros(Ncy+1)
    ε̇xxd       =    zeros(Ncy+1)
    ε̇yyd       =    zeros(Ncy+1)
    ε̇zzd       =    zeros(Ncy+1)
    myz        =    zeros(Ncy+1)
    Rz         =    zeros(Ncy+1)
    Rzc        =    zeros(Ncy+2)
    myz0       =    zeros(Ncy+1)
    Rz0        =    zeros(Ncy+1)
    myz00      =    zeros(Ncy+1)
    Rz00       =    zeros(Ncy+1)
    myz000     =    zeros(Ncy+1)
    Rz000      =    zeros(Ncy+1)
    myz0000    =    zeros(Ncy+1)
    Rz0000     =    zeros(Ncy+1)
    ε̇iiᵉᶠᶠ     =    zeros(Ncy+1)
    τii        =    zeros(Ncy+1)
    η          =    zeros(Ncy+1)
    ΔτV        =    zeros(Ncy)
    ΔτPt       =    zeros(Ncy+1)
    Δτω̇z       =    zeros(Ncy) 
    η_mm       =    zeros(Ncy)
    Vx         =    zeros(Ncy+2);  Vx .= ε0.*yc #.+ rand(size(Vx))
    Vy         =    zeros(Ncy+2)
    ω̇z         =    zeros(Ncy+2)
    ω̇zv        =    zeros(Ncy+1)
    RPt        =    zeros(Ncy+1)
    RVx        =    zeros(Ncy+2)
    RVy        =    zeros(Ncy+2)
    Rω̇z        =    zeros(Ncy+2)
    ∂Pt∂τ      =    zeros(Ncy+1)
    ∂Vx∂τ      =    zeros(Ncy+2)
    ∂Vy∂τ      =    zeros(Ncy+2)
    ∂ω̇z∂τ      =    zeros(Ncy+2)
    σ1         = (x=zeros(size(τxx)), z=zeros(size(τxx)), v=zeros(size(τxx)) )
    σ3         = (x=zeros(size(τxx)), z=zeros(size(τxx)), v=zeros(size(τxx)) )
 
    κ̇yz_pl     =   zeros(Ncy+1)
    ẇz_pl      =   zeros(Ncy+1)
    ε̇xy_pl     =   zeros(Ncy+1)
    ε̇xx_pl     =   zeros(Ncy+1)
    ε̇yy_pl     =   zeros(Ncy+1)
    ε̇zz_pl     =   zeros(Ncy+1)
    ε̇xxd_pl    =   zeros(Ncy+1)
    ε̇yyd_pl    =   zeros(Ncy+1)
    ε̇zzd_pl    =   zeros(Ncy+1)
    ∇v         =   zeros(Ncy+1)
    ∇v_pl      =   zeros(Ncy+1)

    ω̇z        .= -ε0/2

    # Monitoring
    probes    = (t = zeros(Nt), Ẇ0 = zeros(Nt), τxy0 = zeros(Nt), σyy0 = zeros(Nt), Vx0 = zeros(Nt), τii, εyy = zeros(Nt), σxx=zeros(Nt), fric=zeros(Nt), θs3 = zeros(Nt), 
                θs3_out = zeros(Nt), θs3_in = zeros(Nt), fric_in = zeros(Nt), fric_out = zeros(Nt), 
                σxx_in = zeros(Nt), σxx_out = zeros(Nt), εyy_in = zeros(Nt), εyy_out = zeros(Nt), γxy_in = zeros(Nt), γxy_out = zeros(Nt),
                sdotrat_in = zeros(Nt), sdotrat_out = zeros(Nt))
    η        .= μs
   
    # BC
    BC_Vy = :Neumann
    # BC_Vy = :Dirichlet
    VxS   =  ε0*yv[1]
    VxN   =  ε0*yv[end]
    VyS   =  0.0
    VyN   =  0.0

    # PT solver
    niter = 40000
    θVx   = 0.6
    θVy   = 0.6
    θω̇z   = 0.6
    θPt   = 1.0
    nout  = 1000
    ϵ     = 1e-12
    rel   = 1e-2
    errPt, errVx, errVy = 0., 0., 0.

    τxy_v  = zeros(Ncy+1)
    τxy_v0 = zeros(Ncy+1)
    τ̇xy_v  = zeros(Ncy+1)
    τ̇xy_e  = zeros(Ncy+1)

    Δt_ref = Δt
    Δt    /= 1e5
    Δt0    = Δt
    time   = 0.

    # anim = @animate for it=1:Nt
    for it=-1:Nt-1

        # History
        Δt00 = Δt0
        Δt0  = Δt
        if it>0
            Δt   = Δt_ref 
        end
        time += Δt

        @. τxy0000 = τxy000
        @. τxx0000 = τxx000
        @. τyy0000 = τyy000
        @. τzz0000 = τzz000
        @. Pt0000  = Pt000
        @. myz0000 = myz000
        @. Rz0000  = Rz000

        @. τxy000 = τxy00
        @. τxx000 = τxx00
        @. τyy000 = τyy00
        @. τzz000 = τzz00
        @. Pt000  = Pt00
        @. myz000 = myz00
        @. Rz000  = Rz00

        @. τxy00  = τxy0
        @. τxx00  = τxx0
        @. τyy00  = τyy0
        @. τzz00  = τzz0
        @. Pt00   = Pt0
        @. myz00  = myz0
        @. Rz00   = Rz0

        @. τxy0   = τxy
        @. τxx0   = τxx
        @. τyy0   = τyy
        @. τzz0   = τzz
        @. Pt0    = Pt
        @. myz0   = myz
        @. Rz0    = Rz

        @. λ̇      = 0.0
        @. λ̇rel   = 0.0
        @. τxy_v0 = τxy_v

        if params.bdf==1 || it<1
            a = 1/Δt
            b =-1/Δt
            c = 0.
            d = 0.
            e = 0.
        end

        if params.bdf==2 && it>0
            a, b, c = bdf2(Δt0, Δt)
            d = 0.
            e = 0.
        end

        if params.bdf==3 && it>0
            a, b, c, d = bdf3(Δt00, Δt0, Δt)
            e = 0.
        end

        ηe  .= ( (G/a).^(-1) .+ 1.0./η ).^(-1)
        ηve .= ηe

        # nout= 1
        @views for iter=1:niter

            # Kinematics
            # Vx BC
            Vx[1]   = -Vx[2]     + 2VxS
            Vx[end] = -Vx[end-1] + 2VxN
            # Vy BC
            Vy[1]   = - Vy[2]     + 2VyS
            if BC_Vy == :Dirichlet
                Vy[end] = - Vy[end-1] + 2VyN
            elseif BC_Vy == :Neumann   
                #  Vy[end] = (3.0 * Kb .* Vy[end-1] .* Δt .* ηe[end] + 3.0 * Kb .* Δt .* Δy .* ηe[end] .* ∇v_pl[end] + 3.0 * Pt0[end] .* Δy .* ηe[end] + 4.0 * Vy[end-1] .* ηe[end] .* ηve[end] + 6.0 * Δy .* ε̇yy_pl[end] .* ηe[end] .* ηve[end] + 3.0 * Δy .* ηe[end] .* σyyi - 3.0 * Δy .* ηve[end] .* τyy0[end]) ./ (ηe[end] .* (3.0 * Kb .* Δt + 4.0 * ηve[end]))
                #  Vy[end] = (3.0 * G .* Kb .* Vy[end-1] + 3.0 * G .* Kb .* Δy .* ∇v_pl[end] - 3.0 * G .* Pt00[end] .* c .* Δy - 3.0 * G .* Pt0[end] .* b .* Δy + 4.0 * G .* Vy[end-1] .* a .* ηve[end] + 6.0 * G .* a .* Δy .* ε̇yy_pl[end] .* ηve[end] + 3.0 * G .* a .* Δy .* σyyi + 3.0 * a .* b .* Δy .* ηve[end] .* τyy0[end] + 3.0 * a .* c .* Δy .* ηve[end] .* τyy00[end]) ./ (G .* (3.0 * Kb + 4.0 * a .* ηve[end]))
                Vy[end] = (3.0 * G .* Kb .* Vy[end-1] + 2.0 * G .* Kb .* Δy .* ∇v_pl[end] - 2.0 * G .* Pt000[end] .* d .* Δy - 2.0 * G .* Pt00[end] .* c .* Δy - 2.0 * G .* Pt0[end] .* b .* Δy + 2.0 * G .* Vy[end-1] .* a .* ηve[end] + 4.0 * G .* a .* Δy .* ε̇yy_pl[end] .* ηve[end] + 2.0 * G .* a .* Δy .* σyyi + 2.0 * a .* b .* Δy .* ηve[end] .* τyy0[end] + 2.0 * a .* c .* Δy .* ηve[end] .* τyy00[end] + 2.0 * a .* d .* Δy .* ηve[end] .* τyy000[end]) ./ (G .* (3.0 * Kb + 2.0 * a .* ηve[end]))
            end
            ω̇z[1]   = ω̇z[2]
            ω̇z[end] = ω̇z[end-1]
           
            # Total strain rates
            @. ε̇xy          =  0.5*(Vx[2:end] - Vx[1:end-1])/Δy
            @. ε̇yy          =      (Vy[2:end] - Vy[1:end-1])/Δy # total
            @. ε̇zz          = 1/2*(ε̇xx + ε̇yy)
            @. ∇v           = ε̇xx + ε̇yy + ε̇zz
            @. ω̇zv          =  0.5*(ω̇z[2:end] + ω̇z[1:end-1])
            @. Ẇz           =  0.5*(Vx[2:end] - Vx[1:end-1])/Δy + ω̇zv
            @. κ̇yz          =      (ω̇z[2:end] - ω̇z[1:end-1])/Δy
    
            # Deviatoric strain rates
            @. ε̇xxd = ε̇xx - 1/3*∇v
            @. ε̇yyd = ε̇yy - 1/3*∇v
            @. ε̇zzd = ε̇zz - 1/3*∇v

            # Stress
            @. τxx     =      2*ηve * ( ε̇xxd     - (b*τxx0 + c*τxx00 + d*τxx000 + e*τxx0000)/(2*G))
            @. τyy     =      2*ηve * ( ε̇yyd     - (b*τyy0 + c*τyy00 + d*τyy000 + e*τyy0000)/(2*G))
            @. τzz     =      2*ηve * ( ε̇zzd     - (b*τzz0 + c*τzz00 + d*τzz000 + e*τzz0000)/(2*G))
            @. τxy     =      2*ηve * ( ε̇xy      - (b*τxy0 + c*τxy00 + d*τxy000 + e*τxy0000)/(2*G)) 
            @. Rz      =      2*ηve * ( -Ẇz      - (b* Rz0 + c* Rz00 + d* Rz000 + e* Rz0000)/(2*G))
            @. myz     = lc^2*2*ηve * ( κ̇yz      - (b*myz0 + c*myz00 + d*myz000 + e*myz0000)/(2*G*lc^2)) 
            @. τii     = sqrt(τxy^2 + coss*Rz^2 + 0.5*(τyy^2 + τxx^2 + τzz^2 + coss*(myz/lc)^2))

            @. Ptc     =       Kb/a * (-(∇v - ∇v_pl) - (b* Pt0 + c*Pt00  + d*Pt000 + e*Pt0000)/Kb )
            # @show (-Ptc[end] + τyy[end])*sc.σ

            @. ηvep = ηve
            @. Ptc  = Pt
            # # Plasticity
            # @. F    = τii - Coh*cos(ϕ) - Pt*sin(ϕ)
            # @. ispl = 0
            # @. ispl[F>=0] = 1
            # @. ε̇iiᵉᶠᶠ   = sqrt( (ε̇xy + τxy0/2/ηe)^2 + 0.5*( (ε̇xxd + τxx0/(2*ηe))^2 + ((ε̇yyd + τyy0/(2*ηe))).^2 + ((ε̇zzd + τzz0/(2*ηe))).^2 ) ) 
            # for itpl=1:500  
            #     @. ηvep   = (Coh*cos(ϕ) + Ptc*sin(ϕ) + ηvp*λ̇rel) / 2.0 / ε̇iiᵉᶠᶠ
            #     @. ε̇xxd_pl = λ̇rel*(τxx/2/τii)
            #     @. ε̇yyd_pl = λ̇rel*(τyy/2/τii)
            #     @. ε̇zzd_pl = λ̇rel*(τzz/2/τii) 
            #     @. ε̇xy_pl  = λ̇rel*(τxy/2/τii)
            #     @. ẇz_pl   = λ̇rel*( Rz/2/τii) 
            #     @. κ̇yz_pl  = λ̇rel*(myz/2/τii)/lc^2
            #     @. ∇v_pl   = sin(ψ)*λ̇rel
            #     if params.oop == :Vermeer1990
            #         @. ε̇zz_pl  = λ̇rel*(τxx/2/τii + τyy/2/τii)/2   # dqdτzz*λ̇rel
            #         @. ∇v_pl   = 3/2*sin(ψ)*λ̇rel
            #     end
            #     @. Ptc     =       Kb/a * (-(∇v - ∇v_pl) - (b* Pt0 + c*Pt00  + d*Pt000 + e*Pt0000)/Kb )
            #     @. τxx     =      2*ηve * (  ε̇xxd     - (b*τxx0 + c*τxx00 + d*τxx000 + e*τxx0000)/(2*G)      - ε̇xxd_pl)
            #     @. τyy     =      2*ηve * (  ε̇yyd     - (b*τyy0 + c*τyy00 + d*τyy000 + e*τyy0000)/(2*G)      - ε̇yyd_pl)
            #     @. τzz     =      2*ηve * (  ε̇zzd     - (b*τzz0 + c*τzz00 + d*τzz000 + e*τzz0000)/(2*G)      - ε̇zzd_pl)
            #     @. τxy     =      2*ηve * (  ε̇xy      - (b*τxy0 + c*τxy00 + d*τxy000 + e*τxy0000)/(2*G)      - ε̇xy_pl ) 
            #     @. Rz      =      2*ηve * (  -Ẇz      - (b* Rz0 + c* Rz00 + d* Rz000 + e* Rz0000)/(2*G)      + ẇz_pl  )
            #     @. myz     = lc^2*2*ηve * (  κ̇yz      - (b*myz0 + c*myz00 + d*myz000 + e*myz0000)/(2*G*lc^2) - κ̇yz_pl )
            #     @. τii    = sqrt(τxy^2 + coss*Rz^2 + 0.5*(τyy^2 + τxx^2 + τzz^2 + coss*(myz/lc)^2))
            #     @. Fc     = τii - Coh*cos(ϕ) - Ptc*sin(ϕ) - ηvp*λ̇rel
            #     @. λ̇rel  += (F.>0) .* Fc / (ηvp + ηve + Kb*Δt*sin(ϕ)*sin(ψ)) 
            #     if maximum(Fc) < ϵ break end 
            # end 

            @. τxy_v    = ηvp*ε̇xy_pl   
            @. τ̇xy_e    = 2*ηve * ε̇xy #+ τxy0/(2*ηe)
            @. τ̇xy_v    = (τxy_v - τxy_v0)/Δt
            
            # PT time steps
            @. η_mm  = min.(ηve[1:end-1], ηve[2:end]); 
            @. ΔτV   = Δy^2/(η_mm)/2.1 /4
            @. ΔτPt  = 3/2
            @. Δτω̇z  = Δy/η_mm/4.1/4000
            
            # Residuals
            @. RPt          =  (- Kb*∇v/a - (Pt + b/a*Pt0 + c/a*Pt00 + d/a*Pt000 + e/a*Pt0000))
            @. RVx[2:end-1] =  ((τxy[2:end] - τxy[1:end-1])/Δy - coss*( Rz[2:end] -  Rz[1:end-1])/Δy )
            @. RVy[2:end-1] =  ((τyy[2:end] - τyy[1:end-1])/Δy - (Ptc[2:end] - Ptc[1:end-1])/Δy)
            @. Rzc[2:end-1] = 0.5*(Rz[1:end-1]) + 0.5*(Rz[2:end-0])
            @. Rω̇z[2:end-1] =  coss*( (myz[2:end] - myz[1:end-1])/Δy + 2*Rzc[2:end-1])

            # Damp residuals
            @. ∂Vx∂τ = RVx + (1.0 - θVx)*∂Vx∂τ
            @. ∂Vy∂τ = RVy + (1.0 - θVy)*∂Vy∂τ
            @. ∂Pt∂τ = RPt + (1.0 - θPt)*∂Pt∂τ
            @. ∂ω̇z∂τ = Rω̇z + (1.0 - θω̇z)*∂ω̇z∂τ

            # Update solutions
            @. Vx[2:end-1] += ΔτV  * ∂Vx∂τ[2:end-1] 
            @. Vy[2:end-1] += ΔτV  * ∂Vy∂τ[2:end-1] 
            @. Pt          += ΔτPt * ∂Pt∂τ 
            @. ω̇z[2:end-1] += Δτω̇z * ∂ω̇z∂τ[2:end-1]

            if mod(iter, nout) == 0 || iter==1
                errPt = norm(RPt)/sqrt(length(RPt))
                errVx = norm(RVx)/sqrt(length(RVx))
                errVy = norm(RVy)/sqrt(length(RVy))
                errω̇z = norm(Rω̇z)/sqrt(length(Rω̇z))
                σyyBC = τyy[end] - Ptc[end]
                if noisy
                    @printf("Iteration %05d --- Time step %4d --- σyyBC = %2.7e --- max(F) = %2.2e --- max(Fc) = %2.2e \n", iter, it, σyyBC.*sc.σ/1e3, maximum(F.*sc.σ), maximum(Fc.*sc.σ) )
                    @printf("fPt = %2.4e\n", errPt)
                    @printf("fVx = %2.4e\n", errVx)
                    @printf("fVy = %2.4e\n", errVy)
                    @printf("fω̇z = %2.4e\n", errω̇z)
                end
                (errPt < ϵ && errω̇z < ϵ && errVx < ϵ && errVy < ϵ) && break 
                (isnan(errPt) || isnan(errVx) || isnan(errVx)) && error("NaNs")        
            end
        end        

        @. Pt   = Ptc
        @. εxy += ε̇xy*Δt 
        @. εyy += ε̇yy*Δt 

        # Show array infos
        if noisy
            @minmax(Pt)
            @minmax(τxx)
            @minmax(τyy)
            @minmax(τzz)
            @minmax(τxy)
        end

        # if (errPt > ϵ || errVx > ϵ || errVx > ϵ) error("non converged") end
        PrincipalStress!(σ1, σ3, τxx, τyy, τzz, τxy, Pt)

        # Probe model state
        _, iA = findmax(Pt)
        _, iB = findmin(Pt)

        if it>0
            it1 = it+1
            probes.t[it1]        = time
            probes.Ẇ0[it1]       = τxy[end]*ε̇xy[end]
            probes.τxy0[it1]     = τxy[end]*sc.σ
            probes.Vx0[it1]      = 0.5*(Vx[end] + Vx[end-1])
            probes.σyy0[it1]     = τyy[end] - Pt[end]
            probes.θs3_out[it1]  = atand(σ3.z[iA] ./ σ3.x[iA])
            probes.θs3_in[it1]   = atand(σ3.z[iB] ./ σ3.x[iB])
            probes.fric_out[it1] = -τxy[iA]./(τyy[iA] .- Pt[iA])
            probes.fric_in[it1]  = -τxy[iB]./(τyy[iB] .- Pt[iB])
            probes.fric[it1]     = -τxy[iB]./(τyy[iB] .- Pt[iB])
            probes.εyy[it1]      = εyy[iB]
            probes.εyy_out[it1]  = εyy[iA]
            probes.εyy_in[it1]   = εyy[iB]
            probes.θs3[it1]      = atand.(σ3.z[iB] ./ σ3.x[iB])
            probes.σxx[it1]      = (τxx[iB]-Pt[iB])*sc.σ
            probes.σxx_out[it1]  = (τxx[iA]-Pt[iA])*sc.σ
            probes.σxx_in[it1]   = (τxx[iB]-Pt[iB])*sc.σ
            probes.γxy_out[it1]  = εxy[iA]
            probes.γxy_in[it1]   = εxy[iB] 
            probes.sdotrat_out[it1] = τ̇xy_v[iA]*sc.σ/sc.t
            probes.sdotrat_in[it1]  = τ̇xy_v[iB]*sc.σ/sc.t
        end


        # Visualisation
        if visu==true && it==Nt-1

                    
            τxy_ana  = G.*ε0.*probes.t
            τxy_ana  = η[1]*ε0*(1 .- exp.(-G/η[1].*probes.t))*sc.σ
            fric_ana = (η[1]*ε0*(1 .- exp.(-G/η[1].*probes.t)))*sc.σ./100e3
    
            p4 = plot( xlabel = "γxy BC [%]", ylabel = "σxy", foreground_color_legend = nothing, background_color_legend = nothing )
            p4 = plot!((1:it)*ε0*Δt*100, probes.fric[1:it], label="σxy", color=:blue )
            p4 = plot!((1:it)*ε0*Δt*100, fric_ana[1:it], label="analytics", color=:red )
            # p4 = plot!(title=@sprintf("-σxy/σyy = %1.4f", probes.fric_in[it]))
            display(plot(p4))

            @show norm(probes.fric[1:it] .- fric_ana[1:it]) / norm(fric_ana)
            @show norm(probes.τxy0[1:it] .- τxy_ana[1:it]) / norm(τxy_ana)

        end
    end
            
    τxy_ana  = G.*ε0.*probes.t
    τxy_ana  = η[1]*ε0*(1 .- exp.(-G/η[1].*probes.t))*sc.σ

    if make_gif gif(anim, "figures/Test1_MohrCircles_1D.gif", fps = 15) end
    return (γxy=(0:Nt-1)*ε0*Δt*100, εyy=probes.εyy.*100, app_fric=probes.fric, σxx=.-probes.σxx./1e3, θ=probes.θs3, σxx_out=probes.σxx_out/1e3, σxx_in=probes.σxx_in/1e3, θ_out=probes.θs3_out, θ_in=probes.θs3_in, γxy_out=probes.γxy_out.*100, γxy_in=probes.γxy_in.*100, εyy_out=probes.εyy_out.*100, εyy_in=probes.εyy_in.*100, srat_in=probes.sdotrat_out, srat_out=probes.sdotrat_out)
end

#  Case B
σi       = (xx = -400e3, yy=-100e3, xy=0.0)
Main_VEP_1D_test_BDF(σi; visu=true);