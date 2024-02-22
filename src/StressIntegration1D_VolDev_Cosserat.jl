
using PlasticityTests, Plots, Printf, LinearAlgebra
function Main_VEP_1D_vdev_coss(σi; params=(
    #---------------#
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
    pl  = true), 
    #---------------#
    visu     = true, 
    make_gif = false,  
    Ncy      = 40   # default parameter set
    #---------------#
    )
    sc = (σ = params.G, L = 1.0, t = 1.0/params.γ̇xy)

    coss  = true

    # Visualisation is important!
    if visu==false 
        make_gif = false 
    end
    
    # Load digitised data from test 1 of Vermeer (1990)
    σxxi     = σi.xx/sc.σ # Courbe A - Vermeer
    σyyi     = σi.yy/sc.σ # Courbe A - Vermeer
    if σi.xx>σi.xx
        Vermeer1990 = ExtractDataCase("CaseA")
    else
        Vermeer1990 = ExtractDataCase("CaseB")
    end

    # Physical parameters
    σzzi       = 0.5*(σxxi + σyyi)
    Pi         = -(σxxi + σyyi)/2.0
    τxxi       = Pi + σxxi
    τyyi       = Pi + σyyi
    τzzi       = Pi + σzzi
    τxyi       = 0.0

    # Dimensionally dependent
    Ly         = 2e4/sc.L
    Ẇ0         = 5e-5/(sc.σ/sc.t)
    ε0         = params.γ̇xy/(1.0/sc.t)
    G          = params.G/sc.σ
    Kb         = params.K/sc.σ
    Coh0       = params.c/sc.σ
    Coh1       =      1.0/sc.σ
    μs         = 1e52/(sc.σ*sc.t)
    ηvp        = params.ηvp/(sc.σ*sc.t)
    ϕ          = params.ϕ
    ψ          = params.ψ   
    Gc         = G
    lc         = 5e2/sc.L

    # Numerical parameters
    Nt         = params.nt-1
    Δy         = Ly/Ncy
    yc         = LinRange(-Ly/2-Δy/2, Ly/2+Δy/2, Ncy+2)
    yv         = LinRange(-Ly/2,      Ly/2,      Ncy+1)
    Δt         = params.Δt/sc.t
    nout_viz   = 10
    if make_gif 
        nout_viz = 1 
    end

    # Allocate arrays
    ηe         = G*Δt*ones(Ncy+1); #ηe[Int64(ceil(Ncy/2))] = G*Δt/5
    Pt         =  Pi*ones(Ncy+1) 
    Ptc        =  Pi*ones(Ncy+1) 
    Pt0        =  Pi*ones(Ncy+1)
    τxx        =  τxxi*ones(Ncy+1)
    τyy        =  τyyi*ones(Ncy+1)    
    τzz        =  τzzi*ones(Ncy+1) 
    τxy        =  τxyi*ones(Ncy+1)
    σxy        =  zeros(Ncy) # will need averaging
    σyx        =  zeros(Ncy)
    τxx0       =  τxxi*ones(Ncy+1)
    τyy0       =  τyyi*ones(Ncy+1)
    τzz0       =  τzzi*ones(Ncy+1)
    τxy0       =  τxyi*ones(Ncy+1)
    Coh        =  ones((Ncy+1)).*Coh1; 
    Coh[Int64(ceil(Ncy/2))] = Coh0  
    F          =  zeros((Ncy+1))
    Fc         =  zeros((Ncy+1))
    λ̇          =  zeros((Ncy+1))
    λ̇rel       =  zeros((Ncy+1))
    ispl       =  zeros(Int, (Ncy+1))
    ηve        =  zeros((Ncy+1)); ηve .= ηe
    ηvep       =  zeros((Ncy+1))
    εxy        =    zeros(Ncy+1)
    # ẇz         =    zeros(Ncy+1)
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
    probes    = (Ẇ0 = zeros(Nt), τxy0 = zeros(Nt), σyy0 = zeros(Nt), Vx0 = zeros(Nt), τii, εyy = zeros(Nt), σxx=zeros(Nt), fric=zeros(Nt), θs3 = zeros(Nt), 
                θs3_out = zeros(Nt), θs3_in = zeros(Nt), fric_in = zeros(Nt), fric_out = zeros(Nt), 
                σxx_in = zeros(Nt), σxx_out = zeros(Nt), γxy_in = zeros(Nt), γxy_out = zeros(Nt))
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
    θω̇z   = 0.6
    θPt   = 1.0
    nout  = 1000
    ϵ     = 1e-10
    rel   = 1e-2
    errPt, errVx, errVy = 0., 0., 0.

    # anim = @animate for it=1:Nt
    for it=1:Nt
        # History
        @. τxy0  = τxy
        @. τxx0  = τxx
        @. τyy0  = τyy
        @. τzz0  = τzz
        @. Pt0   = Pt
        @. myz0  = myz
        @. Rz0   = Rz
        @. λ̇     = 0.0
        @. λ̇rel  = 0.0

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
                Vy[end] = (3.0 * Kb .* Vy[end-1] .* Δt .* ηe[end] + 3.0 * Kb .* Δt .* Δy .* ηe[end] .* ∇v_pl[end] + 3.0 * Pt0[end] .* Δy .* ηe[end] + 4.0 * Vy[end-1] .* ηe[end] .* ηve[end] + 6.0 * Δy .* ε̇yy_pl[end] .* ηe[end] .* ηve[end] + 3.0 * Δy .* ηe[end] .* σyyi - 3.0 * Δy .* ηve[end] .* τyy0[end]) ./ (ηe[end] .* (3.0 * Kb .* Δt + 4.0 * ηve[end]))
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
            @. τxx     =      2*ηve * ( ε̇xxd     + τxx0/(2*ηe))
            @. τyy     =      2*ηve * ( ε̇yyd     + τyy0/(2*ηe))
            @. τzz     =      2*ηve * ( ε̇zzd     + τzz0/(2*ηe))
            @. τxy     =      2*ηve * ( ε̇xy      + τxy0/(2*ηe)) 
            @. Rz      =     -2*ηve * ( Ẇz       -  Rz0/(2*ηe))
            @. myz     = lc^2*2*ηve * ( κ̇yz      + myz0/(2*ηe*lc^2)) 
            @. τii     = sqrt(τxy^2 + coss*Rz^2 + 0.5*(τyy^2 + τxx^2 + τzz^2 + coss*(myz/lc)^2))

            @. ηvep = ηve
            # Plasticity
            @. F    = τii - Coh*cos(ϕ) - Pt*sin(ϕ)
            @. Ptc  = Pt
            @. ispl = 0
            @. ispl[F>=0] = 1
            @. ε̇iiᵉᶠᶠ   = sqrt( (ε̇xy + τxy0/2/ηe)^2 + 0.5*( (ε̇xxd + τxx0/(2*ηe))^2 + ((ε̇yyd + τyy0/(2*ηe))).^2 + ((ε̇zzd + τzz0/(2*ηe))).^2 ) ) 
            for itpl=1:500  
                @. ηvep   = (Coh*cos(ϕ) + Ptc*sin(ϕ) + ηvp*λ̇rel) / 2.0 / ε̇iiᵉᶠᶠ
                @. ε̇xxd_pl = λ̇rel*(τxx/2/τii)
                @. ε̇yyd_pl = λ̇rel*(τyy/2/τii)
                @. ε̇zzd_pl = λ̇rel*(τzz/2/τii) 
                @. ε̇xy_pl  = λ̇rel*(τxy/2/τii)
                @. ẇz_pl   = λ̇rel*( Rz/2/τii) 
                @. κ̇yz_pl  = λ̇rel*(myz/2/τii)/lc^2
                @. ∇v_pl   = sin(ψ)*λ̇rel
                if params.oop == :Vermeer1990
                    @. ε̇zz_pl  = λ̇rel*(τxx/2/τii + τyy/2/τii)/2   # dqdτzz*λ̇rel
                    @. ∇v_pl   = 3/2*sin(ψ)*λ̇rel
                end
                @. Ptc    = Pt0  - Kb*Δt*(∇v - ∇v_pl)
                @. τxx    =  2*ηve * ( ε̇xxd     + τxx0/(2*ηe) -  ε̇xxd_pl)
                @. τyy    =  2*ηve * ( ε̇yyd     + τyy0/(2*ηe) -  ε̇yyd_pl)
                @. τzz    =  2*ηve * ( ε̇zzd     + τzz0/(2*ηe) -  ε̇zzd_pl)
                @. τxy    =  2*ηve * ( ε̇xy      + τxy0/(2*ηe) -  ε̇xy_pl ) 
                @. Rz     = -2*ηve * ( Ẇz       -  Rz0/(2*ηe) -  ẇz_pl  )
                @. myz    = lc^2*2*ηve*( κ̇yz + myz0/(2*ηe*lc^2) - κ̇yz_pl)
                @. τii    = sqrt(τxy^2 + coss*Rz^2 + 0.5*(τyy^2 + τxx^2 + τzz^2 + coss*(myz/lc)^2))
                @. Fc     = τii - Coh*cos(ϕ) - Ptc*sin(ϕ) - ηvp*λ̇rel
                @. λ̇rel  += (F.>0) .* Fc / (ηvp + ηve + Kb*Δt*sin(ϕ)*sin(ψ)) 
                if maximum(Fc) < ϵ break end 
            end 
          
            # PT time steps
            @. η_mm  = min.(ηve[1:end-1], ηve[2:end]); 
            @. ΔτV   = Δy^2/(η_mm)/2.1 /4 
            @. ΔτPt  = 3/2
            @. Δτω̇z  = Δy/η_mm/4.1/4000
            
            # Residuals
            @. RPt          =  (- Kb*Δt*∇v - (Pt - Pt0))
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
                @printf("Iteration %05d --- Time step %4d --- σyyBC = %2.7e --- max(F) = %2.2e --- max(Fc) = %2.2e \n", iter, it, σyyBC.*sc.σ/1e3, maximum(F.*sc.σ), maximum(Fc.*sc.σ) )
                @printf("fPt = %2.4e\n", errPt)
                @printf("fVx = %2.4e\n", errVx)
                @printf("fVy = %2.4e\n", errVy)
                @printf("fω̇z = %2.4e\n", errω̇z)
                (errPt < ϵ && errω̇z < ϵ && errVx < ϵ && errVy < ϵ) && break 
                (isnan(errPt) || isnan(errVx) || isnan(errVx)) && error("NaNs")        
            end
        end        

        @. Pt   = Ptc
        @. εxy += ε̇xy*Δt 
        @. εyy += ε̇yy*Δt 

        # Show array infos
        @minmax(Pt)
        @minmax(τxx)
        @minmax(τyy)
        @minmax(τzz)
        @minmax(τxy)

        # if (errPt > ϵ || errVx > ϵ || errVx > ϵ) error("non converged") end
        PrincipalStress!(σ1, σ3, τxx, τyy, τzz, τxy, Pt)

        # Probe model state
        _, iA = findmax(Pt)
        _, iB = findmin(Pt)

        probes.Ẇ0[it]       = τxy[end]*ε̇xy[end]
        probes.τxy0[it]     = τxy[end]
        probes.Vx0[it]      = 0.5*(Vx[end] + Vx[end-1])
        probes.σyy0[it]     = τyy[end] - Pt[end]
        probes.θs3_out[it]  = atand(σ3.z[iA] ./ σ3.x[iA])
        probes.θs3_in[it]   = atand(σ3.z[iB] ./ σ3.x[iB])
        probes.fric_out[it] = -τxy[iA]./(τyy[iA] .- Pt[iA])
        probes.fric_in[it]  = -τxy[iB]./(τyy[iB] .- Pt[iB])
        probes.εyy[it]      = εyy[iB]
        probes.fric[it]     = -τxy[iB]./(τyy[iB] .- Pt[iB])
        probes.θs3[it]      = atand.(σ3.z[iB] ./ σ3.x[iB])
        probes.σxx[it]      = (τxx[iB]-Pt[iB])*sc.σ
        probes.σxx_out[it]  = (τxx[iA]-Pt[iA])*sc.σ
        probes.σxx_in[it]   = (τxx[iB]-Pt[iB])*sc.σ
        probes.γxy_out[it]  = εxy[iA]
        probes.γxy_in[it]   = εxy[iB] 
        
        # Visualisation
        if visu==true && (mod(it, nout_viz)==0 || it==1 || it==Nt)
    
            θ   = LinRange(-π, 0, 100)
            σMC = LinRange(-500, 0, 100 ) .*1e3
            τMC = -σMC.*tan(ϕ) 

            PA    = ((σ1.v[iA] + σ3.v[iA])/2)*sc.σ
            τA    = ((σ1.v[iA] - σ3.v[iA])/2)*sc.σ
            τB    = ((σ1.v[iB] - σ3.v[iB])/2)*sc.σ
            PB    = ((σ1.v[iB] + σ3.v[iB])/2)*sc.σ
            
            yield = (x = σMC./1e3, y = τMC./1e3)
            MC_A  = (x = (τA.*cos.(θ) .+ PA)./1e3, y = (τA.*sin.(θ))./1e3) 
            MC_B  = (x = (τB.*cos.(θ) .+ PB)./1e3, y = (τB.*sin.(θ))./1e3)

            p1 = plot( title = "γxy [%]", xlabel = "γxy [%]", ylabel = "y [-]", legend=:none )
            p1 = plot!(2 .* εxy*100,          yv.*sc.L, label="γxy" )
            p1 = scatter!(2 .* εxy[ispl.==1]*100, yv[ispl.==1].*sc.L, linewidth=5 )

            p2 = plot(title="Mohr circles", ylabel="τ [kPa]", xlabel="σₙ [kPa]", size=(300,300), aspect_ratio=1, xlim=(-500,0), ylim=(0,400))
            p2 = plot!( MC_A... , color=:blue, label="out" )
            p2 = plot!( MC_B...,  color=:green, label="in"  )
            p2 = plot!( yield..., color=:red, label="Yield"  )
            p2 = plot!(ylabel="τ [kPa]", xlabel="σₙ [kPa]", foreground_color_legend = nothing, background_color_legend = nothing, legend=:topright)

            p3 = plot(title = "Stress orientation", ylabel = "θ σ3 [ᵒ]", xlabel = "γxy BC [%]", xlims=(0,8), foreground_color_legend = nothing, background_color_legend = nothing  )
            p3 = scatter!(Vermeer1990.θ.x, Vermeer1990.θ.y, label="Vermeer (1990)", legend=:topright)
            p3 = plot!((1:it)*ε0*Δt*100, probes.θs3_out[1:it], label="out", color=:blue  )
            p3 = plot!((1:it)*ε0*Δt*100, probes.θs3_in[1:it],  label="in" , color=:green )

            p4 = plot( xlabel = "γxy BC [%]", ylabel = "-σxy/σyy", foreground_color_legend = nothing, background_color_legend = nothing )
            scatter!(Vermeer1990.Friction.x, Vermeer1990.Friction.y, label="Vermeer (1990)", legend=:bottomright)
            p4 = plot!((1:it)*ε0*Δt*100, probes.fric_out[1:it], label="out", color=:blue,  xlims=(0,8) )
            p4 = plot!((1:it)*ε0*Δt*100, probes.fric_in[1:it],  label="in" , color=:green, xlims=(0,8) )
            p4 = plot!(title=@sprintf("-σxy/σyy = %1.4f", probes.fric_in[it]))
            display(plot(p1,p2,p3,p4))
            
        end
    end
    if make_gif gif(anim, "figures/Test1_MohrCircles_1D.gif", fps = 15) end
    return (γxy=(0:Nt-1)*ε0*Δt*100, εyy=probes.εyy.*100, app_fric=probes.fric, σxx=.-probes.σxx./1e3, θ=probes.θs3)
end

 # Case B
 σi       = (xx = -400e3, yy=-100e3)

Main_VEP_1D_vdev_coss(σi; visu=true)