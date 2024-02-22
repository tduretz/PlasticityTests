function Main_VEP_1D_tot(σi; params=(
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
    Ncy      = 10   # default parameter set
    #---------------#
    )

    # Visualisation is important!
    if visu==false 
        make_gif = false 
    end

    sc = (σ = params.G, L = 1.0, t = 1.0/params.γ̇xy)

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
    σxyi       = 0.0
 
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

    # Numerical parameters
    Nt         = params.nt-1
    Δy         = Ly/Ncy
    yc         = LinRange(-Ly/2-Δy/2, Ly/2+Δy/2, Ncy+2)
    yv         = LinRange(-Ly/2,      Ly/2,      Ncy+1)
    Δt         = params.Δt/sc.t
    ηe         = G*Δt
    nout_viz   = 10
    if make_gif 
        nout_viz = 1 
    end

    # Allocate arrays
    Pt         =  Pi*ones(Ncy+1) 
    Ptc        =  Pi*ones(Ncy+1) 
    Pt0        =  Pi*ones(Ncy+1)
    τxx        =  τxxi*ones(Ncy+1)
    τxy        =  τxyi*ones(Ncy+1)
    τyy        =  τyyi*ones(Ncy+1)    
    τzz        =  τzzi*ones(Ncy+1)
    σxx        =  σxxi*ones(Ncy+1)
    σxy        =  σxyi*ones(Ncy+1)
    σyy        =  σyyi*ones(Ncy+1)    
    σzz        =  σzzi*ones(Ncy+1)  
    τxx0       =  τxxi*ones(Ncy+1)
    τxy0       =  τxyi*ones(Ncy+1)
    τyy0       =  τyyi*ones(Ncy+1)
    τzz0       =  τzzi*ones(Ncy+1)
    σxx0       =  σxxi*ones(Ncy+1)
    σxy0       =  σxyi*ones(Ncy+1)
    σyy0       =  σyyi*ones(Ncy+1)
    σzz0       =  σzzi*ones(Ncy+1)
    Coh        =  ones((Ncy+1)).*Coh1; 
    Coh[Int64(Ncy/2)] = Coh0  
    F          =  zeros((Ncy+1))
    Fc         =  zeros((Ncy+1))
    λ̇          =  zeros((Ncy+1))
    λ̇rel       =  zeros((Ncy+1))
    ispl       =  zeros(Int, (Ncy+1))
    ηve        =  zeros((Ncy+1)); ηve .= ηe
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
    ε̇xxd_pl    =   zeros(Ncy+1)
    ε̇yyd_pl    =   zeros(Ncy+1)
    ε̇zzd_pl    =   zeros(Ncy+1)
    ∇v         =   zeros(Ncy+1)
    ∇v_pl      =   zeros(Ncy+1)

    # Monitoring
    probes    = (Ẇ0 = zeros(Nt), τxy0 = zeros(Nt), σyy0 = zeros(Nt), Vx0 = zeros(Nt), τii, θs3 = zeros(Nt), θs3_out = zeros(Nt), θs3_in = zeros(Nt), fric_in = zeros(Nt), fric_out = zeros(Nt), εyy = zeros(Nt), σxx=zeros(Nt), fric=zeros(Nt))
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

    # anim = @animate for it=1:Nt
    for it=1:Nt
        # History
        @. τxy0  = τxy
        @. τxx0  = τxx
        @. τyy0  = τyy
        @. τzz0  = τzz

        @. σxy0  = σxy
        @. σxx0  = σxx
        @. σyy0  = σyy
        @. σzz0  = σzz

        @. Pt0   = Pt
        @. λ̇     = 0.0
        @. λ̇rel  = 0.0

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
            
            # Total strains
            @. ε̇xy  =  0.5*(Vx[2:end] - Vx[1:end-1])/Δy
            @. ε̇yy  =      (Vy[2:end] - Vy[1:end-1])/Δy # total
            @. ε̇zz  = 1/2*(ε̇xx + ε̇yy)
            @. ∇v   = ε̇xx + ε̇yy + ε̇zz

            # Stress
            @. σxx     =  2*ηve * (ε̇xx  + σxx0/(2*ηe))
            @. σyy     =  2*ηve * (ε̇yy  + σyy0/(2*ηe))
            @. σzz     =  2*ηve * (ε̇zz  + σzz0/(2*ηe))
            @. σxy     =  2*ηve * (ε̇xy  + σxy0/(2*ηe))
            @. Pt      = -1/3*(σxx + σyy + σzz)
            @. τxy     = σxy
            @. τxx     = σxx + Pt
            @. τyy     = σyy + Pt
            @. τzz     = σzz + Pt
            @. τii     = sqrt(τxy^2 + 0.5*(τyy.^2 + τxx.^2 + τzz.^2))

            # Plasticity
            @. F    = τii - Coh*cos(ϕ) - Pt*sin(ϕ)
            @. Ptc  = Pt
            @. ηvep = ηve
            @. ispl = 0
            @. ispl[F>=0] = 1
            @. ε̇iiᵉᶠᶠ   = sqrt( (ε̇xy + τxy0/2/ηe)^2 + 0.5*( (ε̇xxd + τxx0/(2*ηe))^2 + ((ε̇yyd + τyy0/(2*ηe))).^2 + ((ε̇zzd + τzz0/(2*ηe))).^2 ) )    
            
            for it=1:50  
                @. ηvep    = (Coh*cos(ϕ) + Ptc*sin(ϕ) + ηvp*λ̇rel) / 2.0 / ε̇iiᵉᶠᶠ
                @. ε̇xxd_pl = λ̇rel*(τxx/2/τii)
                @. ε̇yyd_pl = λ̇rel*(τyy/2/τii)
                @. ε̇zzd_pl = λ̇rel*(τzz/2/τii)  
                @. ε̇xy_pl  = λ̇rel*(τxy/2/τii)
                @. ∇v_pl   = sin(ψ)*λ̇rel
                if params.oop == :Vermeer1990
                    @. ε̇zz_pl = λ̇rel*(τxx/2/τii + τyy/2/τii)/2   # dqdτzz*λ̇rel
                    @. ∇v_pl  = 3/2*sin(ψ)*λ̇rel
                end
      
                @. σxx     =  2*ηve * (ε̇xx  + σxx0/(2*ηe) - (ε̇xxd_pl + 1/3*∇v_pl))
                @. σyy     =  2*ηve * (ε̇yy  + σyy0/(2*ηe) - (ε̇yyd_pl + 1/3*∇v_pl))
                @. σzz     =  2*ηve * (ε̇zz  + σzz0/(2*ηe) - (ε̇zzd_pl + 1/3*∇v_pl))
                @. σxy     =  2*ηve * (ε̇xy  + σxy0/(2*ηe) -  ε̇xy_pl)
                @. Pt      = -1/3*(σxx + σyy + σzz)
                @. τxy     = σxy
                @. τxx     = σxx + Pt
                @. τyy     = σyy + Pt
                @. τzz     = σzz + Pt
                @. τii     = sqrt(τxy^2 + 0.5*(τyy.^2 + τxx.^2 + τzz.^2))

                @. Fc     = τii - Coh*cos(ϕ) - Pt*sin(ϕ) - ηvp*λ̇rel
                @. λ̇rel  += (F.>0) .* Fc / (ηvp + ηve + Kb*Δt*sin(ϕ)*sin(ψ))
                if maximum(Fc) < ϵ break end 
            end

            # PT time steps
            @. η_mm = min.(ηvep[1:end-1], ηvep[2:end]); 
            @. ΔτV   = Δy^2/(η_mm)/2.1 /4 /10
            @. ΔτPt  = ηvep/Δy/G/Δt/3/10
            
            # Residuals
            @. RVx[2:end-1] =  ((σxy[2:end] - σxy[1:end-1])/Δy )
            @. RVy[2:end-1] =  ((σyy[2:end] - σyy[1:end-1])/Δy )
            
            # Damp residuals
            @. ∂Vx∂τ = RVx + (1.0 - θVx)*∂Vx∂τ
            @. ∂Vy∂τ = RVy + (1.0 - θVy)*∂Vy∂τ
            @. ∂Pt∂τ = RPt + (1.0 - θPt)*∂Pt∂τ

            # Update solutions
            @. Vx[2:end-1] += ΔτV * ∂Vx∂τ[2:end-1] 
            @. Vy[2:end-1] += ΔτV * ∂Vy∂τ[2:end-1] 

            if mod(iter, nout) == 0 || iter==1
                errVx = norm(RVx)/sqrt(length(RVx))
                errVy = norm(RVy)/sqrt(length(RVy))
                σyyBC = τyy[end] - Ptc[end]
                @printf("Iteration %05d --- Time step %4d --- σyyBC = %2.7e --- max(F) = %2.2e --- max(Fc) = %2.2e \n", iter, it, σyyBC.*sc.σ/1e3, maximum(F.*sc.σ), maximum(F.*sc.σ) )
                @printf("fVx = %2.4e\n", errVx)
                @printf("fVy = %2.4e\n", errVy)
                (errVx < ϵ && errVy < ϵ) && break 
                (isnan(errVx) || isnan(errVx)) && error("NaNs")        
            end
        end

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
        _, iA = findmin(σ3.v)
        _, iB = findmax(σ1.v)

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