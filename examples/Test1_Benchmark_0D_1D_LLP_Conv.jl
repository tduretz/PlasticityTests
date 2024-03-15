using PlasticityTests, Plots, Interpolations, LinearAlgebra, Statistics

function main()

    # Vermeer (1990)
    Gv  = 10e6
    Kv  = 2/3*Gv
    σiA = (xx =  -25e3, yy=-100e3, xy=0.0) # Case A
    σiB = (xx = -400e3, yy=-100e3, xy=0.0) # Case B

    # Resolution
    Ncy     = 10
    fac_vec = [1 2 4]
    Δt_vec  = 20 ./fac_vec 
    L2_bdf1 = zeros(size(Δt_vec))
    L2_bdf2 = zeros(size(Δt_vec))
    L2_bdf3 = zeros(size(Δt_vec))
    L2_bdf4 = zeros(size(Δt_vec))
    nt_ref  = 50

    # Compute analytics with high resolution
    dt_fact = 1000
    params  = (
        K    = Kv,
        G    = Gv,
        c    = 0.0,
        ϕ    = 40/180*π,
        ψ    = 10/180*π,
        θt   = 25/180*π,
        ηvp  = 0*1e7,
        lc   = 1e2,
        γ̇xy  = 0.00001,
        Δt   = 20/dt_fact,
        nt   = nt_ref*dt_fact,
        law  = :DruckerPrager,
        coss = false,
        oop  = :Vermeer1990,
        pl   = true
    ) # default parameter set

    # Analytics. Achtung! Correct thickness for Cosserat
    CaseB_LLP = Vermeer3_ana_llp2013(σiB, params, 10)

    # --------------------- Numerics --------------------- #

    for idt in eachindex(fac_vec)

        # ------------------- BDF 1 ------------------- #
        params = (
            K     = Kv,
            G     = Gv,
            η     = 1e52,
            c     = 0.0,
            ϕ     = 40/180*π,
            ψ     = 10/180*π,
            θt    = 25/180*π,
            ηvp   = 0*1e7,
            lc    = 1e2,
            γ̇xy   = 0.00001,
            Δt    = 20/fac_vec[idt],
            nt    = nt_ref*fac_vec[idt],
            law   = :DruckerPrager,
            coss  = false,
            oop   = :Vermeer1990,
            pl    = true, 
            noisy = false,
            bdf   = 1) # default parameter set

        # Numerics
        @info "BDF1 with dt_fact #$(idt)"
        CaseB_1D  = Main_VEP_1D_vdev_coss_BDF(σiB; params, visu=false, Ncy=Ncy)

        # Interpolate analytics on numerics strain axis
        stp = 1 # scatter step
        app_fric_ana = -CaseB_LLP.σ_in[ 3,1:stp:end] ./ CaseB_LLP.σ_in[ 2,1:stp:end]
        strain_ana   = CaseB_LLP.γ_bulk[1:stp:end]*100
        Itp_fric     = linear_interpolation(strain_ana, app_fric_ana, extrapolation_bc=Line())
        app_fric_ana = Itp_fric.(CaseB_1D.γxy)

        L2_bdf1[idt] = norm(app_fric_ana .- CaseB_1D.app_fric) / norm(app_fric_ana) 

        # ------------------- BDF 2 ------------------- #
        params = (
            K     = Kv,
            G     = Gv,
            η     = 1e52,
            c     = 0.0,
            ϕ     = 40/180*π,
            ψ     = 10/180*π,
            θt    = 25/180*π,
            ηvp   = 0*1e7,
            lc    = 1e2,
            γ̇xy   = 0.00001,
            Δt    = 20/fac_vec[idt],
            nt    = nt_ref*fac_vec[idt],
            law   = :DruckerPrager,
            coss  = false,
            oop   = :Vermeer1990,
            pl    = true, 
            noisy = false,
            bdf   = 2) # default parameter set

        # Numerics
        @info "BDF2 with dt_fact #$(idt)"
        CaseB_1D  = Main_VEP_1D_vdev_coss_BDF(σiB; params, visu=false, Ncy=Ncy)

        # Interpolate analytics on numerics strain axis
        stp = 1 # scatter step
        app_fric_ana = -CaseB_LLP.σ_in[ 3,1:stp:end] ./ CaseB_LLP.σ_in[ 2,1:stp:end]
        strain_ana   = CaseB_LLP.γ_bulk[1:stp:end]*100
        Itp_fric     = linear_interpolation(strain_ana, app_fric_ana, extrapolation_bc=Line())
        app_fric_ana = Itp_fric.(CaseB_1D.γxy)

        L2_bdf2[idt] = norm(app_fric_ana .- CaseB_1D.app_fric) / norm(app_fric_ana)  

        # ------------------- BDF 3 ------------------- #
        params = (
            K     = Kv,
            G     = Gv,
            η     = 1e52,
            c     = 0.0,
            ϕ     = 40/180*π,
            ψ     = 10/180*π,
            θt    = 25/180*π,
            ηvp   = 0*1e7,
            lc    = 1e2,
            γ̇xy   = 0.00001,
            Δt    = 20/fac_vec[idt],
            nt    = nt_ref*fac_vec[idt],
            law   = :DruckerPrager,
            coss  = false,
            oop   = :Vermeer1990,
            pl    = true, 
            noisy = false,
            bdf   = 3) # default parameter set

        # Numerics
        @info "BDF3 with dt_fact #$(idt)"
        CaseB_1D  = Main_VEP_1D_vdev_coss_BDF(σiB; params, visu=false, Ncy=Ncy)

        # Interpolate analytics on numerics strain axis
        stp = 1 # scatter step
        app_fric_ana = -CaseB_LLP.σ_in[ 3,1:stp:end] ./ CaseB_LLP.σ_in[ 2,1:stp:end]
        strain_ana   = CaseB_LLP.γ_bulk[1:stp:end]*100
        Itp_fric     = linear_interpolation(strain_ana, app_fric_ana, extrapolation_bc=Line())
        app_fric_ana = Itp_fric.(CaseB_1D.γxy)
        
        L2_bdf3[idt] = norm(app_fric_ana .- CaseB_1D.app_fric) / norm(app_fric_ana)   

        # ------------------- VISU ------------------- #
  
        # Panel (1,1) - Stress ratio
        p1 = plot(title="Stress ratio", ylabel="-σxy/σyy")
        p1 = scatter!(CaseB_1D.γxy, CaseB_1D.app_fric, color=:green, label="numerics")
        # p1 = plot!(CaseB_LLP.γ_bulk[1:stp:end]*100, -CaseB_LLP.σ_in[ 3,1:stp:end] ./ CaseB_LLP.σ_in[ 2,1:stp:end], label="LLP",  color=:red)
        
        # p1 = scatter!(CaseB_1D.γxy, params.G.*CaseB_1D.γxy./100e5)
        p1 = scatter!(CaseB_1D.γxy, app_fric_ana)

        display(plot(p1))
    end

    @info "Convergence of -τxy/σxy"
    @show L2_bdf1[1:end-1]./L2_bdf1[2:end-0]
    @show L2_bdf2[1:end-1]./L2_bdf2[2:end-0]
    @show L2_bdf3[1:end-1]./L2_bdf3[2:end-0] 

end

main()