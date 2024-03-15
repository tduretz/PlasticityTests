using PlasticityTests, Plots, Interpolations, LinearAlgebra, Statistics

function main()

    # Vermeer (1990)
    Gv  = 10e6
    Kv  = 2/3*Gv
    σiB = (xx = -400e3, yy=-100e3, xy=0.0) # Case B

    # Resolution
    Ncy     = 10
    fac_vec = [1 2 4]
    Δt_vec  = 20 ./fac_vec 
    nt_ref  = 50

    L2_bdf1_τxy  = zeros(size(Δt_vec))
    L2_bdf2_τxy  = zeros(size(Δt_vec))
    L2_bdf3_τxy  = zeros(size(Δt_vec))
    L2_bdf1_fric = zeros(size(Δt_vec))
    L2_bdf2_fric = zeros(size(Δt_vec))
    L2_bdf3_fric = zeros(size(Δt_vec))

    # --------------------- Numerics --------------------- #

    for idt in eachindex(fac_vec)

        # ------------------- BDF 1 ------------------- #
        params = (
            K     = Kv,
            G     = Gv,
            η     = 1e10,
            c     = 1e10,
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
        CaseB_1D  = Main_VEP_1D_vdev_coss_BDF(σiB; params, visu=false, Ncy=Ncy)

        # Analytics
        t_ana   = CaseB_1D.t
        app_fric_ana = params.η*params.γ̇xy*(1 .- exp.(-params.G/params.η.*t_ana))./100e3 
        τxy_fric_ana = params.η*params.γ̇xy*(1 .- exp.(-params.G/params.η.*t_ana)) 

        L2_bdf1_τxy[idt]  = norm(τxy_fric_ana .- CaseB_1D.τxy) / norm(τxy_fric_ana)  
        L2_bdf1_fric[idt] = norm(app_fric_ana .- CaseB_1D.app_fric) / norm(app_fric_ana)  

        # ------------------- BDF 2 ------------------- #
        params = (
            K     = Kv,
            G     = Gv,
            η     = 1e10,
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
        CaseB_1D  = Main_VEP_1D_vdev_coss_BDF(σiB; params, visu=false, Ncy=Ncy)

        # Analytics
        t_ana   = CaseB_1D.t
        app_fric_ana = params.η*params.γ̇xy*(1 .- exp.(-params.G/params.η.*t_ana))./100e3 
        τxy_fric_ana = params.η*params.γ̇xy*(1 .- exp.(-params.G/params.η.*t_ana))

        L2_bdf2_τxy[idt]  = norm(τxy_fric_ana .- CaseB_1D.τxy) / norm(τxy_fric_ana)  
        L2_bdf2_fric[idt] = norm(app_fric_ana .- CaseB_1D.app_fric) / norm(app_fric_ana)
        
        # ------------------- BDF 3 ------------------- #
        params = (
            K     = Kv,
            G     = Gv,
            η     = 1e10,
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
        CaseB_1D  = Main_VEP_1D_vdev_coss_BDF(σiB; params, visu=false, Ncy=Ncy)

        # Analytics
        t_ana   = CaseB_1D.t
        app_fric_ana = params.η*params.γ̇xy*(1 .- exp.(-params.G/params.η.*t_ana))./100e3 
        τxy_fric_ana = params.η*params.γ̇xy*(1 .- exp.(-params.G/params.η.*t_ana)) 

        L2_bdf3_τxy[idt]  = norm(τxy_fric_ana .- CaseB_1D.τxy) / norm(τxy_fric_ana)  
        L2_bdf3_fric[idt] = norm(app_fric_ana .- CaseB_1D.app_fric) / norm(app_fric_ana)  

        # ------------------- VISU ------------------- #
  
        # Panel (1,1) - Stress ratio
        p1 = plot(title="Stress ratio", ylabel="-σxy/σyy")
        p1 = scatter!(CaseB_1D.γxy, CaseB_1D.app_fric, color=:green, label="numerics")
        # p1 = plot!(CaseB_LLP.γ_bulk[1:stp:end]*100, -CaseB_LLP.σ_in[ 3,1:stp:end] ./ CaseB_LLP.σ_in[ 2,1:stp:end], label="LLP",  color=:red)
        
        # p1 = scatter!(CaseB_1D.γxy, params.G.*CaseB_1D.γxy./100e5)
        p1 = scatter!(CaseB_1D.γxy, app_fric_ana)

        display(plot(p1))
    end

    @info "Convergence of τxy"
    @show L2_bdf1_τxy[1:end-1]./L2_bdf1_τxy[2:end-0]
    @show L2_bdf2_τxy[1:end-1]./L2_bdf2_τxy[2:end-0]
    @show L2_bdf3_τxy[1:end-1]./L2_bdf3_τxy[2:end-0] 

    @info "Convergence of τxy/σxy"
    @show L2_bdf1_fric[1:end-1]./L2_bdf1_fric[2:end-0]
    @show L2_bdf2_fric[1:end-1]./L2_bdf2_fric[2:end-0]
    @show L2_bdf3_fric[1:end-1]./L2_bdf3_fric[2:end-0] 

end

main()