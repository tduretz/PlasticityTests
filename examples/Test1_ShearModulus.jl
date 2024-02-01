using PlasticityTests, Plots

function main()

    # Load digitised data from test 1 of Vermeer (1990)
    CaseB = ExtractDataCase("CaseB")

    # Case B - Gref/3.5
    params=(
        K   = 0.,
        G   = 2.857142857142857e6,
        c   = 0.0,
        ϕ   = 40/180*π,
        ψ   = 10/180*π,
        θt  = 25/180*π,
        ηvp = 0.,
        γ̇xy = 0.00001,
        Δt  = 20,
        nt  = 600,
        law = :DruckerPrager,
        el  = :Vermeer1990,
        pl  = true) # default parameter set
    σi       = (xx = -400e3, yy=-100e3)
    CaseB_0D_00 = Vermeer1990_StressIntegration_vdev(σi; params)
    CaseB_1D_00 = Main_VEP_1D(σi; params, visu=false)

    # Case B - Gref
    params=(
        K   = 0.,
        G   = 10e6,
        c   = 0.0,
        ϕ   = 40/180*π,
        ψ   = 10/180*π,
        θt  = 25/180*π,
        ηvp = 0.,
        γ̇xy = 0.00001,
        Δt  = 20,
        nt  = 600,
        law = :DruckerPrager,
        el  = :Vermeer1990,
        pl  = true) # default parameter set
    σi       = (xx = -400e3, yy=-100e3)
    CaseB_0D_10 = Vermeer1990_StressIntegration_vdev(σi; params)
    CaseB_1D_10 = Main_VEP_1D(σi; params, visu=false)

    # Case B - Gref*3.5
    params=(
        K   = 0.,
        G   = 35e6,
        c   = 0.0,
        ϕ   = 40/180*π,
        ψ   = 10/180*π,
        θt  = 25/180*π,
        ηvp = 0.,
        γ̇xy = 0.00001,
        Δt  = 20,
        nt  = 600,
        law = :DruckerPrager,
        el  = :Vermeer1990,
        pl  = true) # default parameter set
    σi       = (xx = -400e3, yy=-100e3)
    CaseB_0D_20 = Vermeer1990_StressIntegration_vdev(σi; params)
    CaseB_1D_20 = Main_VEP_1D(σi; params, visu=false)

    #------------------------------#
    # Panel (1,1) - Stress ratio
    p1 = plot(title="Stress ratio", ylabel="-σxy/σyy", legend=:none)
    _,imax = findmax(CaseB_0D_00.app_fric)
    p1 = plot!(CaseB_0D_00.γxy, CaseB_0D_00.app_fric, color=:blue)
    p1 = plot!(CaseB_1D_00.γxy, CaseB_1D_00.app_fric, color=:blue, linestyle=:dashdot)
    p1 = scatter!([CaseB_0D_00.γxy[imax]], [CaseB_0D_00.app_fric[imax]], marker=:star, color=:blue)

    _,imax = findmax(CaseB_0D_10.app_fric)
    p1 = plot!(CaseB_0D_10.γxy, CaseB_0D_10.app_fric, color=:green)
    p1 = plot!(CaseB_1D_10.γxy, CaseB_1D_10.app_fric, color=:green, linestyle=:dashdot)
    p1 = scatter!([CaseB_0D_10.γxy[imax]], [CaseB_0D_10.app_fric[imax]], marker=:star, color=:green)
    
    _,imax = findmax(CaseB_0D_20.app_fric)
    p1 = plot!(CaseB_0D_20.γxy, CaseB_0D_20.app_fric, color=:orange)
    p1 = plot!(CaseB_1D_20.γxy, CaseB_1D_20.app_fric, color=:orange, linestyle=:dashdot)
    p1 = scatter!([CaseB_0D_20.γxy[imax]], [CaseB_0D_20.app_fric[imax]], marker=:star, color=:orange)
    #------------------------------#
    # Panel (1,2) - Horizontal stress
    p2 = plot(title="Horizontal stress", ylabel="-σxx [kPa]")
    p2 = plot!(CaseB_0D_00.γxy,  CaseB_0D_00.σxx, label=:none, color=:blue)
    p2 = plot!(CaseB_1D_00.γxy,  CaseB_1D_00.σxx, label=:none, color=:blue, linestyle=:dashdot)
    p2 = plot!(CaseB_0D_10.γxy,  CaseB_0D_10.σxx, label=:none, color=:green)
    p2 = plot!(CaseB_1D_10.γxy,  CaseB_1D_10.σxx, label=:none, color=:green, linestyle=:dashdot)
    p2 = plot!(CaseB_0D_20.γxy,  CaseB_0D_20.σxx, label=:none, color=:orange)
    p2 = plot!(CaseB_1D_20.γxy,  CaseB_1D_20.σxx, label=:none, color=:orange, linestyle=:dashdot)
    #------------------------------#
    # Panel (2,1) - Volume change
    p3 = plot(title="Volume change", ylabel="εyy [%]", xlabel="γxy [%]", foreground_color_legend = nothing, background_color_legend = nothing)
    p3 = plot!(CaseB_0D_00.γxy, CaseB_0D_00.εyy, color=:blue, label="G = 2.85 MPa")
    p3 = plot!(CaseB_1D_00.γxy, CaseB_1D_00.εyy, color=:blue, label=:none, linestyle=:dashdot)
    p3 = plot!(CaseB_0D_10.γxy, CaseB_0D_10.εyy, color=:green, label="G = 10 MPa")
    p3 = plot!(CaseB_1D_10.γxy, CaseB_1D_10.εyy, color=:green, label=:none, linestyle=:dashdot)
    p3 = plot!(CaseB_0D_20.γxy, CaseB_0D_20.εyy, color=:orange, label="G = 35 MPa")
    p3 = plot!(CaseB_1D_20.γxy, CaseB_1D_20.εyy, color=:orange, label=:none, linestyle=:dashdot)
    # p3 = scatter!(CaseB.εyy.x, CaseB.εyy.y, color=:orange)
    # Panel (2,2) - Stress orientation
    p4 = plot(title="Stress orientation", ylabel="θ [ᵒ]", xlabel="γxy [%]", legend=:none)
    p4 = plot!(CaseB_0D_00.γxy, CaseB_0D_00.θ, color=:blue)
    p4 = plot!(CaseB_1D_00.γxy, CaseB_1D_00.θ, color=:blue, linestyle=:dashdot)
    p4 = plot!(CaseB_0D_10.γxy, CaseB_0D_10.θ, color=:green)
    p4 = plot!(CaseB_1D_10.γxy, CaseB_1D_10.θ, color=:green, linestyle=:dashdot)
    p4 = plot!(CaseB_0D_20.γxy, CaseB_0D_20.θ, color=:orange)
    p4 = plot!(CaseB_1D_20.γxy, CaseB_1D_20.θ, color=:orange, linestyle=:dashdot)
    # p4 = scatter!(CaseB.θ.x, CaseB.θ.y, color=:orange)
    #------------------------------# 
    display(plot(p1, p2, p3, p4, layout=(2,2)))

end

main()