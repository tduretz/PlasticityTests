using PlasticityTests, Plots

function main()

    # Load digitised data from test 1 of Vermeer (1990)
    CaseA = ExtractDataCase("CaseA")
    CaseB = ExtractDataCase("CaseB")

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
        nt  = 200,
        law = :DruckerPrager,
        el  = :Vermeer1990,
        pl  = true) # default parameter set

    # Case A
    σi       = (xx = -25e3, yy=-100e3)
    CaseA_0D = Vermeer1990_StressIntegration_vdev(σi)
    CaseA_1D = Main_VEP_1D(σi; visu=false)

    # Case B
    σi       = (xx = -400e3, yy=-100e3)
    CaseB_0D = Vermeer1990_StressIntegration_vdev(σi; params)
    CaseB_1D = Main_VEP_1D(σi; params, visu=false)

    #------------------------------#
    # Panel (1,1) - Stress ratio
    p1 = plot(title="Stress ratio", ylabel="-σxy/σyy", legend=:none)
    p1 = plot!(CaseA_0D.γxy, CaseA_0D.app_fric, label="0D", color=:blue)
    p1 = plot!(CaseA_1D.γxy, CaseA_1D.app_fric, label="1D", color=:green)    
    p1 = scatter!(CaseA.Friction.x, CaseA.Friction.y, label="V90", color=:orange)
    p1 = plot!(CaseB_0D.γxy, CaseB_0D.app_fric, color=:blue)
    p1 = plot!(CaseB_1D.γxy, CaseB_1D.app_fric, color=:green)
    p1 = scatter!(CaseB.Friction.x, CaseB.Friction.y, color=:orange)
    #------------------------------#
    # Panel (1,2) - Horizontal stress
    p2 = plot(title="Horizontal stress", ylabel="-σxx [kPa]")
    p2 = plot!(CaseA_0D.γxy,  CaseA_0D.σxx, label="0D", color=:blue)
    p2 = plot!(CaseA_1D.γxy,  CaseA_1D.σxx, label="1D", color=:green)
    p2 = scatter!(CaseA.σxx.x, CaseA.σxx.y, label="V90", color=:orange)
    p2 = plot!(CaseB_0D.γxy,  CaseB_0D.σxx, label=:none, color=:blue)
    p2 = plot!(CaseB_1D.γxy,  CaseB_1D.σxx, label=:none, color=:green)
    p2 = scatter!(CaseB.σxx.x, CaseB.σxx.y, label=:none, color=:orange)
    #------------------------------#
    # Panel (2,1) - Volume change
    p3 = plot(title="Volume change", ylabel="εyy [%]", xlabel="γxy [%]", legend=:none)
    p3 = plot!(CaseA_0D.γxy, CaseA_0D.εyy, label="0D", color=:blue)
    p3 = plot!(CaseA_1D.γxy,  CaseA_1D.εyy, label="1D", color=:green)
    p3 = scatter!(CaseA.εyy.x, CaseA.εyy.y, label="V90", color=:orange)
    p3 = plot!(CaseB_0D.γxy, CaseB_0D.εyy, color=:blue)
    p3 = plot!(CaseB_1D.γxy, CaseB_1D.εyy, color=:green)
    p3 = scatter!(CaseB.εyy.x, CaseB.εyy.y, color=:orange)
    # Panel (2,2) - Stress orientation
    p4 = plot(title="Stress orientation", ylabel="θ [ᵒ]", xlabel="γxy [%]", legend=:none)
    p4 = plot!(CaseA_0D.γxy, CaseA_0D.θ, label="0D", color=:blue)
    p4 = plot!(CaseA_1D.γxy, CaseA_1D.θ, label="1D", color=:green)
    p4 = scatter!(CaseA.θ.x, CaseA.θ.y, label="V90", color=:orange)
    p4 = plot!(CaseB_0D.γxy, CaseB_0D.θ, color=:blue)
    p4 = plot!(CaseB_1D.γxy, CaseB_1D.θ, color=:green)
    p4 = scatter!(CaseB.θ.x, CaseB.θ.y, color=:orange)
    #------------------------------# 
    display(plot(p1, p2, p3, p4, layout=(2,2)))

end

main()