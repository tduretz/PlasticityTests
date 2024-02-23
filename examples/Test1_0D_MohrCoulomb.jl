using PlasticityTests, Plots

function main()

    # Load digitised data from test 1 of Vermeer (1990)
    CaseA = ExtractDataCase("CaseA")
    CaseB = ExtractDataCase("CaseB")

    # Formulation of stress-strain rate relationship (:voldev or :total)
    formulation = :voldev 
    # formulation = :total

    # Case A
    σi       = (xx = -25e3, yy=-100e3, xy=0.0)
    if formulation==:voldev
        CaseA_0D = Vermeer1990_StressIntegration_vdev(σi)
    elseif formulation==:total
        CaseA_0D = Vermeer1990_StressIntegration_tot(σi)
    end

    # Case B
    σi       = (xx = -400e3, yy=-100e3, xy=0.0)
    if formulation==:voldev
        CaseB_0D = Vermeer1990_StressIntegration_vdev(σi)
    elseif formulation==:total
        CaseB_0D = Vermeer1990_StressIntegration_tot(σi)
    end

     #------------------------------#
    # Panel (1,1) - Stress ratio
    p1 = plot(title="Stress ratio", ylabel="-σxy/σyy", legend=:none)
    p1 = plot!(CaseA_0D.γxy, CaseA_0D.app_fric, label="A: 0D", color=:blue)
    p1 = scatter!(CaseA.Friction.x, CaseA.Friction.y, label="A: V90", color=:blue)
    p1 = plot!(CaseB_0D.γxy, CaseB_0D.app_fric, label="B: 0D", color=:green)
    p1 = scatter!(CaseB.Friction.x, CaseB.Friction.y, label="B: V90", color=:green)
    #------------------------------#
    # Panel (1,2) - Horizontal stress
    p2 = plot(title="Horizontal stress", ylabel="-σxx [kPa]")
    p2 = plot!(CaseA_0D.γxy, CaseA_0D.σxx, label="A: 0D", color=:blue)
    p2 = scatter!(CaseA.σxx.x, CaseA.σxx.y, label="A: V90", color=:blue)
    p2 = plot!(CaseB_0D.γxy,  CaseB_0D.σxx, label="B: 0D", color=:green)
    p2 = scatter!(CaseB.σxx.x, CaseB.σxx.y, label="B: V90", color=:green)
    #------------------------------#
    # Panel (2,1) - Volume change
    p3 = plot(title="Volume change", ylabel="εyy [%]", xlabel="γxy [%]", legend=:none)
    p3 = plot!(CaseA_0D.γxy, CaseA_0D.εyy, label="A: 0D", color=:blue)
    p3 = scatter!(CaseA.εyy.x, CaseA.εyy.y, label="A: V90", color=:blue)
    p3 = plot!(CaseB_0D.γxy, CaseB_0D.εyy, label="B: 0D", color=:green)
    p3 = scatter!(CaseB.εyy.x, CaseB.εyy.y, label="B: V90", color=:green)
    # Panel (2,2) - Stress orientation
    p4 = plot(title="Stress orientation", ylabel="θ [ᵒ]", xlabel="γxy [%]", legend=:none)
    p4 = plot!(CaseA_0D.γxy, CaseA_0D.θ, label="A: 0D", color=:blue)
    p4 = scatter!(CaseA.θ.x, CaseA.θ.y, label="A: V90", color=:blue)
    p4 = plot!(CaseB_0D.γxy, CaseB_0D.θ, label="B: 0D", color=:green)
    p4 = scatter!(CaseB.θ.x, CaseB.θ.y, label="B: V90", color=:green)
    #------------------------------# 
    display(plot(p1, p2, p3, p4, layout=(2,2)))


end

main()