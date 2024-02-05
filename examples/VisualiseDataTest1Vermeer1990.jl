using PlasticityTests, Plots

function main()

    # Load digitised data from test 1 of Vermeer (1990)
    CaseA = ExtractDataCase("CaseA")
    CaseB = ExtractDataCase("CaseB")

    # Integration case A
    σi = (xx = -25e3, yy=-100e3)
    CaseA_Components = Vermeer1990_Test1_Components(σi)
    CaseA_MatVec     = Vermeer1990_Test1_MatVec(σi)

    σi = (xx = -400e3, yy=-100e3)
    CaseB_Components = Vermeer1990_Test1_Components(σi)
    CaseB_MatVec     = Vermeer1990_Test1_MatVec(σi)

    #------------------------------#
    # Panel (1,1) - Stress ratio
    p1 = plot(title="Stress ratio", ylabel="-σxy/σyy")
    p1 = scatter!(CaseA.Friction.x, CaseA.Friction.y, label="A: Paper scan")
    p1 = scatter!(CaseB.Friction.x, CaseB.Friction.y, label="B: Paper scan")

    p1 = plot!(CaseA_Components.γxy, CaseA_Components.app_fric, label="Components")
    p1 = plot!(CaseA_Components.γxy, CaseB_Components.app_fric, label=:none)

    p1 = plot!(CaseA_MatVec.γxy, CaseA_MatVec.app_fric, label="Matrix-Vector")
    p1 = plot!(CaseB_MatVec.γxy, CaseB_MatVec.app_fric, label=:none)

    #------------------------------#
    # Panel (1,2) - Horizontal stress
    p2 = plot(title="Horizontal stress", ylabel="-σxx [kPa]")
    p2 = scatter!(CaseA.σxx.x, CaseA.σxx.y, label="A")
    p2 = scatter!(CaseB.σxx.x, CaseB.σxx.y, label="B")
    #------------------------------#
    # Panel (2,1) - Volume change
    p3 = plot(title="Volume change", ylabel="εyy [%]", xlabel="γxy [%]")
    p3 = scatter!(CaseA.εyy.x, CaseA.εyy.y, label="A")
    p3 = scatter!(CaseB.εyy.x, CaseB.εyy.y, label="B")
    # Panel (2,2) - Stress orientation
    p4 = plot(title="Stress orientation", ylabel="θ [ᵒ]", xlabel="γxy [%]")
    p4 = scatter!(CaseA.θ.x, CaseA.θ.y, label="A")
    p4 = scatter!(CaseB.θ.x, CaseB.θ.y, label="B")
    #------------------------------# 
    display(plot(p1, p2, p3, p4, layout=(2,2)))

end

main()