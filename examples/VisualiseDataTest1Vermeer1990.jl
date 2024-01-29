using PlasticityTests, Plots

function main()

    # Load digitised data from test 1 of Vermeer (1990)
    CaseA = ExtractDataCase("CaseA")
    CaseB = ExtractDataCase("CaseB")

    #------------------------------#
    # Panel (1,1) - Stress ratio
    p1 = plot(title="Stress ratio", ylabel="-σxy/σyy")
    p1 = plot!(CaseA.Friction.x, CaseA.Friction.y, label="A")
    p1 = plot!(CaseB.Friction.x, CaseB.Friction.y, label="B")
    #------------------------------#
    # Panel (1,2) - Horizontal stress
    p2 = plot(title="Horizontal stress", ylabel="-σxx [kPa]")
    p2 = plot!(CaseA.Sxx.x, CaseA.Sxx.y, label="A")
    p2 = plot!(CaseB.Sxx.x, CaseB.Sxx.y, label="B")
    #------------------------------#
    # Panel (2,1) - Volume change
    p3 = plot(title="Volume change", ylabel="εyy [%]", xlabel="γxy [%]")
    p3 = plot!(CaseA.Eyy.x, CaseA.Eyy.y, label="A")
    p3 = plot!(CaseB.Eyy.x, CaseB.Eyy.y, label="B")
    # Panel (2,2) - Stress orientation
    p4 = plot(title="Stress orientation", ylabel="θ [ᵒ]", xlabel="γxy [%]")
    p4 = plot!(CaseA.θ.x, CaseA.θ.y, label="A")
    p4 = plot!(CaseB.θ.x, CaseB.θ.y, label="B")
    #------------------------------# 
    display(plot(p1, p2, p3, p4, layout=(2,2)))

end

main()