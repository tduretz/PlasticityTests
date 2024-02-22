using PlasticityTests, Plots

function main()

    # Load digitised data from test 1 of Vermeer (1990)
    CaseA = ExtractDataCase("CaseA")
    CaseB = ExtractDataCase("CaseB")

    # MC params
    ϕ  = 40.0*π/180.
    ψ  = 10.0*π/180.
    c  = 0.
    θt = 5.0*π/180

    # Vermeer (1990)
    Gv = 10e6
    Kv = 2/3*Gv

    # Case A - MC
    σi       = (xx = -25e3, yy=-100e3)
    CaseA_0D = Vermeer1990_StressIntegration_vdev(σi)

    # Case A - Drucker-Prager 
    params=(
        K   = Kv,
        G   = 10e6,
        c   = c,
        ϕ   = ϕ,
        ψ   = ψ,
        θt  = 25/180*π,
        ηvp = 0.,
        γ̇xy = 0.00001,
        Δt  = 20,
        nt  = 400,
        law = :DruckerPrager,
        oop = :Vermeer1990,
        pl  = true)
    CaseA_0D_DP0 = Vermeer1990_StressIntegration_vdev(σi; params)

    # # Case A - Drucker-Prager: compression fit - Circumbscribed
    # α = 2*sin(ϕ) / (sqrt(3)*(3 - sin(ϕ)))
    # β = 2*sin(ψ) / (sqrt(3)*(3 - sin(ψ)))
    # H = c/tan(ϕ)
    # params=(
    #     K   = Kv,
    #     G   = 10e6,
    #     c   = 3*α*H/cos(asin(3*α)),
    #     ϕ   = 3*α,
    #     ψ   = 3*β,
    #     θt  = 25/180*π,
    #     ηvp = 0.,
    #     γ̇xy = 0.00001,
    #     Δt  = 20,
    #     nt  = 400,
    #     law = :DruckerPrager,
    #     oop = :Vermeer1990,
    #     pl  = true)
    # CaseA_0D_DP1 = Vermeer1990_StressIntegration_vdev(σi; params)

    # # Case A - Drucker-Prager: extension fit - Middle circumscribed
    # α = 2*sin(ϕ) / (sqrt(3)*(3 + sin(ϕ)))
    # β = 2*sin(ψ) / (sqrt(3)*(3 + sin(ψ)))
    # H = c/tan(ϕ)
    # params=(
    #     K   = Kv,
    #     G   = 10e6,
    #     c   = 3*α*H/cos(asin(3*α)),
    #     ϕ   = 3*α,
    #     ψ   = 3*β,
    #     θt  = 25/180*π,
    #     ηvp = 0.,
    #     γ̇xy = 0.00001,
    #     Δt  = 20,
    #     nt  = 400,
    #     law = :DruckerPrager,
    #     oop = :Vermeer1990,
    #     pl  = true)
    # CaseA_0D_DP2 = Vermeer1990_StressIntegration_vdev(σi; params)

    # # Case A - Drucker-Prager: inscribed fit
    # α = 1*sin(ϕ) / (sqrt(9 + 3*sin(ϕ)^2))
    # β = 1*sin(ψ) / (sqrt(9 + 3*sin(ψ)^2))
    # H = c/tan(ϕ)
    # params=(
    #     K   = Kv,
    #     G   = 10e6,
    #     c   = 3*α*H/cos(asin(3*α)),
    #     ϕ   = 3*α,
    #     ψ   = 3*β,
    #     θt  = 25/180*π,
    #     ηvp = 0.,
    #     γ̇xy = 0.00001,
    #     Δt  = 20,
    #     nt  = 400,
    #     law = :DruckerPrager,
    #     oop = :Vermeer1990,
    #     pl  = true)
    # CaseA_0D_DP3 = Vermeer1990_StressIntegration_vdev(σi; params)

    # Case A - Mohr-Coulomb from Abbo and Sloan (1995) 
    params=(
        K   = Kv,
        G   = 10e6,
        c   = c,
        ϕ   = ϕ,
        ψ   = ψ,
        θt  = θt,
        ηvp = 0.,
        γ̇xy = 0.00001,
        Δt  = 20,
        nt  = 400,
        law = :MC_AS95,
        oop = :Vermeer1990,
        pl  = true)
    CaseA_0D_MC_AB = Vermeer1990_StressIntegration_vdev(σi; params)

    # Case A - Mohr-Coulomb from de Borst (1990) 
    params=(
        K   = Kv,
        G   = 10e6,
        c   = c,
        ϕ   = ϕ,
        ψ   = ψ,
        θt  = θt,
        ηvp = 0.,
        γ̇xy = 0.00001,
        Δt  = 20,
        nt  = 140,
        law = :MC_deBorst90,
        oop = :Vermeer1990,
        pl  = true)
    CaseA_0D_MC_dB = Vermeer1990_StressIntegration_vdev(σi; params)

    #-----------------------------------------------------------------#

    # Case B
    σi       = (xx = -400e3, yy=-100e3)
    CaseB_0D = Vermeer1990_StressIntegration_vdev(σi)

    # Case B - Drucker-Prager 
    params=(
        K   = Kv,
        G   = 10e6,
        c   = c,
        ϕ   = ϕ,
        ψ   = ψ,
        θt  = 25/180*π,
        ηvp = 0.,
        γ̇xy = 0.00001,
        Δt  = 20,
        nt  = 400,
        law = :DruckerPrager,
        oop = :Vermeer1990,
        pl  = true)
    CaseB_0D_DP0 = Vermeer1990_StressIntegration_vdev(σi; params)

    # # Case B - Drucker-Prager: compression fit
    # α = 2*sin(ϕ) / (sqrt(3)*(3 - sin(ϕ)))
    # β = 2*sin(ψ) / (sqrt(3)*(3 - sin(ψ)))
    # H = c/tan(ϕ)
    # params=(
    #     K   = Kv,
    #     G   = 10e6,
    #     c   = 3*α*H/cos(asin(3*α)),
    #     ϕ   = asin(3*α),
    #     ψ   = asin(3*β),
    #     θt  = 25/180*π,
    #     ηvp = 0.,
    #     γ̇xy = 0.00001,
    #     Δt  = 20,
    #     nt  = 400,
    #     law = :DruckerPrager,
    #     oop = :Vermeer1990,
    #     pl  = true)
    # CaseB_0D_DP1 = Vermeer1990_StressIntegration_vdev(σi; params)

    # # Case B - Drucker-Prager: extension fit
    # α = 2*sin(ϕ) / (sqrt(3)*(3 + sin(ϕ)))
    # β = 2*sin(ψ) / (sqrt(3)*(3 + sin(ψ)))
    # H = c/tan(ϕ)
    # params=(
    #     K   = Kv,
    #     G   = 10e6,
    #     c   = 3*α*H/cos(asin(3*α)),
    #     ϕ   = asin(3*α),
    #     ψ   = asin(3*β),
    #     θt  = 25/180*π,
    #     ηvp = 0.,
    #     γ̇xy = 0.00001,
    #     Δt  = 20,
    #     nt  = 400,
    #     law = :DruckerPrager,
    #     oop = :Vermeer1990,
    #     pl  = true)
    # CaseB_0D_DP2 = Vermeer1990_StressIntegration_vdev(σi; params)

    # # Case B - Drucker-Prager: inscribed fit
    # α = 1*sin(ϕ) / (sqrt(9 + 3*sin(ϕ)^2))
    # β = 1*sin(ψ) / (sqrt(9 + 3*sin(ψ)^2))
    # H = c/tan(ϕ)
    # params=(
    #     K   = Kv,
    #     G   = 10e6,
    #     c   = 3*α*H/cos(asin(3*α)),
    #     ϕ   = asin(3*α),
    #     ψ   = asin(3*β),
    #     θt  = 25/180*π,
    #     ηvp = 0.,
    #     γ̇xy = 0.00001,
    #     Δt  = 20,
    #     nt  = 400,
    #     law = :DruckerPrager,
    #     oop = :Vermeer1990,
    #     pl  = true)
    # CaseB_0D_DP3 = Vermeer1990_StressIntegration_vdev(σi; params)

    # Case B - Mohr-Coulomb from Abbo and Sloan (1995) 
    params=(
        K   = Kv,
        G   = 10e6,
        c   = c,
        ϕ   = ϕ,
        ψ   = ψ,
        θt  = θt,
        ηvp = 0.,
        γ̇xy = 0.00001,
        Δt  = 20,
        nt  = 400,
        law = :MC_AS95,
        oop = :Vermeer1990,
        pl  = true)
    CaseB_0D_MC_AB = Vermeer1990_StressIntegration_vdev(σi; params)

    # Case B - Mohr-Coulomb from de Borst (1990) 
    params=(
        K   = Kv,
        G   = 10e6,
        c   = c,
        ϕ   = ϕ,
        ψ   = ψ,
        θt  = θt,
        ηvp = 0.,
        γ̇xy = 0.00001,
        Δt  = 20,
        nt  = 10,
        law = :MC_deBorst90,
        oop = :Vermeer1990,
        pl  = true)
    CaseB_0D_MC_dB = Vermeer1990_StressIntegration_vdev(σi; params)

    # ------------------------------ #
    # Panel (1,1) - Stress ratio
    p1 = plot(title="Stress ratio", ylabel="-σxy/σyy", legend=:none, size=(300,300))
    # ---- #
    p1 = plot!(CaseA_0D.γxy, CaseA_0D.app_fric, label="A: 0D", color=:blue)
    p1 = plot!(CaseA_0D.γxy, CaseA_0D_DP0.app_fric, label="A: 0D DP0", color=:blue, linestyle=:dash)
    # p1 = plot!(CaseA_0D.γxy, CaseA_0D_DP1.app_fric, label="A: 0D DP1", color=:blue, linestyle=:dot)
    # p1 = plot!(CaseA_0D.γxy, CaseA_0D_DP2.app_fric, label="A: 0D DP2", color=:blue, linestyle=:dashdotdot)
    # p1 = plot!(CaseA_0D.γxy, CaseA_0D_DP3.app_fric, label="A: 0D DP3", color=:blue, linestyle=:dashdot)
    p1 = scatter!(CaseA_0D_MC_AB.γxy[1:10:end], CaseA_0D_MC_AB.app_fric[1:10:end], label="S: 0D MC AB95", color=:blue, marker=:cross)
    p1 = scatter!(CaseA_0D_MC_dB.γxy[5:10:end], CaseA_0D_MC_dB.app_fric[5:10:end], label="S: 0D MC dB90", color=:blue, marker=:star)
    p1 = scatter!(CaseA.Friction.x, CaseA.Friction.y, label="A: V90", color=:blue)
    # ---- #
    p1 = plot!(CaseB_0D.γxy, CaseB_0D.app_fric, label="B: 0D", color=:green)
    p1 = plot!(CaseB_0D.γxy, CaseB_0D_DP0.app_fric, label="B: 0D DP0", color=:green, linestyle=:dash)
    # p1 = plot!(CaseB_0D.γxy, CaseB_0D_DP1.app_fric, label="B: 0D DP1", color=:green, linestyle=:dot)
    # p1 = plot!(CaseB_0D.γxy, CaseB_0D_DP2.app_fric, label="B: 0D DP2", color=:green, linestyle=:dashdotdot)
    # p1 = plot!(CaseB_0D.γxy, CaseB_0D_DP3.app_fric, label="B: 0D DP3", color=:green, linestyle=:dashdot)
    p1 = scatter!(CaseB_0D_MC_AB.γxy[1:10:end], CaseB_0D_MC_AB.app_fric[1:10:end], label="B: 0D MC AB95", color=:green, marker=:cross)
    p1 = scatter!(CaseB_0D_MC_dB.γxy[5:10:end], CaseB_0D_MC_dB.app_fric[3:10:end], label="B: 0D MC dB90", color=:green, marker=:star)
    p1 = scatter!(CaseB.Friction.x, CaseB.Friction.y, label="B: V90", color=:green)
    # ------------------------------ #
    # Panel (1,2) - Horizontal stress
    p2 = plot(title="Horizontal stress", ylabel="-σxx [kPa]")
    # # ---- #
    p2 = plot!(CaseA_0D.γxy, CaseA_0D.σxx, label="A: 0D", color=:blue)
    p2 = plot!(CaseA_0D.γxy, CaseA_0D_DP0.σxx, label="A: 0D DP0", color=:blue, linestyle=:dash)
    # p2 = plot!(CaseA_0D.γxy, CaseA_0D_DP1.σxx, label="A: 0D DP1", color=:blue, linestyle=:dot)
    # p2 = plot!(CaseA_0D.γxy, CaseA_0D_DP2.σxx, label="A: 0D DP2", color=:blue, linestyle=:dashdotdot)
    # p2 = plot!(CaseA_0D.γxy, CaseA_0D_DP3.σxx, label="A: 0D DP3", color=:blue, linestyle=:dashdot)
    p2 = scatter!(CaseA_0D_MC_AB.γxy[1:10:end], CaseA_0D_MC_AB.σxx[1:10:end], label="A: 0D MC AB95", color=:blue, marker=:cross)
    p2 = scatter!(CaseA_0D_MC_dB.γxy[5:10:end], CaseA_0D_MC_dB.σxx[5:10:end], label="A: 0D MC dB90", color=:blue, marker=:star)    
    p2 = scatter!(CaseA.σxx.x, CaseA.σxx.y, label="A: V90", color=:blue)
    # ---- #
    p2 = plot!(CaseB_0D.γxy,  CaseB_0D.σxx, label="B: 0D", color=:green)
    p2 = plot!(CaseB_0D.γxy, CaseB_0D_DP0.σxx, label="B: 0D DP0", color=:green, linestyle=:dash)
    # p2 = plot!(CaseB_0D.γxy, CaseB_0D_DP1.σxx, label="B: 0D DP1", color=:green, linestyle=:dot)
    # p2 = plot!(CaseB_0D.γxy, CaseB_0D_DP2.σxx, label="B: 0D DP2", color=:green, linestyle=:dashdotdot)
    # p2 = plot!(CaseB_0D.γxy, CaseB_0D_DP3.σxx, label="B: 0D DP3", color=:green, linestyle=:dashdot)
    p2 = scatter!(CaseB_0D_MC_AB.γxy[1:10:end], CaseB_0D_MC_AB.σxx[1:10:end], label="B: 0D MC AB95", color=:green, marker=:cross)
    p2 = scatter!(CaseB_0D_MC_dB.γxy[5:10:end], CaseB_0D_MC_dB.σxx[3:10:end], label="B: 0D MC dB90", color=:green, marker=:star)
    p2 = scatter!(CaseB.σxx.x, CaseB.σxx.y, label="B: V90", color=:green)
    # ------------------------------ #
    # Panel (2,1) - Volume change
    p3 = plot(title="Volume change", ylabel="εyy [%]", xlabel="γxy [%]", legend=:none)
    # ---- #
    p3 = plot!(CaseA_0D.γxy, CaseA_0D.εyy, label="A: 0D", color=:blue)
    p3 = plot!(CaseA_0D.γxy, CaseA_0D_DP0.εyy, label="A: 0D DP0", color=:blue, linestyle=:dash)
    # p3 = plot!(CaseA_0D.γxy, CaseA_0D_DP1.εyy, label="A: 0D DP1", color=:blue, linestyle=:dashdot)
    # p3 = plot!(CaseA_0D.γxy, CaseA_0D_DP2.εyy, label="A: 0D DP2", color=:blue, linestyle=:dashdotdot)
    # p3 = plot!(CaseA_0D.γxy, CaseA_0D_DP3.εyy, label="A: 0D DP3", color=:blue, linestyle=:dashdot)
    p3 = scatter!(CaseA_0D_MC_AB.γxy[1:10:end], CaseA_0D_MC_AB.εyy[1:10:end], label="A: 0D MC AB95", color=:blue, marker=:cross)
    p3 = scatter!(CaseA_0D_MC_dB.γxy[5:10:end], CaseA_0D_MC_dB.εyy[3:10:end], label="A: 0D MC dB90", color=:blue, marker=:star)
    p3 = scatter!(CaseA.εyy.x, CaseA.εyy.y, label="A: V90", color=:blue)
    # ---- #
    p3 = plot!(CaseB_0D.γxy, CaseB_0D.εyy, label="B: 0D", color=:green)
    p3 = plot!(CaseB_0D.γxy, CaseB_0D_DP0.εyy, label="B: 0D DP0", color=:green, linestyle=:dash)
    # p3 = plot!(CaseB_0D.γxy, CaseB_0D_DP1.εyy, label="B: 0D DP1", color=:green, linestyle=:dot)
    # p3 = plot!(CaseB_0D.γxy, CaseB_0D_DP2.εyy, label="B: 0D DP2", color=:green, linestyle=:dashdotdot)
    # p3 = plot!(CaseB_0D.γxy, CaseB_0D_DP3.εyy, label="B: 0D DP3", color=:green, linestyle=:dashdot)
    p3 = scatter!(CaseB_0D_MC_AB.γxy[1:10:end], CaseB_0D_MC_AB.εyy[1:10:end], label="B: 0D MC AB95", color=:green, marker=:cross)
    p3 = scatter!(CaseB_0D_MC_dB.γxy[5:10:end], CaseB_0D_MC_dB.εyy[3:10:end], label="B: 0D MC dB90", color=:green, marker=:star)    
    p3 = scatter!(CaseB.εyy.x, CaseB.εyy.y, label="B: V90", color=:green)
    # ------------------------------ # 
    # Panel (2,2) - Stress orientation
    p4 = plot(title="Stress orientation", ylabel="θ [ᵒ]", xlabel="γxy [%]", legend=:none)
    # ---- #
    p4 = plot!(CaseA_0D.γxy, CaseA_0D.θ, label="A: 0D", color=:blue)
    p4 = plot!(CaseA_0D.γxy, CaseA_0D_DP0.θ, label="A: 0D DP0", color=:blue, linestyle=:dash)
    # p4 = plot!(CaseA_0D.γxy, CaseA_0D_DP1.θ, label="A: 0D DP1", color=:blue, linestyle=:dot)
    # p4 = plot!(CaseA_0D.γxy, CaseA_0D_DP2.θ, label="A: 0D DP2", color=:blue, linestyle=:dashdotdot)
    # p4 = plot!(CaseA_0D.γxy, CaseA_0D_DP3.θ, label="A: 0D DP3", color=:blue, linestyle=:dashdot)
    p4 = scatter!(CaseA_0D_MC_AB.γxy[1:10:end], CaseA_0D_MC_AB.θ[1:10:end], label="A: 0D MC AB95", color=:blue, marker=:cross)
    p4 = scatter!(CaseA_0D_MC_dB.γxy[5:10:end], CaseA_0D_MC_dB.θ[3:10:end], label="A: 0D MC dB90", color=:blue, marker=:star)
    p4 = scatter!(CaseA.θ.x, CaseA.θ.y, label="A: V90", color=:blue)
    # ---- #
    p4 = plot!(CaseB_0D.γxy, CaseB_0D.θ, label="B: 0D", color=:green)
    p4 = plot!(CaseB_0D.γxy, CaseB_0D_DP0.θ, label="B: 0D DP0", color=:green, linestyle=:dash)
    # p4 = plot!(CaseB_0D.γxy, CaseB_0D_DP1.θ, label="B: 0D DP1", color=:green, linestyle=:dot)
    # p4 = plot!(CaseB_0D.γxy, CaseB_0D_DP2.θ, label="B: 0D DP2", color=:green, linestyle=:dot) 
    # p4 = plot!(CaseB_0D.γxy, CaseB_0D_DP3.θ, label="B: 0D DP3", color=:green, linestyle=:dashdot)
    p4 = scatter!(CaseB_0D_MC_AB.γxy[1:10:end], CaseB_0D_MC_AB.θ[1:10:end], label="B: 0D MC AB95", color=:green, marker=:cross)
    p4 = scatter!(CaseB_0D_MC_dB.γxy[5:10:end], CaseB_0D_MC_dB.θ[3:10:end], label="B: 0D MC dB90", color=:green, marker=:star)    
    p4 = scatter!(CaseB.θ.x, CaseB.θ.y, label="B: V90", color=:green)
    # ------------------------------ # 
    display(plot(p1, p2, p3, p4, layout=(2,2)))


end

main()