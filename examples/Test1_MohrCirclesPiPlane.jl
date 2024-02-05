using PlasticityTests, Plots
import CairoMakie
import Makie.SpecApi as S

function main()

    # Load digitised data from test 1 of Vermeer (1990)
    CaseA = ExtractDataCase("CaseA")
    CaseB = ExtractDataCase("CaseB")

    # MC params
    ϕ  = 40.0*π/180.
    ψ  = 10.0*π/180.
    c  = 0.
    θt = 2.0*π/180

    # Vermeer (1990)
    Gv = 10e6
    Kv = 2/3*Gv

    # Case A - MC
    σi       = (xx = -25e3, yy=-100e3)
    CaseA_0D = Vermeer1990_StressIntegration_vdev(σi)

    # Case A - Drucker-Prager 
    params=(
        K   = Kv,
        G   = Gv,
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

    #-----------------------------------------------------------------#

    # Case B
    σi       = (xx = -400e3, yy=-100e3)
    CaseB_0D = Vermeer1990_StressIntegration_vdev(σi)

    # Case B - Drucker-Prager 
    params=(
        K   = Kv,
        G   = Gv,
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

    #-----------------------------------------------------------------#
    θ   = LinRange(0, π, 100)
    σMC = LinRange(-500, 0, 100 ) .*1e3
    τMC = -σMC.*tan(ϕ) .+ c
    θl  = -((1/6)*π):0.1:((1/6)*π)

    J2A = zero(θl)
    J2B = zero(θl)
        
    τA  = abs.(CaseA_0D_DP0.σ1 - CaseA_0D_DP0.σ3)/2
    PA  = (CaseA_0D_DP0.σ1 + CaseA_0D_DP0.σ3)/2
    τB  = abs.(CaseB_0D_DP0.σ1 - CaseB_0D_DP0.σ3)/2
    PB  = (CaseB_0D_DP0.σ1 + CaseB_0D_DP0.σ3)/2
    θlA = CaseA_0D_DP0.θl
    θlB = CaseB_0D_DP0.θl


    # fig, _, spec = plot(S.GridLayout()) # create empty spec plot 

    # anim = @animate for it=2:params.nt
    # record(fig, "anim_MWE.gif", 1:10) do it
    for it=400:400#params.nt
       
        yield = (x = σMC./1e3, y = τMC./1e3)
        MC_A  = (x = (τA[it].*cos.(θ).+PA[it])./1e3, y = (τA[it].*sin.(θ))./1e3) 
        MC_B  = (x = (τB[it].*cos.(θ).+PB[it])./1e3, y = (τB[it].*sin.(θ))./1e3)

        # p1 = plot(title="Test 1 Mohr circles", ylabel="τ", xlabel="σₙ", size=(300,300), aspect_ratio=1)
        # p1 = plot!( MC_A... , color=:blue, label="Case A" )
        # p1 = plot!( MC_B...,  color=:green, label="Case B"  )
        # p1 = plot!( yield..., color=:red, label="Yield"  )

        fig = CairoMakie.Figure(size = (600, 600))

        # spec[1] = S.GridLayout(S.Axis(title="Animation", plots=[
        #     S.Lines(1:10, (1:10).^2),
        #     S.Scatter([it], [it^2], )
        # ]))

        ax1 = CairoMakie.Axis(fig[1, 1], title = "Mohr circles", xlabel="σₙ [kPa]", ylabel="τ [kPa]")
        ax1.aspect = 500/400
        CairoMakie.lines!(ax1, MC_A..., label="Case A")
        CairoMakie.lines!(ax1, MC_B..., label="Case B")
        CairoMakie.lines!(ax1, yield..., color=:red, label="Yield")
        CairoMakie.ylims!(0, 400)
        CairoMakie.xlims!(-500, 0)
        CairoMakie.axislegend(position = :rt)

        # Plot in th π-plane. 
        # Radial coordinate:       sqrt.(2 .*J2) 
        # Longitudinal coordinate: θl

        # See notebook
        @. J2A = 9 * (PA[it] .* sin(ϕ) + c .* cos(ϕ)) .^ 2 ./ (sqrt(3) * sin(θl) .* sin(ϕ) - 3 * cos(θl)) .^ 2
        @. J2B = 9 * (PB[it] .* sin(ϕ) + c .* cos(ϕ)) .^ 2 ./ (sqrt(3) * sin(θl) .* sin(ϕ) - 3 * cos(θl)) .^ 2

        ax2 = CairoMakie.PolarAxis(fig[1:1, 2], title = "π-plane", thetalimits = (-pi/6, pi/6))
        CairoMakie.lines!(ax2, θl, sqrt.(2 .*J2A)./1e3)
        CairoMakie.scatter!(ax2, θlA[it], sqrt(2).*τA[it]./1e3)
        CairoMakie.lines!(ax2, θl, sqrt.(2 .*J2B)./1e3)
        CairoMakie.scatter!(ax2, θlB[it], sqrt(2).*τB[it]./1e3)
    
        display(fig)
        # sleep(0.001)
    end
    # gif(anim, "figures/Test1_MohrCircles_v2.gif", fps = 15)

end

main()