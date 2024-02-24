using PlasticityTests, Plots

function main()

    # Vermeer (1990)
    Gv  = 10e6
    Kv  = 2/3*Gv
    σiA = (xx =  -25e3, yy=-100e3, xy=0.0) # Case A
    σiB = (xx = -400e3, yy=-100e3, xy=0.0) # Case B
     
    # Resolution
    Ncy     = 10
    dt_fact = 1  

    # 1D model
    params = (
        K    = Kv,
        G    = Gv,
        c    = 0.0,
        ϕ    = 40/180*π,
        ψ    = 10/180*π,
        θt   = 25/180*π,
        ηvp  = 0e7,
        lc   = 1e2,
        γ̇xy  = 0.00001,
        Δt   = 20/dt_fact,
        nt   = 800*dt_fact,
        law  = :DruckerPrager,
        coss = false,
        oop  = :Vermeer1990,
        pl   = true) # default parameter set
    
    # CaseA_1D  = Main_VEP_1D_vdev_coss(σiA; params, visu=false, Ncy=Ncy)
    CaseB_1D  = Main_VEP_1D_vdev(σiB; params, visu=true, Ncy=Ncy)

    # LLP solution
    params = (
        K    = Kv,
        G    = Gv,
        c    = 0.0,
        ϕ    = 40/180*π,
        ψ    = 10/180*π,
        θt   = 25/180*π,
        ηvp  = 0.,
        lc   = 1e3,
        γ̇xy  = 0.00001,
        Δt   = 2,
        nt   = 4000,
        law  = :DruckerPrager,
        coss = true,
        oop  = :Vermeer1990,
        pl   = true) # default parameter set

    CaseA_LLP = Vermeer3_ana_llp2013(σiA, params, Ncy)
    CaseB_LLP = Vermeer3_ana_llp2013(σiB, params, Ncy)

    #------------------------------#
    stp = 1 # scatter step
    # Panel (1,1) - Stress ratio
    p1 = plot(title="Stress ratio", ylabel="-σxy/σyy", legend=:none)
    # p1 = plot!(CaseA_1D.γxy, CaseA_1D.app_fric, label="1D", color=:green)    
    # p1 = plot!(CaseA_LLP.ε_out[3,1:stp:end]*100, -CaseA_LLP.σ_in[ 3,1:stp:end] ./ CaseA_LLP.σ_in[ 2,1:stp:end], label="sxy/syy_in",  color=:red)
    # p1 = plot!(CaseA_LLP.ε_out[3,1:stp:end]*100, -CaseA_LLP.σ_out[3,1:stp:end] ./ CaseA_LLP.σ_out[2,1:stp:end], label="sxy/syy_out", color=:blue)
    
    p1 = plot!(CaseB_1D.γxy, CaseB_1D.app_fric, color=:green)
    p1 = plot!(CaseB_LLP.γ_bulk[1:stp:end]*100, -CaseB_LLP.σ_in[ 3,1:stp:end] ./ CaseB_LLP.σ_in[ 2,1:stp:end], label="sxy/syy_in",  color=:red)
    p1 = plot!(CaseB_LLP.γ_bulk[1:stp:end]*100, -CaseB_LLP.σ_out[3,1:stp:end] ./ CaseB_LLP.σ_out[2,1:stp:end], label="sxy/syy_out", color=:red)

    #------------------------------#
    # Panel (1,2) - Horizontal stress
    p2 = plot(title="Horizontal stress", ylabel="-σxx [kPa]", legend=:none)
    # p2 = plot!(CaseA_1D.γxy,  -CaseA_1D.σxx_in, label="1D in", color=:green)
    # p2 = plot!(CaseA_LLP.ε_out[3,1:stp:end]*100, -CaseA_LLP.σ_in[1, 1:stp:end]./1e3, label=:none, color=:blue)
    # p2 = plot!(CaseA_LLP.ε_out[3,1:stp:end]*100, -CaseA_LLP.σ_out[1, 1:stp:end]./1e3, label=:none, color=:blue)

    p2 = plot!(CaseB_1D.γxy,  -CaseB_1D.σxx_in, color=:green)
    p2 = plot!(CaseB_1D.γxy,  -CaseB_1D.σxx_out, label=:none, color=:green)
    p2 = plot!(CaseB_LLP.γ_bulk[1:stp:end]*100, -CaseB_LLP.σ_in[ 1, 1:stp:end]./1e3, color=:red)
    p2 = plot!(CaseB_LLP.γ_bulk[1:stp:end]*100, -CaseB_LLP.σ_out[1, 1:stp:end]./1e3, color=:red)

    #------------------------------#
    # Panel (2,1) - Volume change
    p3 = plot(title="Volume change", ylabel="εyy [%]", xlabel="γxy [%]")
    # p3 = plot!(CaseA_1D.γxy,  CaseA_1D.εyy, label="1D", color=:green)
    # p3 = plot!(CaseA_LLP.ε_out[3,1:stp:end]*100, CaseA_LLP.ε_in[2, 1:stp:end]*100, label=:none, color=:blue)

    p3 = plot!(CaseB_1D.γxy, CaseB_1D.εyy_in, color=:green, label="1D in")
    p3 = plot!(CaseB_1D.γxy, CaseB_1D.εyy_out, color=:green, label="1D out")

    p3 = plot!(CaseB_LLP.γ_bulk[1:stp:end]*100, CaseB_LLP.ε_in[ 2, 1:stp:end]*100, label="LLP in" , color=:red)
    p3 = plot!(CaseB_LLP.γ_bulk[1:stp:end]*100, CaseB_LLP.ε_out[2, 1:stp:end]*100, label="LLP out", color=:red)


    # Panel (2,2) - Stress orientation
    p4 = plot(title="Stress orientation", ylabel="θ [ᵒ]", xlabel="γxy [%]", legend=:none)
    # p4 = plot!(CaseA_1D.γxy, CaseA_1D.θ, label="1D", color=:green)
    # p4 = plot!(CaseA_LLP.ε_out[3,1:stp:end]*100, CaseA_LLP.θ_in[1:stp:end], label=:none, color=:blue)
    
    p4 = plot!(CaseB_1D.γxy, CaseB_1D.θ_in, color=:green)
    p4 = plot!(CaseB_1D.γxy, CaseB_1D.θ_out, color=:green)
    p4 = plot!(CaseB_LLP.γ_bulk[1:stp:end]*100, CaseB_LLP.θ_in[1:stp:end], label=:none, color=:red)
    p4 = plot!(CaseB_LLP.γ_bulk[1:stp:end]*100, CaseB_LLP.θ_out[1:stp:end], label=:none, color=:red)


    #------------------------------# 
    display(plot(p1, p2, p3, p4, layout=(2,2)))
    # display(plot(p2))

end

main()