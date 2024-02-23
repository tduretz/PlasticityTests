using PlasticityTests,Plots 

function Vermeer3_ana()

    # Vermeer (1990)
    Gv = 10e6
    Kv = 2/3*Gv

    params=(
        K   = Kv,
        G   = Gv,
        c   = 0.0,
        ϕ   = 40/180*π,
        ψ   = 10/180*π,
        θt  = 25/180*π,
        ηvp = 0.,
        γ̇xy = 0.00001,
        Δt  = 20,
        nt  = 400,
        law = :DruckerPrager,
        oop = :Vermeer1990,
        pl  = true) # default parameter set

    σi       = (xx = -400e3, yy=-100e3, xy=0.0)    
    shearband_thickness = 1.
    model_thickness     = 10.
    factor = model_thickness/shearband_thickness
    LLP    = Vermeer3_ana_llp2013(σi, params, factor)
  
    # Plotting 
    p1=plot( LLP.ε_out[3, :], -LLP.σ_in[3, :] ./ LLP.σ_in[2, :], label="sxy/syy_in", color=:red)
    p1=plot!(LLP.ε_out[3, :], -LLP.σ_out[3, :] ./ LLP.σ_out[2, :], label="sxy/syy_out", color=:blue)

    p2=plot( LLP.ε_out[3, :], LLP.θ_in, label="θ_in", color=:red)
    p2=plot!(LLP.ε_out[3, :], LLP.θ_out, label="θ_out", color=:blue)

    p3=plot( LLP.ε_out[3, :], LLP.ε_in[2, :], label="eyy_in", color=:red)
    p3=plot!(LLP.ε_out[3, :], LLP.ε_out[2, :], label="eyy_out", color=:blue)

    p4=plot( LLP.ε_out[3, :], -LLP.σ_in[1, :], label="sxx_in", color=:red)
    p4=plot!(LLP.ε_out[3, :], -LLP.σ_out[1, :], label="sxx_out", color=:blue)

    display(plot(p1,p2,p3,p4))
end



Vermeer3_ana()