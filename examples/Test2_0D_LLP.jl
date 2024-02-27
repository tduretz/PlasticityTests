using PlasticityTests, Plots

function Vermeer2_analytics_strain_LLP2013()
    degrad      = pi/180
    phi         = 40*degrad 
    psi         = 10*degrad 
    
    # Vermeer angle in LLP 2013
    # alpha_0     = atan(sin(phi)*cos(psi)/(1-sin(phi)*sin(psi))) 
    sigh        = -1.e5 
    G           = 1.e7
    nu          = 0.0 
    gammatot    = 0.1 
    nstep       = 201
    d_gxy       = gammatot/nstep; # strain increment
    alpha_0     = 39.999999*degrad
    Coulomb     = Vermeer2_analytics(phi,psi,nu,alpha_0,sigh,G,nstep,d_gxy)

    alpha_0     = 10. *degrad
    Roscoe      = Vermeer2_analytics(phi,psi,nu,alpha_0,sigh,G,nstep,d_gxy)

    alpha_0     = (phi+psi)/2 
    Arthur      = Vermeer2_analytics(phi,psi,nu,alpha_0,sigh,G,nstep,d_gxy)

    Test2 = ExtractDataTest2()
    p1 = plot()
    p1 = plot!(Arthur.gamma_xy,Arthur.sigma_v/sigh, label="Arthur")
    p1 = plot!(Coulomb.gamma_xy,Coulomb.sigma_v/sigh, label="Coulomb")
    p1 = plot!(Roscoe.gamma_xy,Roscoe.sigma_v/sigh, label="Roscoe")
    p1 = scatter!(Test2.StressRatioArthur.x./100,Test2.StressRatioArthur.y)
    p1 = scatter!(Test2.StressRatioCoulomb.x./100,Test2.StressRatioCoulomb.y,xlims=(-0.0,0.08),ylims=(1,5))

    p2 = plot()
    p2 = plot!(Arthur.gamma_xy, Arthur.sxx./Arthur.sxx0, label="Arthur")
    p2 = plot!(Coulomb.gamma_xy, Coulomb.sxx./Coulomb.sxx0, label="Coulomb")
    p2 = plot!(Roscoe.gamma_xy, Roscoe.sxx./Roscoe.sxx0, label="Roscoe",xlims=(-0.0,0.08),ylims=(0,1))
    
    p3 = plot()
    p3 = plot!(Arthur.gamma_xy, Arthur.theta, label="Arthur")
    p3 = plot!(Coulomb.gamma_xy, Coulomb.theta, label="Coulomb")
    p3 = plot!(Roscoe.gamma_xy, Roscoe.theta, label="Roscoe",xlims=(-0.0,0.08),ylims=(45,65))

    display(plot(p1,p2,p3, layout=(3,1)))

end

Vermeer2_analytics_strain_LLP2013()