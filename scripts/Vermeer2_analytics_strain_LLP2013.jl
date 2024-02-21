using PlasticityTests,Plots
function Vermeer2_analytics_strain_LLP2013()
    degrad= pi/180
    phi         = 40*degrad 
    psi         = 10*degrad 
    alpha_0     = 39.999999*degrad
    
    # Vermeer angle in LLP 2013
    # alpha_0     = atan(sin(phi)*cos(psi)/(1-sin(phi)*sin(psi))) 
    sigh        = -1.e5 
    G           = 1.e7
    nu          = 0.0 
    gammatot    = 0.1 
    nstep       = 201
    d_gxy       = gammatot/nstep; # strain increment
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

function  Vermeer2_analytics(phi,psi,nu,alpha0,sigh,G,nstep,d_gxy)
    sigma_v  = zeros(nstep) # vertical stress outside the band
    beta     = zeros(nstep) # orientation of stress in the shear band
    gamma_xy = zeros(nstep) # engineering shear strain in the band
    sxx      = zeros(nstep) # band-parallel normal stress in
    sxx0     = zeros(nstep) # band-parallel normal stress out
    syy      = zeros(nstep) # band-normal  normal stress 
    sxy      = zeros(nstep) # band-parallel shear stress 
    #initial values
    S_ini       = (1+sin(phi))/(1-sin(phi));
    sigma_v[1]   = S_ini*sigh;
    beta[1]      = alpha0;  
    gamma_xy[1]  = 0;
    R            = sigma_v[1]-sigh;
    # 4 independant stress components 
    sxy[1]       = -R*cos(alpha0)/2;
    syy[1]       =  (sigma_v[1]+sigh-R*sin(alpha0))/2;
    sxx0[1]      =  (sigma_v[1]+sigh+R*sin(alpha0))/2;
    sxx[1]       =  sxx0[1];

    for n=1:nstep-1
        f              = evalF(sxx[n],syy[n],sxy[n],phi,0)
        gamma_xy[n+1]  = gamma_xy[n]+d_gxy
        # tangent modulus
                    A = (sin(beta[n])-sin(psi))*(sin(phi)-sin(beta[n]))
                    B = 2*(-cos(alpha0-beta[n])+sin(phi))*(1-nu)
        Gamma         = -A/(B*cos(beta[n])+A*cos(alpha0))*2*G

        
        # compute new sigma_v  
        sigma_v[n+1]  = sigma_v[n]+Gamma*d_gxy;
        # update stress outside the band 
        R = sigma_v[n+1]-sigh;
        sxy[n+1]      = -R*cos(alpha0)/2
        syy[n+1]      =  (sigma_v[n+1]+sigh-R*sin(alpha0))/2
        sxx0[n+1]     =  (sigma_v[n+1]+sigh+R*sin(alpha0))/2
        # update stress inside the band 
        # using consistency condition to solve for d_sxx (eq. 32 paper)
        d_syy = syy[n+1]-syy[n]
        d_sxy = sxy[n+1]-sxy[n] 
        d_sxx = -(d_syy*(sin(beta[n])+sin(phi))+2*d_sxy*cos(beta[n])+2*f)/
                    (-sin(beta[n])+sin(phi));
        sxx[n+1]=sxx[n]+d_sxx 
        # compute new orientation of stress 
        beta[n+1]   = atan(1/2*(-sxx[n+1]+syy[n+1])/sxy[n+1]);
        
        
    end
    theta = (pi/4 .+ beta./2) .*180 ./ pi
   return (sigma_v=sigma_v,gamma_xy=gamma_xy,sxx=sxx,sxx0=sxx0,sxy=sxy,syy=syy,theta=theta)
end

function evalF(sxx,syy,sxy,phi,Co)
    taustar = sqrt((sxx-syy)^2/4+sxy^2);
    sstar   = (sxx+syy)/2;
    F       = taustar+sstar*sin(phi)+Co*cos(phi);
    return F
end

Vermeer2_analytics_strain_LLP2013()