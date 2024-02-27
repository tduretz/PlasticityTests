using PlasticityTests, Plots

function Vermeer2_analytics(phi,psi,nu,alpha0,sigh,G,nstep,d_gxy)
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