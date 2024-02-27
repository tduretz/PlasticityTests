using PlasticityTests, GeoParams, Plots, Printf, MathTeXEngine, LinearAlgebra, StaticArrays, Statistics
import LinearAlgebra:norm

function main_Yury()

    # Physical parameters
    # dimenzionally independent
    G        = 1e7
    Ly       = 1e+0
    eps0     = 1e-4
    # nondimentional
    phi      = 40.0/180.0*π
    psi      = 10.0/180.0*π
    sigh_G   = -1e-2
    # dimenzionally dependent
    Coh      = 0*G
    sigh     = sigh_G*G
    sigv     = (1.0 + sin(phi))/(1.0 - sin(phi))*sigh/1.0
    # Numerical parameters
    nt       = 150
    ny       = 10
    dt       = 1e-4/eps0/1
    etavp    = 0.0.*G/eps0
    # PT solver
    niter    = 1e5
    thetaVx  = 6/ny
    thetaVy  = 6/ny
    thetaPt  = 1.0
    nout     = 1000
    eps2     = 1e-8*G/Ly
    errPt    = 0.
    errVx    = 0.
    errVy    = 0.
    relVy    = 30*Ly*eps0/G
    # preprocessing
    dy       = Ly/ny
    yc       = LinRange(-Ly/2-dy/2, Ly/2+dy/2, ny+2)
    yv       = LinRange(-Ly/2,      Ly/2,      ny+1)
    theta_A  = pi/4 + 0.25*(phi + psi)
    theta_C  = pi/4 + 0.5*(phi)
    theta_R  = pi/4 + 0.5*(psi)
    theta_SB = theta_C
    Dev      = I(3)-1/3*ones(3,3)

    alpha_0     = 39.999999*(π/180)
    nstep       = 1000
    d_gxy       = 2*0.08/nstep
    nu          = 0.
    Coulomb     = Vermeer2_analytics(phi,psi,nu,alpha_0,-100e3,G,nstep,d_gxy)

    # to Cartesian
    sigxxi   = 1/2*(sigh + sigv) +  1/2*(sigh - sigv)*cos(2*theta_SB)
    sigyyi   = 1/2*(sigh + sigv) -  1/2*(sigh - sigv)*cos(2*theta_SB)
    sigxyi   = 1/2*(sigh - sigv)*sin(2*theta_SB)
    sigzzi   = 0.5*(sigxxi + sigyyi)
    Pri      = -(sigxxi + sigyyi)/2.0
    tauxxi   = Pri + sigxxi
    tauyyi   = Pri + sigyyi
    tauzzi   = Pri + sigzzi
    tauxyi   = copy(sigxyi)
    etae     = G*dt
    etave    = etae;#1.0./(1.0/mus  + 1.0/etae);
    # initial conditions
    Pt        =  Pri*ones(ny+1); Pt[Int64(ny/2)] = Pri/1.001
    Ptc       =  copy(Pt)
    Pt0       =  copy(Pt)
    tauxx     =  tauxxi*ones(ny+1); tauxx0 = copy(tauxx)
    tauxy     =  tauxyi*ones(ny+1); tauxy0 = copy(tauxy)
    tauyy     =  tauyyi*ones(ny+1); tauyy0 = copy(tauyy)
    tauzz     =  tauzzi*ones(ny+1); tauzz0 = copy(tauzz)
    tau_new   =  zeros(3, ny+1)
    tauii     =  zeros(ny+1)
    F         =  zeros(ny+1)
    Fc        =  zeros(ny+1)
    Vx        =  collect(eps0.*yc)
    Vy        =  zeros(ny+2)
    divV      =  zeros(ny+1)
    epsxx     =  zeros(ny+1)
    epsyy     =  zeros(ny+1)
    epszz     =  zeros(ny+1)
    epsxy     =  zeros(ny+1)
    sigxx_dot =  zeros(ny+1)
    theta     =  zeros(ny+1)
    lamrel    =  zeros(ny+1)
    epsxx_pl  =  zeros(ny+1)
    epsyy_pl  =  zeros(ny+1)
    epszz_pl  =  zeros(ny+1)
    epsxy_pl  =  zeros(ny+1)
    divV_pl   =  zeros(ny+1)
    epsxyt    =  zeros(ny+1)
    epsyyt    =  zeros(ny+1)
    RVx       =  zeros(ny+2)
    RVy       =  zeros(ny+2)
    RPt       =  zeros(ny+1)
    dPtdtau   =  zeros(ny+1)
    dVxdtau   =  zeros(ny+2)
    dVydtau   =  zeros(ny+2)
    # BC
    VxS       =  Vx[1]
    VxN       =  Vx[end]
    VyS       =  0.0
    VyN       =  0.0
    sigv_evol = zeros(nt)
    sigh_evol = zeros(nt)
    gamxy     = zeros(nt)
    theta_i   = zeros(nt)
    theta_o   = zeros(nt)
    # action
    for it=1:nt
        # History
        @. lamrel   = 0.
        @. tauxy0  = tauxy
        @. tauxx0  = tauxx
        @. tauyy0  = tauyy
        @. tauzz0  = tauzz
        @. Pt0     = Pt
        sigv0      = sigv
        sigv_evol[it] = tauxx[end]*sin(theta_SB).^2 + tauyy[end]*cos(theta_SB).^2 - tauxy[end]*sin(2*theta_SB) - Pt[end]
        sigh_evol[it] = tauxx[end]*cos(theta_SB).^2 + tauyy[end]*sin(theta_SB).^2 + tauxy[end]*sin(2*theta_SB) - Pt[end]
        for iter=1:niter
            errBC       = tauxx[end]*sin(theta_SB).^2 + tauyy[end]*cos(theta_SB).^2 - tauxy[end]*sin(2*theta_SB) - Pt[end] - sigv
            sigh_it     = tauxx[end]*cos(theta_SB).^2 + tauyy[end]*sin(theta_SB).^2 + tauxy[end]*sin(2*theta_SB) - Pt[end]
            sigv        += 5e-1*(sigh - sigh_it)
            Vy[end]     = Vy[end] - relVy*errBC
            Vx[1]       = -Vx[2]      + 2*VxS
            Vx[end]     = -Vx[end-1]  + 2*VxN
            Vy[1]       = - Vy[2]     + 2*VyS
            @. theta     = 0.5*acos((tauxx-tauyy)/2/tauxy)
            @. sigxx_dot = 1/2*(1-cos(2*theta)) * (sigv - sigv_evol[it])/dt 
            @. epsxx    = -(sigxx_dot) /2/G 
            @. epsxy    =  0.5*(Vx[2:end] - Vx[1:end-1])/dy
            @. epsyy    =      (Vy[2:end] - Vy[1:end-1])/dy
            @. epszz    = 1/2*(epsxx + epsyy);
            @. divV     = epsxx + epsyy + epszz;
            # Stress
            tau_new .= etave./etae.*[tauxx0'; tauyy0'; tauzz0'] + 2.0.*etave.*Dev*[epsxx'; epsyy'; epszz' ]
            @. tauxx       = tau_new[1,:]
            @. tauyy       = tau_new[2,:]
            @. tauzz       = tau_new[3,:]
            @. tauxy       = tauxy0 + 2 * etave * epsxy;
            @. tauii       = sqrt(tauxy.^2 + 0.5*(tauyy.^2 + tauxx.^2 + tauzz.^2));

            # Plasticity
            @. F           = tauii - Coh*cos(phi) - Pt*sin(phi)
            @. Ptc         = Pt

            for itplast = 1:50
                @. epsxx_pl = lamrel.*(tauxx/2/tauii)
                @. epsyy_pl = lamrel.*(tauyy/2/tauii)
                @. epszz_pl = lamrel.*(tauxx/2/tauii + tauyy/2/tauii)/2
                @. epsxy_pl = lamrel.*(tauxy/2/tauii)
                @. divV_pl = 2/3*sin(psi)*lamrel
                @. Ptc      = Pt0 - 3/2*G*dt*(divV - divV_pl)
                tau_new .= etave./etae.*[tauxx0'; tauyy0'; tauzz0'] +    2.0.*etave.*Dev*[epsxx' - epsxx_pl'; epsyy' - epsyy_pl'; epszz' - epszz_pl']
                @. tauxx       = tau_new[1,:]
                @. tauyy       = tau_new[2,:]
                @. tauzz       = tau_new[3,:]
                @. tauxy       = tauxy0 +  2 * etave * (epsxy - epsxy_pl)
                @. tauii  = sqrt(tauxy.^2 + 0.5*(tauyy.^2 + tauxx.^2 + tauzz.^2))
                @. Fc     = tauii - Coh*cos(phi) - Ptc*sin(phi) - etavp*lamrel
                @. lamrel += (F>0) .* Fc / (etavp + etave + 2/3*G*dt*sin(phi)*sin(psi))
                if (maximum(Fc) < eps2) break; end
            end
            # PT time steps
            deltatauV    = dy.^2/etae/2.1 /4
            deltatauPt   = 3/2      # (2/3*G*dt*deltatauV/dy^2)*deltatauPt = 1,  = 2.1*4*etae/(2/3*G*dt)
            # @show deltatauV, deltatauPt
            # Residuals
            @. RPt          =  (- 2/3*G*dt*divV - (Pt - Pt0))
            @. RVx[2:end-1] =  ((tauxy[2:end] - tauxy[1:end-1])/dy )
            @. RVy[2:end-1] =  ((tauyy[2:end] - tauyy[1:end-1])/dy - (Ptc[2:end] - Ptc[1:end-1])/dy)
            # Damp residuals
            @. dVxdtau      = RVx + (1.0 - thetaVx)*dVxdtau
            @. dVydtau      = RVy + (1.0 - thetaVy)*dVydtau
            @. dPtdtau      = RPt + (1.0 - thetaPt)*dPtdtau 
            # Update solutions
            @. Vx[2:end-1]  += deltatauV  * dVxdtau[2:end-1]
            @. Vy[2:end-1]  += deltatauV  * dVydtau[2:end-1]
            @. Pt           += deltatauPt * dPtdtau 

             sigv = sigv_evol[it]
             sigh = sigh_evol[it]

            if mod(iter, nout) == 0 || iter==1
                @show maximum(abs.(RVx))
                @show maximum(abs.(RVy))
                errPt = maximum(abs.(RPt))
                errVx = maximum(abs.(RVx))
                errVy = maximum(abs.(RVy))
                σyyBC = tauyy[end] - Ptc[end]
                @printf("Iteration %05d --- Time step %4d --- εxy = %2.2e --- σyyBC = %2.7e --- max(F) = %2.2e --- max(Fc) = %2.2e \n", iter, it, eps0*it*dt, σyyBC*G, maximum(F)*G, maximum(Fc)*G )
                @printf("fPt = %2.4e\n", errPt)
                @printf("fVx = %2.4e\n", errVx)
                @printf("fVy = %2.4e\n", errVy)
                σh = (tauxx[end]-Pt[end])*cos(theta_SB)^2 + (tauyy[end]-Pt[end])*sin(theta_SB)^2 + (tauxy[end])*sin(2*theta_SB)
                @printf("σh = %2.4e\n", σh/1e3)
                (errVx < eps2 && errVy < eps2) && break 
                (isnan(errPt) || isnan(errVx) || isnan(errVx)) && error("NaNs") 
            end
        end
        @. Pt      = Ptc
        @. epsxyt += epsxy*dt
        @. epsyyt += epsyy*dt
        gamxy[it]  = maximum(epsxyt*2)

        @show epsxyt

        sigv_evol[it] = sigv
        sigh_evol[it] = sigh

        # method 2 - angle outside
        sigxx_o  = 1/2*(sigh + sigv) +  1/2*(sigh - sigv)*cos(2*theta_SB)
        sigyy_o  = 1/2*(sigh + sigv) -  1/2*(sigh - sigv)*cos(2*theta_SB)
        sigxy_o  = 1/2*(sigh - sigv)*sin(2*theta_SB)
        sig_o_v2 = [sigxx_o; sigyy_o; sigxy_o]
        tau_o_v2 = sqrt(0.25*(sig_o_v2[1]-sig_o_v2[2])^2 + sig_o_v2[3]^2)
        theta_o[it] = 0.5*acos((sig_o_v2[1]-sig_o_v2[2])/2/tau_o_v2)
        
        # method 2 - angle inside
        _, imin = findmin(Pt)
        sigxx_i = -Pt[imin] + tauxx[imin]
        sig_i = [sigxx_i; sigyy_o; sigxy_o]
        tau_i = sqrt(0.25*(sig_i[1]-sig_i[2])^2 + sig_i[3]^2)
        theta_i[it]  = 0.5*acos((sig_i[1]-sig_i[2])/2/tau_i)

        if mod(it,10)==0
            p1 = plot(Vx, yc)
            p2 = plot(xlabel="strain", ylabel="sv/sh")
            p2 = plot!(gamxy,sigv_evol./sigh_evol)
            p2 = plot!(Coulomb.gamma_xy, Coulomb.sigma_v/(-100e3))
            p3 = plot(xlabel="strain", ylabel="angle")
            p3 = plot!(gamxy[1:it], theta_o[1:it]*180/pi, label="out")
            p3 = plot!(gamxy[1:it], theta_i[1:it]*180/pi, label="in")
            display(plot(p1, p2, p3))
        end
    end
end

main_Yury()