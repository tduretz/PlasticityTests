function Vermeer3_ana_llp2013_v2(σi, params, factor)
    
    sxx = σi.xx
    syy = σi.yy
    sxy = σi.xy

    sxx_in = σi.xx
    syy_in = σi.yy
    sxy_in = σi.xy
    
    beenhere         = 0
    dt               = params.Δt
    eta_vp           = params.ηvp
    nstep            = params.nt
    G                = params.G
    K                = params.K 
    c                = params.c
    phi              = params.ϕ
    psi              = params.ψ
    dsyy       = 0.
    dexx       = 0.
    dγxy_in  = params.γ̇xy*params.Δt
    dγxy_out = params.γ̇xy*params.Δt

    sigma_out = zeros(3,nstep)
    sigma_out[:,1]=[sxx, syy, sxy]'
    sigma_in = zeros(3,nstep)
    sigma_in[:,1] = [sxx, syy, sxy]'
    epsilon_out = zeros(3,nstep)
    epsilon_in  = zeros(3,nstep)
    alpha       = zeros(nstep)
    beta        = zeros(nstep)
    theta       = zeros(nstep) 
    theta_out   = zeros(nstep)
    gamma_bulk  = zeros(nstep)  
    dsxy_v      = zeros(nstep)

    nu = (3K-2G)/(6K+2G)
    L  = 2 * G * nu / (1 - 2 * nu)
    D  = [L + 2 * G L 0; L L + 2 * G 0; 0 0 G]
    @show D

    if sxx>syy
        theta[1]     = 0.0
    else
        theta[1]     = 90.0
    end
    theta_out[1] = theta[1]
    beta[1]      = π / 2    
    alpha[1]     = beta[1]

    sxy_in_v = 0.
    dγxyp    = 0.
    dσe      = 0.

    for i = 1:nstep-1
        M = copy(D)
        dγxyp0    = dγxyp
        dσe0      = dσe
        sxy_in_v0 = sxy_in_v
        sxx_in0   = sxx_in
        syy_in0   = syy_in
        sxy_in0   = sxy_in

        # Trial stress
        deyy_in = 0.0 
        dsxx_in = 0.0
        dsxy_in =  D[3, 3] * (dγxy_in)
        dsyy_in = dsyy 

        sxx_in = sxx_in0 + dsxx_in
        syy_in = syy_in0 + dsyy_in
        sxy_in = sxy_in0 + dsxy_in

        σ = [sxx_in; syy_in; sxy_in] 

        dγ = 0
        F  = f(σ, phi, c, dγ, dt, eta_vp)
        if F  >= 0
    
            for it=1:50
                tau_star = sqrt(1/4*(σ[1] - σ[2])^2 + σ[3]^2)
                dexxp = dγ*(1/4*(σ[1]-σ[2])/tau_star + 1/2*sin(psi))
                deyyp = dγ*(1/4*(σ[2]-σ[1])/tau_star + 1/2*sin(psi))
                dγxyp = dγ*(σ[3]/tau_star)

                deyy_in =  deyyp
                dsxx_in = -D[1,1] * dexxp
                dsyy_in = dsyy 
                dsxy_in =  D[3,3] * (dγxy_in - dγxyp)

                sxx_in = sxx_in0 + dsxx_in
                syy_in = syy_in0 + dsyy_in
                sxy_in = sxy_in0 + dsxy_in
        
                σ = [sxx_in; syy_in; sxy_in] 

                F  = f(σ, phi, c, dγ, dt, eta_vp)

                dγ += F/(G + K*sin(phi)*sin(psi) + eta_vp/dt)

                if abs(F)< 1e-5 break; end 
            end
        end

        sxy_in_v  = eta_vp*(dγxyp/2)/ dt 
        dsxy_in_v = (sxy_in_v - sxy_in_v0)/dt

        dσv = (dγxyp - dγxyp0)/dt
        dσe = D[3, 3] * (dγxy_in-dγxyp) / dt

        dsxy_in_v = (dσe - dσe0)/ dt

        if dsxy_in_v <0 && dσe>0 # dsxy_in_v >= 0
            deyy_out = deyy_in
            dsxx_out = dsxx_in
            dsyy_in  = dsyy 
            dsxy_out = dsxy_in
        else
            if beenhere == 0
                dγxy_in *= factor
                beenhere = 1
            end
            deyy_out   = 0 # elastic deformation with constant sigma_yy and sigma_xy
            dsxx_out = 0 # epsilon xx is 0 in 1D simple shear and dsyy = 0 
            dsxy_out = dsxy_in # continuity 
            dγxy_out = dsxy_in/G # elastic deformation
        end
        gamma_bulk[i+1]       = gamma_bulk[i]  +  dγxy_in/factor+dγxy_out*(1-1/factor)     

        sigma_in[:, i + 1]    = sigma_in[:, i]   + [dsxx_in, dsyy, dsxy_in]
        epsilon_in[:, i + 1]  = epsilon_in[:, i] + [dexx, deyy_in, dγxy_in]

        sigma_out[:, i + 1]   = sigma_out[:, i]  + [dsxx_out, dsyy, dsxy_out]
        epsilon_out[:, i + 1] = epsilon_out[:, i]+ [dexx, deyy_out, dγxy_out]
 
        tau_star         = 1 / 2 * sqrt((sigma_in[1, i + 1] - sigma_in[2, i + 1])^2 + 4 * sigma_in[3, i + 1]^2)
        beta[i + 1]      = asin((sigma_in[2, i + 1] - sigma_in[1, i + 1]) / 2 / tau_star)
        theta[i + 1]     = 45 + beta[i + 1] * 90 / π
        R                = 1 / 2 * sqrt((sigma_out[1, i + 1] - sigma_out[2, i + 1])^2 + 4 * sigma_out[3, i + 1]^2)
        alpha[i + 1]     = asin((sigma_out[2, i + 1] - sigma_out[1, i + 1]) / 2 / R)
        theta_out[i + 1] = 45 + alpha[i + 1] * 90 / π
        dsxy_v[i+1] =  dsxy_in_v
    end
    return (σ_out=sigma_out, σ_in=sigma_in, θ_out=theta_out, θ_in=theta, ε_out=epsilon_out, ε_in=epsilon_in,γ_bulk= gamma_bulk, dsxy_v=dsxy_v ) 
end

function f(sigma, phi, c, dγ, dt, eta_vp)
    # tau_star = 1 / 2 * sqrt((sigma[1] - sigma[2])^2 + 4 * sigma[3]^2)
    tau_star = sqrt(1/4*(sigma[1] - sigma[2])^2 + sigma[3]^2)
    sigma_star = (sigma[1] + sigma[2]) / 2
    return tau_star + sigma_star * sin(phi) - c * cos(phi) - dγ/dt*eta_vp
end