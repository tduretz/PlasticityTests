function Vermeer3_ana_llp2013(σi, params, factor)
    
    sxx = σi.xx
    syy = σi.yy
    sxy = σi.xy
    
    beenhere         = 0
    nstep            = params.nt
    G                = params.G
    K                = params.K 
    c                = params.c
    phi              = params.ϕ
    psi              = params.ψ
    sigma_dot_yy     = 0.
    eps_dot_xx       = 0.
    gamma_dot_xy_in  = params.γ̇xy*params.Δt
    gamma_dot_xy_out = params.γ̇xy*params.Δt

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

    nu = (3K-2G)/(6K+2G)
    L  = 2 * G * nu / (1 - 2 * nu)
    D  = [L + 2 * G L 0; L L + 2 * G 0; 0 0 G]

    if sxx>syy
        theta[1]     = 0.0
    else
        theta[1]     = 90.0
    end
    theta_out[1] = theta[1]
    beta[1]      = π / 2
    alpha[1]     = beta[1]

    for i = 1:nstep-1
        M = copy(D)

        F = f(sigma_in[:, i], phi, c)

        if F >= 0
            dFdsig = [(sin(phi) - sin(beta[i])) / 2;
                      (sin(phi) + sin(beta[i])) / 2;
                      cos(beta[i])]

            dGdsig = [(sin(psi) - sin(beta[i])) / 2;
                      (sin(psi) + sin(beta[i])) / 2;
                      cos(beta[i])]

            a = D * dGdsig
            b = D * dFdsig

            d = dFdsig' * D * dGdsig
            M = D - 1 / d * a * b'
        end

        eps_dot_yy_in   = -M[2, 3] / M[2, 2] * gamma_dot_xy_in
        sigma_dot_xx_in = (-M[1, 2] * M[2, 3] + M[1, 3] * M[2, 2]) / M[2, 2] * gamma_dot_xy_in
        sigma_dot_xy_in = (-M[3, 2] * M[2, 3] + M[3, 3] * M[2, 2]) / M[2, 2] * gamma_dot_xy_in

        if sigma_dot_xy_in >= 0
            eps_dot_yy_out = -M[2, 3] / M[2, 2] * gamma_dot_xy_out
            sigma_dot_xx_out = (-M[1, 2] * M[2, 3] + M[1, 3] * M[2, 2]) / M[2, 2] * gamma_dot_xy_out
            sigma_dot_xy_out = (-M[3, 2] * M[2, 3] + M[3, 3] * M[2, 2]) / M[2, 2] * gamma_dot_xy_out
        else
            if beenhere == 0
                gamma_dot_xy_in *= factor
                beenhere = 1
            end

            eps_dot_yy_out   = 0 # elastic deformation with constant sigma_yy and sigma_xy
            sigma_dot_xx_out = 0 # epsilon xx is 0 in 1D simple shear and sigma_dot_yy = 0 
            sigma_dot_xy_out = sigma_dot_xy_in # continuity 
            gamma_dot_xy_out = sigma_dot_xy_in/G # elastic deformation
        end
        gamma_bulk[i+1]      = gamma_bulk[i]  +  gamma_dot_xy_in/factor+gamma_dot_xy_out*(1-1/factor)     

        sigma_in[:, i + 1]    = sigma_in[:, i]   + [sigma_dot_xx_in, sigma_dot_yy, sigma_dot_xy_in]
        epsilon_in[:, i + 1]  = epsilon_in[:, i] + [eps_dot_xx, eps_dot_yy_in, gamma_dot_xy_in]

        sigma_out[:, i + 1]   = sigma_out[:, i]  + [sigma_dot_xx_out, sigma_dot_yy, sigma_dot_xy_out]
        epsilon_out[:, i + 1] = epsilon_out[:, i]+ [eps_dot_xx, eps_dot_yy_out, gamma_dot_xy_out]

        tau_star = 1 / 2 * sqrt((sigma_in[1, i + 1] - sigma_in[2, i + 1])^2 + 4 * sigma_in[3, i + 1]^2)
        beta[i + 1] = asin((sigma_in[2, i + 1] - sigma_in[1, i + 1]) / 2 / tau_star)
        theta[i + 1] = 45 + beta[i + 1] * 90 / π
        R = 1 / 2 * sqrt((sigma_out[1, i + 1] - sigma_out[2, i + 1])^2 + 4 * sigma_out[3, i + 1]^2)
        alpha[i + 1] = asin((sigma_out[2, i + 1] - sigma_out[1, i + 1]) / 2 / R)
        theta_out[i + 1] = 45 + alpha[i + 1] * 90 / π
    end
    return (σ_out=sigma_out, σ_in=sigma_in, θ_out=theta_out, θ_in=theta, ε_out=epsilon_out, ε_in=epsilon_in,γ_bulk= gamma_bulk ) 
end

function f(sigma, phi, c)
    tau_star = 1 / 2 * sqrt((sigma[1] - sigma[2])^2 + 4 * sigma[3]^2)
    sigma_star = (sigma[1] + sigma[2]) / 2
    return tau_star + sigma_star * sin(phi) - c * cos(phi)
end