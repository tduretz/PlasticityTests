function Vermeer3_ana_llp2013(phi, psi, nu, sxx, syy, sxy, G, nstep, gamma_dot_xy,factor)
    beenhere = 0
    c = 0
    sigma_dot_yy = 0
    eps_dot_xx = 0
    gamma_dot_xy_in = gamma_dot_xy
    gamma_dot_xy_out = gamma_dot_xy

    sigma_out = zeros(3,nstep)
    sigma_out[:,1]=[sxx, syy, sxy]'
    sigma_in = zeros(3,nstep)
    sigma_in[:,1] = [sxx, syy, sxy]'
    epsilon_out = zeros(3,nstep)
    epsilon_in = zeros(3,nstep)
    alpha = zeros(nstep)
    beta  = zeros(nstep)
    theta = zeros(nstep) 
    theta_out = zeros(nstep)  
    

    L = 2 * G * nu / (1 - 2 * nu)
    D = [L + 2 * G L 0; L L + 2 * G 0; 0 0 G]

    theta[1] = 90
    theta_out[1] = theta[1]
    beta[1] = π / 2
    alpha[1] = beta[1]

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

        eps_dot_yy_in = -M[2, 3] / M[2, 2] * gamma_dot_xy_in
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

            detA = (sin(beta[i]) - sin(phi)) * (sin(beta[i]) - sin(psi)) / 4
            Z1 = -2 * (nu - 1) * (-sin(alpha[i]) * sin(phi) + cos(alpha[i] - beta[i])) * cos(beta[i])
            Z2 = (sin(beta[i]) - sin(phi)) * (sin(beta[i]) - sin(psi)) * cos(alpha[i])
            Gamma = detA / (Z2 + Z1)
            alpha_dot = 2 * G * Gamma * gamma_dot_xy_in

            sigma_dot_xx_in = 4 * alpha_dot * cos(alpha[i]) * cos(beta[i]) / (sin(beta[i]) - sin(phi))
            sigma_dot_xy_in = 2 * alpha_dot * cos(alpha[i])

            Z1 = (1 - 2 * nu) * sin(alpha[i]) * sin(psi) * (sin(beta[i]) - sin(phi))
            Z2 = (1 - nu) * cos(alpha[i]) * cos(beta[i]) * (sin(beta[i]) + sin(psi))
            Z3 = sin(alpha[i]) * sin(beta[i]) * (sin(beta[i]) - sin(phi))
            eps_dot_yy_in = (Z1 + Z2 + Z3) / detA / 2 / G * alpha_dot

            eps_dot_yy_out = -2 * alpha_dot * nu * sin(alpha[i]) / G
            sigma_dot_xx_out = 4 * alpha_dot * sin(alpha[i])
            sigma_dot_xy_out = sigma_dot_xy_in
        end

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
    return sigma_out, sigma_in,theta_out,theta, epsilon_out, epsilon_in 
end

function f(sigma, phi, c)
    tau_star = 1 / 2 * sqrt((sigma[1] - sigma[2])^2 + 4 * sigma[3]^2)
    sigma_star = (sigma[1] + sigma[2]) / 2
    return tau_star + sigma_star * sin(phi) - c * cos(phi)
end