using PlasticityTests,Plots 

function Vermeer3_ana()
    degrad = π / 180
    G = 10e6
    phi = 40 * degrad
    psi = 10 * degrad
    c = 0
    gamma_tot = 0.08
    nu = 0.1
    # IC
    sxx = -400e3
    syy = -100e3
    sxy = 0
    # CL
    dt = 1
    gamma_dot_xy = 1e-4
    nstep = Int64(gamma_tot / (gamma_dot_xy * dt))
    shearband_thickness = 1.
    model_thickness     = 10.
    factor = model_thickness/shearband_thickness
    sigma_out, sigma_in, theta_out,theta,  epsilon_out, epsilon_in = Vermeer3_ana_llp2013(phi, psi, nu, sxx, syy, sxy, G, nstep, gamma_dot_xy,factor)
# Plotting
    
p1=plot(epsilon_out[3, :], -sigma_in[3, :] ./ sigma_in[2, :], label="sxy/syy", color=:red)
p1=plot!(epsilon_out[3, :], -sigma_out[3, :] ./ sigma_out[2, :], label="sxy/syy", color=:blue)

p2=plot(epsilon_out[3, :], theta, label="theta", color=:red)
p2=plot!(epsilon_out[3, :], theta_out, label="theta", color=:blue)

p3=plot(epsilon_out[3, :], epsilon_in[2, :], label="epsilon yy", color=:red)
p3=plot!(epsilon_out[3, :], epsilon_out[2, :], label="epsilon yy", color=:blue)

p4=plot(epsilon_out[3, :], -sigma_in[1, :], label="sxx", color=:red)
p4=plot!(epsilon_out[3, :], -sigma_out[1, :], label="sxx", color=:blue)

display(plot(p1,p2,p3,p4))
end



Vermeer3_ana()