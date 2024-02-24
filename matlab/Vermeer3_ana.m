function Vermeer3_ana()
degrad = pi/180;
% parametre;
G         = 10e6;
phi       = 40*degrad;
psi       = 10*degrad;
c         =  0;
gamma_tot = 0.08;
nu        = 0. ;
col       = 'b';
ifig      = 8;
%%%%  Condition Initiale
sxx = -400e3;
syy = -100e3;
sxy = 0 ;
%%%%  CL
dt           = 1;
gamma_dot_xy = 1e-4;
eps_dot_xx   = 0;
sigma_dot_yy = 0;

nstep = gamma_tot/(gamma_dot_xy*dt);
thick  = 1/10; %one cell thick in a 10 cell model / ratio of band thickness to model thickness
[sigma_out,sigma_in,beta,alpha,epsilon_out,epsilon_in] = Vermeer3_ana_llp2013(phi,psi,nu,sxx,syy,sxy,G,nstep,gamma_dot_xy,thick);
end

function [sigma_out,sigma_in,beta,alpha,epsilon_out,epsilon_in] = Vermeer3_ana_llp2013(phi,psi,nu,sxx,syy,sxy,G,nstep,gamma_dot_xy,thick)
beenhere = 0;
%elasticity
L       = 2*G*nu/(1-2*nu);
D = [L+2*G L 0 ;L L+2*G 0; 0 0 G];
% cohesion
c= 0;
% forcing 
sigma_dot_yy =0;
eps_dot_xx   = 0;
gamma_dot_xy_in = gamma_dot_xy;
gamma_dot_xy_out = gamma_dot_xy;
%initialize arrays
sigma_out   = [sxx,syy,sxy]';
sigma_in   = [sxx,syy,sxy]';
epsilon_out   = [0 0 0]';
epsilon_in    = [0 0 0]';
gamma_bulk    = 0; 
% 
theta      = 90;
theta_out  = theta;

beta       = pi/2; 
alpha      = beta;

for i = 1:nstep
    M = D;
    F = f(sigma_in(:,i),phi,c);
    if F>=0
    dFdsig = [(sin(phi)-sin(beta(i)))/2;...
              (sin(phi)+sin(beta(i)))/2;...
              cos(beta(i))];
    dGdsig = [(sin(psi)-sin(beta(i)))/2;...
              (sin(psi)+sin(beta(i)))/2;...
              cos(beta(i))];
    a               = D*dGdsig;
    b               = D*dFdsig;
    d = dFdsig'*D*dGdsig;
    M = D-1/d*a*b';
    end
    
    eps_dot_yy_in   = - M(2,3)/M(2,2)*gamma_dot_xy_in;
    sigma_dot_xx_in = (-M(1,2)*M(2,3) + M(1,3)*M(2,2))/M(2,2)*gamma_dot_xy_in;
    sigma_dot_xy_in = (-M(3,2)*M(2,3) + M(3,3)*M(2,2))/M(2,2)*gamma_dot_xy_in;
    
    %if beta(i) >= phi*0.95
    if sigma_dot_xy_in >=0    
    %if beta(i)>atan(sin(phi)*cos(psi)/(1-(sin(phi)*sin(psi))))
    %if beta(i)>0.6283    
        eps_dot_yy_out   = - M(2,3)/M(2,2)*gamma_dot_xy_out;
        sigma_dot_xx_out = (-M(1,2)*M(2,3) + M(1,3)*M(2,2))/M(2,2)*gamma_dot_xy_out;
        sigma_dot_xy_out = (-M(3,2)*M(2,3) + M(3,3)*M(2,2))/M(2,2)*gamma_dot_xy_out;
    else
        if beenhere == 0 
            gamma_dot_xy_in=gamma_dot_xy/thick;
            beenhere = 1
        end  
          eps_dot_yy_out   = 0;
          sigma_dot_xx_out = 0; 
          sigma_dot_xy_out = sigma_dot_xy_in;
          gamma_dot_xy_out = sigma_dot_xy_in/G;
     end
    % increment stress and strain 
    sigma_in(:,i+1)   = sigma_in(:,i)     +  [sigma_dot_xx_in; sigma_dot_yy; sigma_dot_xy_in];
    epsilon_in(:,i+1) = epsilon_in(:,i)   +  [eps_dot_xx; eps_dot_yy_in; gamma_dot_xy_in];
    
    sigma_out(:,i+1)   = sigma_out(:,i)   +  [sigma_dot_xx_out; sigma_dot_yy; sigma_dot_xy_out];
    epsilon_out(:,i+1) = epsilon_out(:,i) +  [eps_dot_xx; eps_dot_yy_out; gamma_dot_xy_out];
    gamma_bulk(i+1)    = gamma_bulk(i)    +  gamma_dot_xy_in*thick+gamma_dot_xy_out*(1-thick); 

     % update alpha and beta + theta for post processing purpose 
    tau_star   = 1/2*sqrt((sigma_in(1,i+1)-sigma_in(2,i+1))^2+4*sigma_in(3,i+1)^2);
    sigma_star = (sigma_in(1,i+1)+sigma_in(2,i+1))/2;
    beta(i+1)  = asin((sigma_in(2,i+1)-sigma_in(1,i+1))/2/tau_star);
    theta(i+1)  = 45+beta(i+1)*90/pi;
    R          = 1/2*sqrt((sigma_out(1,i+1)-sigma_out(2,i+1))^2+4*sigma_out(3,i+1)^2);
    alpha(i+1)  = asin((sigma_out(2,i+1)-sigma_out(1,i+1))/2/R);
    theta_out(i+1)  = 45+alpha(i+1)*90/pi;
    
end
figure(8);clf;
subplot(221);hold on;
plot(gamma_bulk,-sigma_in(3,:)./sigma_in(2,:),'r',gamma_bulk,-sigma_out(3,:)./sigma_out(2,:),'b');title('sxy/syy');
subplot(224);hold on;
plot(gamma_bulk,theta,'r',gamma_bulk,theta_out,'b');title('theta');
subplot(223);hold on;
plot(gamma_bulk,epsilon_in(2,:),'r',gamma_bulk,epsilon_out(2,:),'b');title('epsilon yy');
subplot(222);hold on;
plot(gamma_bulk,-sigma_in(1,:),'r',gamma_bulk,-sigma_out(1,:),'b');title('sxx');
end

function F = f(sigma,phi,c)
tau_star   = 1/2*sqrt((sigma(1)-sigma(2))^2+4*sigma(3)^2);
sigma_star = (sigma(1)+sigma(2))/2;
F = tau_star + sigma_star*sin(phi)-c*cos(phi);
end