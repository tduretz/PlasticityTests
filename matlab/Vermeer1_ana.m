function Vermeer1_ana()
degrad = pi/180;
% parametre;
G         = 10e6;
phi       = 40;
psi       = 10;
c         =  0;
gamma_tot = 0.08;
nu        = 0.1 ;
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
[sigma,beta,epsilon] = Vermeer1_ana_llp2013(phi,psi,nu,sxx,syy,sxy,G,nstep,gamma_dot_xy)
end

function [sigma,beta,epsilon] = Vermeer1_ana_llp2013(phi,psi,nu,sxx,syy,sxy,G,nstep,gamma_dot_xy)
c= 0;
sigma_dot_yy =0;
eps_dot_xx   = 0;
sigma   = [sxx,syy,sxy]';
%   [exx, eyy, gammaxy]'
epsilon = [0 0 0]';
L       = 2*G*nu/(1-2*nu);
D = [L+2*G L 0 ;L L+2*G 0; 0 0 G];
% rayon du cercle de morh
tau_star   = 1/2*sqrt((sigma(1)-sigma(2))^2+4*sigma(3)^2);
% centre du cercle de morh
sigma_star = (sigma(1)+sigma(2))/2;

cos2theta  = (sigma(1)-sigma(2))/2/tau_star;
sin2theta  = sigma(3)/tau_star;
theta      = 1/2*asind(sin2theta);
beta       = pi/2 ; 

for i = 1:nstep
    
    M = D;
   
    F = f(sigma(:,i),phi,c)
    if F>=0
    dFdsig = [(sind(phi)-sin(beta(i)))/2;...
              (sind(phi)+sin(beta(i)))/2;...
              cos(beta(i))];
    dGdsig = [(sind(psi)-sin(beta(i)))/2;...
              (sind(psi)+sin(beta(i)))/2;...
              cos(beta(i))];
    
    % dGdsig          = [ 1/2*cos2theta + 1/2*sind(psi);...
    %                    -cos2theta/2 + sind(psi)/2;...
    %                     sin2theta];
    a               = D*dGdsig;
    
    
    
    % dFdsig          = [ 1/2*cos2theta + 1/2*sind(phi);...
    %                    -cos2theta/2 + 1/2*sind(phi);...
    %                     sin2theta];
    b               = D*dFdsig;
    
    d = dFdsig'*D*dGdsig;
    
    M = D-1/d*a*b';
    end
    
    eps_dot_yy   = - M(2,3)/M(2,2)*gamma_dot_xy;
    sigma_dot_xx = (-M(1,2)*M(2,3) + M(1,3)*M(2,2))/M(2,2)*gamma_dot_xy;
    sigma_dot_xy = (-M(3,2)*M(2,3) + M(3,3)*M(2,2))/M(2,2)*gamma_dot_xy;
    
  
    
    % a chaque pas de temps j'ajoute l'increment de contrainte et de
    % deformation
    sigma(:,i+1)   = sigma(:,i)   +  [sigma_dot_xx; sigma_dot_yy; sigma_dot_xy];
    epsilon(:,i+1) = epsilon(:,i) +  [eps_dot_xx; eps_dot_yy; gamma_dot_xy];
    
    tau_star   = 1/2*sqrt((sigma(1,i+1)-sigma(2,i+1))^2+4*sigma(3,i+1)^2);
    sigma_star = (sigma(1,i+1)+sigma(2,i+1))/2;
    
    cos2theta  = (sigma(1,i+1)-sigma(2,i+1))/2/tau_star;
    sin2theta  = sigma(3,i+1)/tau_star;
    theta(i+1) = 1/2*acosd(cos2theta);
    beta(i+1)  = asin((sigma(2,i+1)-sigma(1,i+1))/2/tau_star);
    
end
figure(8);
subplot(221);hold on;
plot(epsilon(3,:),-sigma(3,:)./sigma(2,:));title('sxy/syy');
subplot(224);hold on;
plot(epsilon(3,:),theta);title('theta');
subplot(223);hold on;
plot(epsilon(3,:),epsilon(2,:));title('gamma yy');
subplot(222);hold on;
plot(epsilon(3,:),-sigma(1,:));title('sxx');
end

function F = f(sigma,phi,c)
tau_star   = 1/2*sqrt((sigma(1)-sigma(2))^2+4*sigma(3)^2);
sigma_star = (sigma(1)+sigma(2))/2;
F = tau_star + sigma_star*sind(phi)-c*cosd(phi);
end