
function Vermeer2_ana()
degrad = pi/180;
phi    = 40*degrad
psi    = 10*degrad
nu     = 0
alpha0 = phi * 0.999999
sigh   = -1e5
G      = 1e7
nstep  = 101
[sigma_v,gamma_xy] = Vermeer2_ana_LLP2013(phi,psi,nu,alpha0,sigh,G,nstep)
end

function [sigma_v,gamma_xy] = Vermeer2_ana_LLP2013(phi,psi,nu,alpha0,sigh,G,nstep)
S_ini       = (1+sin(phi))/(1-sin(phi));
gammatot    = 0.08;
d_gxy       = gammatot/nstep;

%initial values
sigma_v(1)   = S_ini*sigh;
beta(1)      = alpha0;
gamma_xy(1)  = 0;
R            = sigma_v(1)-sigh;
sxy(1)      = -R*cos(alpha0)/2;
syy(1)      =  (sigma_v+sigh-R*sin(alpha0))/2;
sxx0(1)     =  (sigma_v+sigh+R*sin(alpha0))/2;
sxx(1)      =  sxx0(1);

for n=1:nstep
    [f,taustar,sigmastar]      = evalF([sxx(n),syy(n),sxy(n)],phi,0);
    gamma_xy(n+1) = gamma_xy(n)+d_gxy;
    % tangent modulus
    A = (sin(beta(n))-sin(psi))*(sin(phi)-sin(beta(n)));
    B = 2*(-cos(alpha0-beta(n))+sin(phi))*(1-nu);

    Gamma(n+1)    = -A/(B*cos(beta(n))+A*cos(alpha0))*2*G;

    sigma_v(n+1)  = sigma_v(n)+Gamma(n+1)*d_gxy;
    % radius Morh circle
    R = sigma_v(n+1)-sigh;


    sxy(n+1)      = -(R)*cos(alpha0)/2;
    syy(n+1)      =  (sigma_v(n+1)+sigh-R*sin(alpha0))/2;
    sxx0(n+1)     =  (sigma_v(n+1)+sigh+R*sin(alpha0))/2;

    d_syy = diff(syy(n:n+1));
    d_sxy = diff(sxy(n:n+1));
    %%% consistency condition to solve for sxxd (eq. 32 paper)
    d_sxx = -(d_syy*(sin(beta(n))+sin(phi))+2*d_sxy*cos(beta(n))+2*f)/...
        (-sin(beta(n))+sin(phi));
    sxx(n+1)=sxx(n)+d_sxx;
    %%% compute new orientation of stress
    beta(n+1)   = atan(1/2*(-sxx(n+1)+syy(n+1))/sxy(n+1));
end
figure(1);
plot(gamma_xy,sigma_v/sigh)
end

function [F,taustar,sstar]= evalF(stress,phi,Co)
taustar = sqrt((stress(1)-stress(2))^2/4+stress(3)^2);
sstar   = (stress(1)+stress(2))/2;
F       = taustar+sstar*sin(phi);
end