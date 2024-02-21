using PlasticityTests,Plots 
function vermeer2_analytics_LLP2013()
#this program compute all the prediction of the MC model for bi-axial
#compression test with the analytical solution developped in the joint
#maple file analytics.mws
moviecircle = 1; # (put 1 if you like to see a movie of mohr circle)
#%%%%%%%%%%% physical parameters%%%%%%%%%%%%%%%%%%%%%%%%%
col      = "rbg"
phi      = 40/180*pi
psi      = 10/180*pi
alpha_v  = atan(sin(phi)*cos(psi)/(1-sin(psi)*sin(phi)))
alpha_test  = [phi*0.9999,(phi+psi)/2,psi]
nu       = 0.0
G        = 1e7
p0       = -100e3
nsteps   = 20001; 
p1=plot()
p2=plot()
p3=plot()
Test2 = ExtractDataTest2()

 for i = 1:2
     alpha_0= alpha_test[i];
#     %%%%%%%%%%%%% Initial strength%%%%%%%%%%%%%%%%%%
     S_ini   = (1+sin(phi))/(1-sin(phi))
#     %%% discretisation of the stress drop phase using effective friction %%%%%%%%%
     mu_ini  = (sin(phi)*cos(alpha_0))/(1-sin(phi)*sin(alpha_0))
     mu_res  = (sin(phi)*cos(psi))/(1-sin(phi)*sin(psi))
     phi_eff = atan.(LinRange(mu_ini,mu_res,nsteps))
     alpha   = zeros(nsteps)
     S0       = zeros(nsteps)
     S       = zeros(nsteps)
     Sxx_out = zeros(nsteps-1)
     Sxx_in  = zeros(nsteps-1)
     Gxy_in  = zeros(nsteps)
     Z1      = zeros(nsteps-1)
     Z2      = zeros(nsteps-1)
     Z3      = zeros(nsteps-1)
     dgxy_in = zeros(nsteps-1)
     dsxx_in  = zeros(nsteps-1)
     dsxx_out = zeros(nsteps-1)
     GAMMA    = zeros(nsteps-1)
#     %%%%%%%%%%%%%%% Compute S and alpha from phi_eff%%%%%
 @.  alpha   = acos.(sin.(phi_eff).*cos.(phi_eff)./sin(phi).*
              (1 .+sqrt.((cos.(phi_eff).^2 .- cos(phi)^2)./cos.(phi_eff).^2)));
 @.  S       = (cos(alpha_0)+tan.(phi_eff)*(sin(alpha_0)+1))./
               (cos(alpha_0)+tan.(phi_eff)*(sin(alpha_0)-1));
               @show S[1],S[end]
#     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     S_res   = (mu_res*(sin(alpha_0)+1)+cos(alpha_0))/
               (mu_res*(sin(alpha_0)-1)+cos(alpha_0))
     checkSres = S_res - S[end]
#     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     dsdt  = diff(S)
     @show dsdt
#     %%%%%%%%%%%%%%%%% center the discretisation %%%%%%%%%%%%%%%%
     alpha = alpha[1:end-1] .+ diff(alpha,dims=1) ./ 2
     S     = S[1] .+ cumsum(dsdt,dims=1)
#     %%%%%%%%%%%%%%%%%% Compute Stress out from equation ???? %%%%%%
@.    Sxx_out  = 1/2*p0*(1+S+sin(alpha_0)*S-sin(alpha_0));
#     Syy_out  =-1/2*p0*(-S-1+sin(alpha_0)*S-sin(alpha_0));
#     Sxy_out  =-1/2*p0*(S-1)*cos(alpha_0);
#     %%%%%%%%%%%%%%%%stress rate and strain rate outside%%%%%%%%%%%%
#     dsxx_out = 1/2*p0*dsdt*(1+sin(alpha_0));
#     dsyy_out = 1/2*p0*dsdt*(1-sin(alpha_0));
#     dsxy_out =-1/2*p0*dsdt*cos(alpha_0);
#     dexx_out = 1/4*p0*dsdt*(1+sin(alpha_0)-2*nu)/G;
#     deyy_out =-1/4*p0*dsdt*(2*nu-1+sin(alpha_0))/G;
#     dgxy_out =-1/2*1/G*p0*dsdt*cos(alpha_0);
#     %%%%%%%%%%%%%%%%%%%%%%Fancy formula %%%%%%%%%%%%%%%%%%%%%%%%%%%
    @.  Z1 = 1/2*(cos(alpha_0-alpha)-sin(phi));
    @.  Z2 = (sin(phi)-sin(alpha)).*(sin(alpha)-sin(psi));
    @.  Z3 = 4*Z1.*cos(alpha)*(1-nu)-Z2*cos(alpha_0);
#     %%%%%%analytics from strain formulation on maple%%%%%%%%%%%
    @. GAMMA    = 2*Z2./Z3*G/p0;
    @show GAMMA
#     Lambda   = 2*Z1*(1-nu).*dsdt*p0./Z2/G;


#     %%%%%%%%%%% INCREMENTAL ANALYTICS %%%%%%%%%%%%%%%%%%
#     %%%%%%%%%  strain increments%%%%%%%%
#     dexx_in  = dexx_out; % compatibility
#     deyy_in  = deyy_out+p0/G*Z1.*dsdt.*((1-2*nu)*sin(alpha)+sin(psi))./Z2; %%% here%%
      dgxy_in  = dsdt./GAMMA;
#     %%%%%%%%%% stress increments%%%%%%%%%%%%%
@.    dsxx_in = dsxx_out+2*p0*dsdt.*Z1.*(sin(alpha)-sin(psi))./Z2
#     dsyy_in = dsyy_out; %continuity
#     dsxy_in = dsxy_out; %continuity

#     %%%%%%%%%% MATLAB CUMSUM %%%%%%%%%%%%%%%%%%%%%%%%%
#     %%%%%%%%%%  cumm. strain%%%%%%%%%%%%%%%%%%%%%%
#     Exx_in  = cumsum(dexx_in);
#     Eyy_in  = cumsum(deyy_in);
      Gxy_in  = cumsum(dgxy_in);
#     Exx_out  = cumsum(dexx_out);
#     Eyy_out  = cumsum(deyy_out);
#     Gxy_out = cumsum(dgxy_out);
#     %%%%%%%%%% stress %%%%%%%%%%%%%%%%%%%%%%
      Sxx_in  .= cumsum(dsxx_in) .+ Sxx_out[1]
#     Syy_in  = Syy_out;
#     Sxy_in  = Sxy_out;
#     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     %%%%%%%%%%%%%%%%%% Post Processing %%%%%%%%%%
#     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#     figure(1);hold on;
#     plot(Gxy_in,S,'k');
#     plot(Gxy_in,Sxx_in,'r');
#     hold off;


      p1=plot!(Gxy_in,S,xlims=(-0.01,0.08),ylims=(3.7,5)) 
     # p1=scatter!(Test2.StressRatioArthur.x./100,Test2.StressRatioArthur.y,xlims=(-0.01,0.08),ylims=(3.7,5))
     # p1=scatter!(Test2.StressRatioCoulomb.x./100,Test2.StressRatioCoulomb.y,xlims=(-0.01,0.08),ylims=(3.7,5))
#     figure(8);
#     subplot(211);hold on;
     # p2 = plot(Gxy_in,(pi/4 .+ alpha/2)/pi*180,xlims=(0,0.08))
#     subplot(212);hold on;
     # p3 = plot(Gxy_in,Sxx_in./Sxx_out,xlims=(0,0.08))

      

#     if moviecircle == 1
#         C_in  = (Sxx_in+Syy_in)/2;
#         C_out = (Sxx_out+Syy_out)/2;
#         R_in  = sqrt(((Sxx_in-Syy_in)/2).^2+Sxy_out.^2);
#         R_out = sqrt(((Sxx_out-Syy_out)/2).^2+Sxy_out.^2);
#         theta = linspace(0,pi,100);
#         for i = 1:199
#             xin = C_in(i)+R_in(i)*cos(theta);
#             yin = R_in(i)*sin(theta);
#             xout = C_out(i)+R_out(i)*cos(theta);
#             yout = R_out(i)*sin(theta);
#             figure(2); clf;hold on;
#             plot([0,-5],-tan(phi)*[0,-5],'k');
#             plot(xin,yin,'r');
#             plot(xout,yout,'b');
#             plot(Sxx_out(i),Sxy_out(i),'k+');plot(Syy_out(i),Sxy_out(i),'g+');
#             plot(Sxx_in(i),Sxy_in(i),'m+'); axis equal
#             hold off; pause (0.01);
#         end
#     end
 end
 p1=scatter!(Test2.StressRatioArthur.x./100,Test2.StressRatioArthur.y,xlims=(-0.01,0.16),ylims=(3.7,5))
 p1=scatter!(Test2.StressRatioCoulomb.x./100,Test2.StressRatioCoulomb.y,xlims=(-0.01,0.16),ylims=(3.7,5))
 display(plot(p1))

end
vermeer2_analytics_LLP2013()