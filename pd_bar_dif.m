function [dif] = pc_bar_dif(x,Params)

[gamma,psi,beta,pc,pd,kap_1,kap_0,kap_1m,kap_0m,mu_c,mu_d,rho_x,...
    rho_d,rho_cd,phi_x,phi_c,phi_d,phi_rm,phi_rf,phi_s,sig_h,pi] = unpack(Params);
theta = (1-gamma)/(1-1/psi);

gx_pc   = (1-1/psi)/(1-kap_1*rho_x);
kap_1m  = exp(x)/(1+exp(x));
kap_0m  = log(1+exp(x))-kap_1m*x;
gx_pd   = (rho_d-1/psi)/(1-kap_1m*rho_x);
g0_pd   = (theta*log(beta)+(theta-1-theta/psi)*mu_c+(theta-1)*kap_0+kap_0m+mu_d+...
    (theta-1)*(kap_1-1)*pc+0.5*(((theta-1-theta/psi+rho_cd)*phi_c)^2+(((theta-1)*kap_1*gx_pc+kap_1m*gx_pd)*phi_x)^2+phi_d^2)*(pi*sig_h^2+(1-pi)))/(1-kap_1m);
dif     = g0_pd-x;

end

