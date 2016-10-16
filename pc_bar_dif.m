function [dif] = pc_bar_dif(x,Params)

[gamma,psi,beta,pc,pd,kap_1,kap_0,kap_1m,kap_0m,mu_c,mu_d,rho_x,...
    rho_d,rho_cd,phi_x,phi_c,phi_d,phi_rm,phi_rf,phi_s,sig_h,pi] = unpack(Params);
theta = (1-gamma)/(1-1/psi);

kap_1   = exp(x)/(1+exp(x));
kap_0   = log(1+exp(x))-kap_1*x;
gx_pc   = (1-1/psi)/(1-kap_1*rho_x);
g0_pc   = (log(beta)+kap_0+(1-1/psi)*mu_c+(1/2)*(1/theta)*((theta*(1-1/psi)*phi_c)^2+(kap_1*gx_pc*phi_x)^2)*(pi*sig_h^2+(1-pi)))/(1-kap_1);
dif     = g0_pc-x; 

end

