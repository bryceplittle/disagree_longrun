function [gamma,psi,beta,pc,pd,kap_1,kap_0,kap_1m,kap_0m,mu_c,mu_d,rho_x,...
    rho_d,rho_cd,phi_x,phi_c,phi_d,phi_rm,phi_rf,phi_s,sig_h,pi,vscale] = unpack(Params)

gamma   = Params(1);
psi     = Params(2);
beta    = Params(3);
pc      = Params(4);
pd      = Params(5);
kap_1   = Params(6);
kap_0   = Params(7);
kap_1m  = Params(8);
kap_0m  = Params(9);
mu_c    = Params(10);
mu_d    = Params(11);
rho_x   = Params(12);
rho_d   = Params(13);
rho_cd  = Params(14);
phi_x   = Params(15);
phi_c   = Params(16);
phi_d   = Params(17);
phi_rm  = Params(18);
phi_rf  = Params(19);
phi_s   = Params(20);
sig_h   = Params(21);
pi      = Params(22);
vscale  = Params(23);

end

