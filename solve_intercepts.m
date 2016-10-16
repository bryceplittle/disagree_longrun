function [g0_pc,g0_pd,g0_rf,fail] = solve_intercepts(gx_pc,gx_pd,gx_rf,D,Pz,params,n_hist,jlead1,jlead0)

gamma       = params(1);
psi         = params(2);
beta        = params(3);
pc          = params(4);
pd          = params(5);
kap_1       = params(6);
kap_0       = params(7);
kap_1m      = params(8);
kap_0m      = params(9);
mu_c        = params(10);
mu_d        = params(11);
rho_x       = params(12);
rho_d       = params(13);
rho_cd      = params(14);
phi_x       = params(15);
phi_c       = params(16);
phi_d       = params(17);
phi_rm      = params(18);
phi_rf      = params(19);
phi_s       = params(20);
sig_h       = params(21);
pi          = params(22);
vscale      = params(23);
theta       = (1-gamma)/(1-1/psi);
dif         = inf;
dif_1       = inf;
dif_2       = inf;
dif_3       = inf;
count       = 1;
wt          = 0.9;
g0_pc       = ones(1,1,n_hist)*params(4);
g0_pd       = ones(1,1,n_hist)*params(5);
g0_rf       = zeros(1,1,n_hist);
g0_pc_old   = zeros(1,1,n_hist);
g0_pd_old   = zeros(1,1,n_hist);
g0_rf_old   = zeros(1,1,n_hist);
pc_cvar1    = 0;
pc_cvar0    = 0;
pd_cvar1    = 0;
pd_cvar0    = 0;
rf_cvar1    = 0;
rf_cvar0    = 0;

while dif>1e-6 && count<1e5
    for j=1:n_hist
        pc_cvar1 = [theta*kap_1*gx_pc(:,:,jlead1(j)+1);0;theta*(1-1/psi)*D(2,2,jlead1(j)+1);0;0;0]'*Pz(:,:,jlead1(j)+1)*[theta*kap_1*gx_pc(:,:,jlead1(j)+1);0;theta*(1-1/psi)*D(2,2,jlead1(j)+1);0;0;0];
        pc_cvar0 = [theta*kap_1*gx_pc(:,:,jlead0(j)+1);0;theta*(1-1/psi)*D(2,2,jlead0(j)+1);0;0;0]'*Pz(:,:,jlead0(j)+1)*[theta*kap_1*gx_pc(:,:,jlead0(j)+1);0;theta*(1-1/psi)*D(2,2,jlead0(j)+1);0;0;0];
        g0_pc(1,1,j) = log(beta)+kap_0+(1-1/psi)*mu_c+kap_1*(pi*g0_pc(1,1,jlead1(j)+1)+(1-pi)*g0_pc(1,1,jlead0(j)+1))+...
            (pi/(2*theta))*pc_cvar1+((1-pi)/(2*theta))*pc_cvar0;
        
        rf_cvar1 = [(theta-1)*kap_1*gx_pc(:,:,jlead1(j)+1); 0; (theta-theta/psi-1)*D(2,2,jlead1(j)+1);0;0;phi_rf]'*Pz(:,:,jlead1(j)+1)*[(theta-1)*kap_1*gx_pc(:,:,jlead1(j)+1); 0; (theta-theta/psi-1)*D(2,2,jlead1(j)+1);0;0;phi_rf];
        rf_cvar0 = [(theta-1)*kap_1*gx_pc(:,:,jlead0(j)+1); 0; (theta-theta/psi-1)*D(2,2,jlead0(j)+1);0;0;phi_rf]'*Pz(:,:,jlead0(j)+1)*[(theta-1)*kap_1*gx_pc(:,:,jlead0(j)+1); 0; (theta-theta/psi-1)*D(2,2,jlead0(j)+1);0;0;phi_rf];
        g0_rf(1,1,j) = -log(beta)-(theta-1)*kap_0/theta-(theta-theta/psi-1)*mu_c/theta-...
            (theta-1)*(kap_1*(pi*g0_pc(1,1,jlead1(j)+1)+(1-pi)*g0_pc(1,1,jlead0(j)+1))-g0_pc(1,1,j))/theta+...
            (pi/2)*rf_cvar1/theta+((1-pi)/2)*rf_cvar0/theta;
        
        pd_cvar1 = [(theta-1)*kap_1*gx_pc(:,:,jlead1(j)+1)+kap_1m*gx_pd(:,:,jlead1(j)+1);0;(theta-theta/psi-1-rho_cd)*D(2,2,jlead1(j)+1); D(3,2,jlead1(j)+1); phi_rm; 0]'*Pz(:,:,jlead1(j)+1)*...
            [(theta-1)*kap_1*gx_pc(:,:,jlead1(j)+1)+kap_1m*gx_pd(:,:,jlead1(j)+1);0;(theta-theta/psi-1-rho_cd)*phi_c*D(2,2,jlead1(j)+1); D(3,2,jlead1(j)+1); phi_rm; 0];
        pd_cvar0 = [(theta-1)*kap_1*gx_pc(:,:,jlead0(j)+1)+kap_1m*gx_pd(:,:,jlead0(j)+1);0;(theta-theta/psi-1-rho_cd)*phi_c*D(2,2,jlead0(j)+1); D(3,2,jlead0(j)+1); phi_rm; 0]'*Pz(:,:,jlead0(j)+1)*...
            [(theta-1)*kap_1*gx_pc(:,:,jlead0(j)+1)+kap_1m*gx_pd(:,:,jlead0(j)+1);0;(theta-theta/psi-1-rho_cd)*phi_c*D(2,2,jlead0(j)+1); D(3,2,jlead0(j)+1); phi_rm; 0];
        g0_pd(1,1,j) = theta*log(beta)+(theta-1-theta/psi)*mu_c+(theta-1)*kap_0+kap_0m+mu_d+...
            kap_1m*(pi*g0_pd(1,1,jlead1(j)+1)+(1-pi)*g0_pd(1,1,jlead0(j)+1))+...
            (theta-1)*(kap_1*(pi*g0_pc(1,1,jlead1(j)+1)+(1-pi)*g0_pc(1,1,jlead0(j)+1))-g0_pc(1,1,j))+...
            (pi/2)*pd_cvar1+((1-pi)/2)*pd_cvar0;
    end
    
    dif_1       = max(max(max(abs(g0_pc-g0_pc_old))));
    dif_2       = max(max(max(abs(g0_pd-g0_pd_old))));
    dif_3       = max(max(max(abs(g0_rf-g0_rf_old))));
    dif         = max([dif_1,dif_2,dif_3]);
    g0_pc_old   = wt*g0_pc+(1-wt)*g0_pc_old;
    g0_rf_old   = wt*g0_rf+(1-wt)*g0_rf_old;
    g0_pd_old   = wt*g0_pd+(1-wt)*g0_pd_old;
    
    count       = count+1;
end

if count>1e5
    fail = 1;
else
    fail = 0;
end

end

