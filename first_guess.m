function [A,B,D,P,p] = first_guess(Params,n_hist,n_shk,n_Y,kbar)

[gamma,psi,beta,pc,pd,kap_1,kap_0,kap_1m,kap_0m,mu_c,mu_d,rho_x,...
    rho_d,rho_cd,phi_x,phi_c,phi_d,phi_rm,phi_rf,phi_s,sig_h,pi,phi_ce] = unpack(Params);
theta   = (1-gamma)/(1-1/psi);

%% Pre-Allocate Arrays %%

A       = zeros(kbar+1,kbar+1,n_hist);
B       = zeros(kbar+1,n_shk,n_hist);
P       = zeros(kbar+1,kbar+1,n_hist);
p       = zeros(kbar+1,kbar+1,n_hist);

%% Measurement Equation Covariance Matrix %%

D0              = zeros(n_Y,n_shk+1);
D0(1,n_shk+1)   = phi_s;
D0(2,2)         = phi_c;
D0(3,3)         = phi_d;
D0(3,2)         = rho_cd*phi_c;
D0(4,4)         = phi_rm;
D0(4,3)         = phi_d;
D0(4,2)         = rho_cd*phi_c;
D0(5,5)         = phi_rf;

for j = 1:n_hist
    A(:,:,j)        = eye(kbar+1,kbar+1)*rho_x*0.999;
    B(:,:,j)        = [ones(kbar+1,1)*phi_x, zeros(kbar+1,n_shk-1)];
    D(:,:,j)        = D0;
    P(:,:,j)        = B(:,:,j)*B(:,:,j)';
    p(:,:,j)        = B(:,:,j)*B(:,:,j)';
    
end

for j = 2:2:n_hist
    B(1,1,j)  = B(1,1,j)*sig_h;
    D(2,2,j)  = D(2,2,j)*sig_h;
    D(3,3,j)  = D(3,3,j)*sig_h;
    D(3,2,j)  = D(3,2,j)*sig_h;
    D(4,3,j)  = D(4,3,j)*sig_h;
    D(4,2,j)  = D(4,2,j)*sig_h;
end

end

