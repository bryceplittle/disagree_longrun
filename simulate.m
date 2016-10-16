clear, clc;
% mainpath = 'C:/Users/bryce/Dropbox/Research/macrofinance/dynhoe/code/main';
mainpath = '/Users/brycelittle/Dropbox/research/macrofinance/dynhoe/code/main';
addpath(mainpath); cd(mainpath); setup;

[A,B,D,P,p] = ...
    first_guess(params,n_hist,n_shk,n_Y,kbar);
[params,fail] = ...
    solve_steady(params);
[A,B,C1,C2,K,gx_pc,gx_pd,gx_rf,fail] = ...
    solve_AB(params,A,B,D,P,p,kbar,n_shk,n_Y,n_hist,jlead1,jlead0,jlag,tol,fail,A_old,B_old,C1,C2,K,K_old,L,gx_pc,gx_pd,gx_rf,gx_pc_old,gx_pd_old,gx_rf_old,A1,A2,A2_big,A3,A3_big,B1,B2,B2_big,e_1,H);
[Pz] = ...
    solve_condvcov(A,B,C1,C2,D,n_shk,n_Y,kbar,n_hist,jlag);
[g0_pc,g0_pd,g0_rf,fail] = ...
    solve_intercepts(gx_pc,gx_pd,gx_rf,D,Pz,params,n_hist,jlead1,jlead0);
[Pi] = ...
    solve_xvar(A,C1,C2,D,K,n_hist,jlag);

gamma   = params(1);
psi     = params(2);
beta    = params(3);
pc      = params(4);
pd      = params(5);
kap_1   = params(6);
kap_0   = params(7);
kap_1m  = params(8);
kap_0m  = params(9);
mu_c    = params(10);
mu_d    = params(11);
rho_x   = params(12);
rho_d   = params(13);
rho_cd  = params(14);
phi_x   = params(15);
phi_c   = params(16);
phi_d   = params(17);
phi_rm  = params(18);
phi_rf  = params(19);
phi_s   = params(20);
sig_h   = params(21);
pi      = params(22);
phi_ce  = params(23);
theta   = (1-gamma)/(1-1/psi);

% T 		= 40; % 1e6;
% sT 		= zeros(T+tau,1);
% eps 	= zeros(T,n_shk);
% Y 		= zeros(T,n_Y-1);
% X 		= zeros(T,kbar+1);
% ret     = zeros(T,1);

%% output kalman gain %%
% xlswrite('kalman_gain_lov.xlsx', K(:,:,1));
% xlswrite('kalman_gain_hiv.xlsx', K(:,:,256));

%% simulate cross-sectional variance %%

T       = 100+tau;
X	    = zeros(T,kbar+1);
n_F     = 50;
Fc      = zeros(T,n_F);
sT      = zeros(T,1); sT(tau+50:end) = 1;
eps     = zeros(T,n_shk);

for tt = tau+1:T
    eps(tt,:) = mvnrnd(zeros(n_shk,1),eye(n_shk));
end

for tt = tau+1:T
    s           = sT(tt-tau+1:tt)';
    s_lag       = sT(tt-tau:tt-1)';
    jj          = binvec2dec(s)+1;
    jj_lag      = binvec2dec(s_lag)+1;
    X(tt,:)     = (A(:,:,jj)*X(tt-1,:)'+B(:,:,jj)*eps(tt,:)')';
end

for tt = tau+1:T
    s           = sT(tt-tau+1:tt)';
    s_lag       = sT(tt-tau:tt-1)';
    jj          = binvec2dec(s)+1;
    jj_lag      = binvec2dec(s_lag)+1;
    
    for ff = 1:n_F
        Fc(tt,ff) = mu_c + e_1'*H*X(tt,:)' + mvnrnd(0,1)*sqrt(e_1'*H*Pi(:,:,jj)*H'*e_1);
    end
end

% xlswrite('XVAR_X_sim.xlsx', X);
% xlswrite('XVAR_Fc_sim.xlsx',Fc);

%% irfs %%

T 		= 100; % 1e6;
X_hiv	= zeros(T,kbar+1);
X_lov	= zeros(T,kbar+1);
pd_hiv  = zeros(T,1);
pd_lov  = zeros(T,1);
pd_fi   = zeros(T,1);
pc_hiv  = zeros(T,1);
pc_lov  = zeros(T,1);
pc_filo = zeros(T,1);
pc_fihi = zeros(T,1);
pd_filo = zeros(T,1);
pd_fihi = zeros(T,1);
A1M     = (rho_d-1/psi)/(1-kap_1m*rho_x);
A1      = (1-1/psi)/(1-kap_1*rho_x);

for tt = 2:T
    if tt == 2
        X_hiv(tt,:) = A(:,:,256)*X_hiv(tt-1,:)' + B(:,:,256)*[2,0,0,0,0]';
        X_lov(tt,:) = A(:,:,1)*X_lov(tt-1,:)' + B(:,:,1)*[2,0,0,0,0]';
    else
        X_hiv(tt,:) = A(:,:,256)*X_hiv(tt-1,:)';
        X_lov(tt,:) = A(:,:,1)*X_lov(tt-1,:)';
    end
    
    pd_hiv(tt)  = gx_pd(:,:,256)'*X_hiv(tt,:)';
    pd_lov(tt)  = gx_pd(:,:,1)'*X_lov(tt,:)';
    pc_hiv(tt)  = gx_pc(:,:,256)'*X_hiv(tt,:)';
    pc_lov(tt)  = gx_pc(:,:,1)'*X_lov(tt,:)';
    pc_filo(tt) = A1*X_lov(tt,1);
    pc_fihi(tt) = A1*X_hiv(tt,1);
    pd_filo(tt) = A1M*X_lov(tt,1);
    pd_fihi(tt) = A1M*X_hiv(tt,1);
end

% plot(X_hiv)
% plot(X_lov)

plot(exp(pc_hiv))
hold on
plot(exp(pc_lov))
plot(exp(pc_filo))
plot(exp(pc_fihi))

plot(exp(pd_hiv))
hold on
plot(exp(pd_lov))
plot(exp(pd_filo))
plot(exp(pd_fihi))

xlswrite('IRF_X_hiv_rfshock.xlsx', X_hiv);
xlswrite('IRF_X_lov_rfshock.xlsx', X_lov);

xlswrite('IRF_pc_hiv.xlsx', pc_hiv);
xlswrite('IRF_pc_lov.xlsx', pc_lov);
xlswrite('IRF_pc_filo.xlsx', pc_filo);
xlswrite('IRF_pc_fihi.xlsx', pc_fihi);

xlswrite('IRF_pd_hiv.xlsx', pd_hiv);
xlswrite('IRF_pd_lov.xlsx', pd_lov);
xlswrite('IRF_pd_filo.xlsx', pd_filo);
xlswrite('IRF_pd_fihi.xlsx', pd_fihi);

% simulate history of vol states and shocks
for tt = tau+1:T
    if rand < params(22)
        sT(tt) = 1;
    end
    eps(tt,:) = mvnrnd(zeros(n_shk,1),eye(n_shk));
end;


%% simulate state X and observables Y %%
for tt = tau+1:T
    s           = sT(tt-tau+1:tt)';
    s_lag       = sT(tt-tau:tt-1)';
    jj          = binvec2dec(s)+1;
    jj_lag      = binvec2dec(s_lag)+1;
    X(tt,:)     = (A(:,:,jj)*X(tt-1,:)'+B(:,:,jj)*eps(tt,:)')';
    mu_Y        = [mu_c; mu_d; kap_0m+kap_1m*g0_pd(jj)-g0_pd(jj_lag)+mu_d; g0_rf(jj)];
    Y(tt,:)     = (mu_Y+C1(2:end,:,jj)*X(tt,:)'+C2(2:end,:,jj_lag)*X(tt-1,:)'+D(2:end,1:end-1,jj)*eps(tt,:)')';
end;


%% compute risk prices %%
lam_e = zeros(n_shk,1,n_hist);
lam_x = zeros(kbar+1,1,n_hist);
lam_s = zeros(1,1,n_hist);

for jj = 1:n_hist
    lam_e(:,:,jj) = -gamma*D(2,2,jj)*e_c'+(theta-1)*kap_1*gx_pc(:,:,jj)'*B(:,:,jj);
    lam_e(1,1,jj) = lam_e(1,1,jj)/B(1,1,jj);
    lam_e(2,1,jj) = lam_e(2,1,jj)/D(2,2,jj);
    lam_e(3,1,jj) = lam_e(3,1,jj)/D(3,3,jj);
    lam_e(4,1,jj) = lam_e(4,1,jj)/D(4,4,jj);
    lam_e(5,1,jj) = lam_e(5,1,jj)/D(5,5,jj);
    lam_x(:,:,jj) = -gamma*(e_1'-e_1'*H)+(theta-1)*kap_1*(gx_pc(:,:,jj)'*A(:,:,jj)-...
        (pi*gx_pc(:,:,jlead1(jlag(jj)+1)+1)'*A(:,:,jlead1(jlag(jj)+1)+1)*H+...
        (1-pi)*gx_pc(:,:,jlead0(jlag(jj)+1)+1)'*A(:,:,jlead0(jlag(jj)+1)+1)*H));
    lam_s(:,:,jj) = (theta-1)*kap_1*(g0_pc(:,:,jj)-pi*g0_pc(:,:,jlead1(jlag(jj)+1)+1)-(1-pi)*g0_pc(:,:,jlead0(jlag(jj)+1)+1));
end


%% observe variation in risk prices %%
s_mat   = [triu(ones(tau, tau)); zeros(1,tau)];
riskp   = [];

for ii = 1:tau+1
    s       = s_mat(ii,:);
    jj      = binvec2dec(s)+1;
    riskp   = [lam_e(1,:,jj), riskp];
end

bar(riskp)

%% compute betas %%
beta_ex_fi  = -kap_1m*(rho_d-1/psi)/(1-kap_1m);
beta_e      = zeros(n_shk,1,n_hist);

for jj = 1:n_hist
    beta_e(:,:,jj) = D(3,3,jj)*e_d'+kap_1m*gx_pd(:,:,jj)'*B(:,:,jj);
    beta_e(1,1,jj) = beta_e(1,1,jj)/B(1,1,jj);
    beta_e(2,1,jj) = beta_e(2,1,jj)/D(2,2,jj);
    beta_e(3,1,jj) = beta_e(3,1,jj)/D(3,3,jj);
    beta_e(4,1,jj) = beta_e(4,1,jj)/D(4,4,jj);
    beta_e(5,1,jj) = beta_e(5,1,jj)/D(5,5,jj);
end


%% compute average risk price %%
% 0*phi : -363.2985
% 1*phi : -187.2069
% 2*phi : -62.3483
% 3*phi : -24.4617
% 4*phi : -12.7712
lam_e_bar = zeros(n_shk,1);
lam_ex_fi = ((gamma-1/psi)*(kap_1))/(1-kap_1*rho_x);

for jj = 1:n_hist
    s_hist(:,:,jj) = dec2binvec((jj-1), tau);
    lam_e_bar = lam_e_bar + bern_pmf(tau,pi,sum(s_hist(:,:,jj)==1))*lam_e(:,:,jj);
end


%% compute the average beta for LRR shock %%
% 0*phi : 15.0480
% 1*phi : 8.0511
% 2*phi : 3.6916
% 3*phi : 2.0475
% 4*phi : 1.3019
beta_e_bar = zeros(n_shk,1);
for jj = 1:n_hist
    s_hist(:,:,jj) = dec2binvec((jj-1), tau);
    beta_e_bar = beta_e_bar + bern_pmf(tau,pi,sum(s_hist(:,:,jj)==1))*beta_e(:,:,jj);
end


%% output pca loadings to excel %%
% [coef, score, latent] = princomp(X);
% xlswrite('X_pca_coef.xlsx',coef);
% xlswrite('X_pca_score.xlsx',score);
% xlswrite('X_sim.xlsx',X);
% round(latent./sum(latent),3)


%% output lrr risk prices to excel %%
% lam_e_out = [];
% for jj = 1:n_hist
%     lam_e_out(jj,:) = lam_e(:,:,jj)';
% end
%
% xlswrite('lam_e_output.xlsx', lam_e_out);


