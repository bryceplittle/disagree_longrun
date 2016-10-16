function [A,B,C1,C2,K,gx_pc,gx_pd,gx_rf,fail] = solve_AB(params,A,B,D,P,p,kbar,n_shk,n_Y,n_hist,jlead1,jlead0,jlag,tol,fail,A_old,B_old,C1,C2,K,K_old,L,gx_pc,gx_pd,gx_rf,gx_pc_old,gx_pd_old,gx_rf_old,A1,A2,A2_big,A3,A3_big,B1,B2,B2_big,e_1,H)

%% Parameters %%

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
nan_1   = 0;
nan_2   = 0;
nan_3   = 0;
nan_4   = 0;
nan_5   = 0;
nan_6   = 0;
dif     = inf;
dif_1   = inf;
dif_2   = inf;
dif_3   = inf;
dif_4   = inf;
dif_5   = inf;
dif_6   = inf;
count   = 1;


%% Main Loop %%

while dif>tol && fail==0
    for jj=1:n_hist
        
        %% Coefficients %%
        
        gx_pc(:,:,jj)  = ((1-1/psi)*e_1'*H+kap_1*(pi*gx_pc(:,:,jlead1(jj)+1)'*A(:,:,jlead1(jj)+1)+(1-pi)*gx_pc(:,:,jlead0(jj)+1)'*A(:,:,jlead0(jj)+1))*H)';
        gx_pd(:,:,jj)  = (theta-1-theta/psi+rho_d)*e_1'*H+...
            (theta-1)*(kap_1*(pi*gx_pc(:,:,jlead1(jj)+1)'*A(:,:,jlead1(jj)+1)+(1-pi)*gx_pc(:,:,jlead0(jj)+1)'*A(:,:,jlead0(jj)+1))*H-gx_pc(:,:,jj)')+...
            kap_1m*(pi*gx_pd(:,:,jlead1(jj)+1)'*A(:,:,jlead1(jj)+1)+(1-pi)*gx_pd(:,:,jlead0(jj)+1)'*A(:,:,jlead0(jj)+1))*H;
        gx_rf(:,:,jj)  = -(theta-theta/psi-1)*e_1'*H-(theta-1)*(kap_1*(pi*gx_pc(:,:,jlead1(jj)+1)'*A(:,:,jlead1(jj)+1)+(1-pi)*gx_pc(:,:,jlead0(jj)+1)'*A(:,:,jlead0(jj)+1))*H-gx_pc(:,:,jj)');
        
        %% Kalman Filter %%
        
        C1(:,:,jj)   = [e_1'; zeros(2,kbar+1); kap_1m*gx_pd(:,:,jj)'; zeros(1,kbar+1)];
        C2(:,:,jj)   = [zeros(1,kbar+1); e_1'; rho_d*e_1'; rho_d*e_1'-gx_pd(:,:,jlag(jj)+1)'; gx_rf(:,:,jlag(jj)+1)'];
        
        P(:,:,jj)    = A(:,:,jj)*p(:,:,jlag(jj)+1)*A(:,:,jj)'+B(:,:,jj)*B(:,:,jj)';
        L(:,:,jj)    = ((C1(:,:,jj)*A(:,:,jj)+C2(:,:,jj))*p(:,:,jlag(jj)+1)*(C1(:,:,jj)*A(:,:,jj)+C2(:,:,jj))'+(C1(:,:,jj)*B(:,:,jj)+D(:,1:n_shk,jj))*(C1(:,:,jj)*B(:,:,jj)+D(:,1:n_shk,jj))'+D(:,n_shk+1:end,jj)*D(:,n_shk+1:end,jj)');
        if sum(sum(isnan(L(:,:,jj)))) == 0 || max(max(abs(L(:,:,jj))))<inf
            L_inv    = pinv(L(:,:,jj));
        else
            L_inv    = nan(size(L(:,:,jj)));
        end
        K(:,:,jj)    = (A(:,:,jj)*p(:,:,jlag(jj)+1)*(C1(:,:,jj)*A(:,:,jj)+C2(:,:,jj))'+B(:,:,jj)*B(:,:,jj)'*C1(:,:,jj)'+B(:,:,jj)*D(:,1:n_shk,jj)')*L_inv;
        p(:,:,jj)    = P(:,:,jj)-K(:,:,jj)*L(:,:,jj)*K(:,:,jj)';
        
        %% Law of Motion %%
        
        A1        = diag(e_1*rho_x);
        A2_big    = [zeros(1,kbar+2); zeros(kbar+1,1), A(:,:,jj)-K(:,:,jj)*(C1(:,:,jj)*A(:,:,jj)+C2(:,:,jj))];
        A2        = A2_big(1:end-1,1:end-1);
        A3_big    = [zeros(1,kbar+1); K(:,:,jj)*(C1(:,:,jj)*A(:,:,jj)+C2(:,:,jj))];
        A3        = A3_big(1:end-1,:);
        A(:,:,jj) = A1+A2+A3;
        
        
        B2_big    = [zeros(1,n_shk); K(:,:,jj)*(C1(:,:,jj)*B(:,:,jj)+D(:,1:n_shk,jj))];
        B2        = B2_big(1:end-1,:);
        B(:,:,jj) = B1(:,:,jj)+B2;
        
    end
    
    %% Compute Step Size %%
    
    dif_1   = max(max(max(abs(A-A_old))));
    dif_2   = max(max(max(abs(B-B_old))));
    dif_3   = max(max(max(abs(K-K_old))));
    dif_4   = max(max(max(abs(gx_pc-gx_pc_old))));
    dif_5   = max(max(max(abs(gx_pd-gx_pd_old))));
    dif_6   = max(max(max(abs(gx_rf-gx_rf_old))));
    dif     = max([dif_1,dif_2,dif_3,dif_4,dif_5,dif_6]);
    
    count   = count+1;
    
    %% Check Model Solution %%
    
    nan_1   = sum(sum(sum(isnan(A))));
    nan_2   = sum(sum(sum(isnan(B))));
    nan_3   = sum(sum(sum(isnan(K))));
    nan_4   = sum(sum(sum(isnan(gx_pc))));
    nan_5   = sum(sum(sum(isnan(gx_pd))));
    nan_6   = sum(sum(sum(isnan(gx_rf))));
    
    if nan_1+nan_2+nan_3+nan_4+nan_5+nan_6>0 || count>1e4
        fail = 1;
    end
    
    if fail == 1
        dif     = 0;
        count   = inf;
    end
    
    %% Store Old Values %%
    
    A_old       = A;
    B_old       = B;
    K_old       = K;
    gx_pc_old   = gx_pc;
    gx_pd_old   = gx_pd;
    gx_rf_old   = gx_rf;
    
end
end

