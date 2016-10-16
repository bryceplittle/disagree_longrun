function [LL] = log_likelihood(Y1,Y2,T,sT,params,kbar,tau,n_Y,n_shk,not0,jlag,e_1,H,A,B,D,g0_pd,g0_rf,gx_pd,gx_rf,Pi)

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

LL      = 0;
X_hat   = zeros(kbar+1,T+1);
Z       = [Y1, Y2];

% Z       = [Y1];
% not0    = not0*0;

try
    for tt = 1:T
        
        s       = sT(tt+1:tt+tau)';
        s_lag   = sT(tt:tt+tau-1)';
        jj      = binvec2dec(s)+1;
        jj_lag  = binvec2dec(s_lag)+1;
		
		if tt == 1
			BB  = B(:,:,jj_lag)*B(:,:,jj_lag)';
			P0  = reshape((eye((kbar+1)^2)-kron(A(:,:,2^tau),A(:,:,2^tau)))*BB(:),kbar+1,kbar+1);
		end
        
        % Construct mu_Z(s^tt), Q1(s^tt), Q2(s^tt), B_tilde(s^tt) %
        mu_Z    = [mu_c; mu_d; kap_0m+kap_1m*g0_pd(jj)-g0_pd(jj_lag)+mu_d; g0_rf(jj); ones(not0(tt),1)*mu_c];
        
        % Q1, Q2 for dC_{t+1|t} %
        Q1      = [zeros(2,kbar+1); kap_1m*gx_pd(:,:,jj)'; zeros(1,kbar+1); repmat(e_1'*H,not0(tt),1)];
        Q2      = [e_1'; rho_d*e_1'; rho_d*e_1'-gx_pd(:,:,jj_lag)'; gx_rf(:,:,jj_lag)'; zeros(not0(tt),kbar+1)];
        
        R       = [D(2:end,1:end-1,jj), zeros(n_Y-1,not0(tt)); zeros(not0(tt),n_shk), sqrt(e_1'*H*Pi(:,:,jj_lag)*H'*e_1)*eye(not0(tt))];
        B_tilde = [B(:,:,jj), zeros(kbar+1,not0(tt))];
        
        % Kalman Filter %
        Z_tilde     = Z(tt,1:4+not0(tt))'-mu_Z-(Q1*A(:,:,jj)+Q2)*X_hat(:,tt);
        Omega       = (Q1*A(:,:,jj)+Q2)*P0*(Q1*A(:,:,jj)+Q2)'+(Q1*B_tilde+R)*(Q1*B_tilde+R)';
		Omega_inv   = eye(size(Z_tilde,1))/Omega;
        K           = (A(:,:,jj)*P0*(Q1*A(:,:,jj)+Q2)'+B_tilde*B_tilde'*Q1'+B_tilde*R')*Omega_inv;
        X_hat(:,tt+1) = A(:,:,jj)*X_hat(:,tt) + K*Z_tilde;
        P1          = A(:,:,jj)*P0*A(:,:,jj)'+B_tilde*B_tilde';
        P0          = P1-K*Omega*K';
        
        % Likelihood
        LL = LL-0.5*(logdet(Omega)+Z_tilde'*Omega_inv*Z_tilde);
        
    end
    
   
    if isnan(LL) == 1;
        disp('LL is nan.');
        LL=-9e200;
    end
    if abs(LL) == Inf;
        disp('LL is inf.');
        LL=-9e200;
    end
    if isreal(LL) == 0;
        disp('LL not real.');
        LL=-9e200;
    end
    
catch
    disp('LL catch err.');
    LL = -9e200;
end
end

