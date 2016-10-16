kbar        = 12;                                       % order of expectations
tau         = 8;                                        % number of periods to track for vol history
n_shk       = 5;                                        % number of aggregate shocks
n_Y         = 5;                                        % number of obs variables
tol         = 1e-10;                                    % model solution tolerance
n_hist      = 2^tau;                                    % total number of histories
bbase       = zeros(1,tau);                             % binary base
e_1	 	    = [1, zeros(1,kbar)]';
temp		= eye(n_shk); e_c = temp(:,2); e_d = temp(:,3); clear temp;
H           = [zeros(kbar,1), eye(kbar); zeros(1,kbar+1)];
mainpath 	= 'C:/Users/bryce/Dropbox/Research/macrofinance/dynhoe/code/main';

for j=1:tau;
    bbase(1,tau-j+1) = 2^(j-1);
end

jlag    = zeros(n_hist,1);                              % index for a given history
jlead1  = zeros(n_hist,1);                              % index for lead into vol_hi given a history
jlead0  = zeros(n_hist,1);                              % index for lead into vol_lo given a history

for j=1:n_hist;
    zvec        = dec2binvec(j-1,tau);                  % volatility history index
    jlag(j)     = bbase*[0 zvec(1,1:end-1)]';
    jlead1(j)   = bbase*[zvec(1,2:end) 1]';
    jlead0(j)   = bbase*[zvec(1,2:end) 0]';
end

%% Data Prep %%

% import_Y1;
% import_Y2;
% Y1      = Y1(:,2:end);
% Y2      = Y2(:,2:end);
% 
% import_Y1_qtr;
% import_Y2_qtr;
% caltime = Y1(1:end,1:2);
% Y1      = Y1(1:end,3:end);
% Y2      = Y2(1:end,3:end);
% 
% T       = size(Y1,1);
% not0    = sum(isnan(Y2) == 0,2);
% not0    = not0.*(not0 > 1);

%% Initialize Parameters %%

[params]                = initial_params();
[params_UB, params_LB]  = param_bounds();

A         = zeros(kbar+1,kbar+1,n_hist);
A_old     = zeros(kbar+1,kbar+1,n_hist);
B         = zeros(kbar+1,n_shk,n_hist);
B_old     = zeros(kbar+1,n_shk,n_hist);
C1        = zeros(n_Y,kbar+1,n_hist);
C2        = zeros(n_Y,kbar+1,n_hist);
D         = zeros(n_Y,n_shk+1);
P         = zeros(kbar+1,kbar+1,n_hist);
p         = zeros(kbar+1,kbar+1,n_hist);
K         = zeros(kbar+1,n_Y,n_hist);
K_old     = zeros(kbar+1,n_Y,n_hist);
L         = zeros(n_Y,n_Y,n_hist);
gx_pc     = zeros(kbar+1,1,n_hist);
gx_pd     = zeros(kbar+1,1,n_hist);
gx_rf     = zeros(kbar+1,1,n_hist);
gx_pc_old = zeros(kbar+1,1,n_hist);
gx_pd_old = zeros(kbar+1,1,n_hist);
gx_rf_old = zeros(kbar+1,1,n_hist);
A1        = zeros(kbar+1,kbar+1);
A2        = zeros(kbar+1,kbar+1);
A2_big    = zeros(kbar+2,kbar+2);
A3        = zeros(kbar+1,kbar+1);
A3_big    = zeros(kbar+2,kbar+2);
B1        = zeros(kbar+1,n_shk,n_hist);
B2        = zeros(kbar+1,n_shk);
B2_big    = zeros(kbar+2,n_shk);
e_1       = [1, zeros(1,kbar)]';
H         = [zeros(kbar,1), eye(kbar); zeros(1,kbar+1)];

for jj = 1:n_hist
    if mod(jj,2) == 0
        B1(:,:,jj) = [e_1, zeros(kbar+1,n_shk-1)]*params(15)*params(21);
    else
        B1(:,:,jj) = [e_1, zeros(kbar+1,n_shk-1)]*params(15);
    end
end
