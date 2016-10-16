function [A,B,C1,C2,D,g0_pc,g0_pd,g0_rf,gx_pc,gx_pd,gx_rf,Pi,params,fail] = solve_model(params,kbar,n_hist,jlag,jlead1,jlead0)

% clear solve_AB_mex solve_condvcov_mex solve_intercepts_mex solve_xvar_mex;
% clear mex;
% warning('off','Coder:MATLAB:singularMatrix');

% setup;

n_shk       = 5;                                        % num aggregate shocks
n_Y         = 5;                                        % num observables
tol         = 1e-6;                                     % convergence criteria
fail        = 0;                                        % model solution failure boolean

%% Compute F.I.R.E. Steady State %%

try
    options = optimset('MaxIter',1e20,'MaxFunEvals',1e20,'TolFun',1e-6,'Display','off');
    
    Fpc_bar_dif = @(x) pc_bar_dif(x,params);
    pc_bar      = fsolve(Fpc_bar_dif,5,options);
    
    if abs(pc_bar) < inf && isnan(pc_bar) == 0
        params(4) = pc_bar;
        params(6) = exp(pc_bar)/(1+exp(pc_bar));
        params(7) = log(1+exp(pc_bar))-params(6)*pc_bar;
    else
        fail = 1;
    end
    
    Fpd_bar_dif = @(x) pd_bar_dif(x,params);
    pd_bar      = fsolve(Fpd_bar_dif,5,options);
    
    if abs(pd_bar) < inf && isnan(pd_bar) == 0
        params(5) = pd_bar;
        params(8) = exp(pd_bar)/(1+exp(pd_bar));
		if params(8) > 0.999
			params(8) = 0.999;
		end
        params(9) = log(1+exp(pd_bar))-params(8)*pd_bar;
    else
        fail = 1;
    end
catch
    fail = 1;
end

%% Model Solution %%

try
    [A,B,D,P,p] = first_guess(params,n_hist,n_shk,n_Y,kbar); % FIRE soln
    [A,B,C1,C2,K,gx_pc,gx_pd,gx_rf,fail] = solve_AB(params,A,B,D,P,p,kbar,n_shk,n_Y,n_hist,jlead1,jlead0,jlag,tol,fail);
	% [A,B,C1,C2,K,gx_pc,gx_pd,gx_rf,fail] = solve_AB_mex(params,A,B,D,P,p,kbar,n_shk,n_Y,n_hist,jlead1,jlead0,jlag,tol,fail);
    
    if fail == 0
		% optional mex functions written in c++
        % [Pz] = solve_condvcov_mex(A,B,C1,C2,D,n_shk,n_Y,kbar,n_hist,jlag);
        % [g0_pc,g0_pd,g0_rf,fail] = solve_intercepts_mex(gx_pc,gx_pd,gx_rf,D,Pz,params,n_hist,jlead1,jlead0);
        % [Pi] = solve_xvar_mex(A,C1,C2,D,K,n_hist,jlag);
        
        [Pz] = solve_condvcov(A,B,C1,C2,D,n_shk,n_Y,kbar,n_hist,jlag);  % conditional variances
        [g0_pc,g0_pd,g0_rf,fail] = solve_intercepts(gx_pc,gx_pd,gx_rf,D,Pz,params,n_hist,jlead1,jlead0);  % intercepts
        [Pi] = solve_xvar(A,C1,C2,D,K,n_hist,jlag);  % cross sectional variance
    end;
    if fail == 1
        A       = nan;
        B       = nan;
        C1      = nan;
        C2      = nan;
        K       = nan;
        gx_pc   = nan;
        gx_pd   = nan;
        gx_rf   = nan;
        g0_pc   = nan;
        g0_pd   = nan;
        g0_rf   = nan;
        gx_pc   = nan;
        gx_pd   = nan;
        gx_rf   = nan;
        Pi      = nan;
    end
catch
    fail = 1;
end

end
