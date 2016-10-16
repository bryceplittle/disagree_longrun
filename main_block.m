clear, clc;
mainpath = 'C:/Users/bryce/Dropbox/Research/macrofinance/dynhoe/code/main';
addpath(mainpath); cd(mainpath); setup;

jj1_max         = 25000;  							% size of MCMC sample
n_block         = 2;  								% number of parameter blocks in sampling
n_p             = size(params,1);  					% number of parameters
block_draw      = zeros(n_p,1);  					% stores most recent vector in markov chain
prob_move       = 4/(tau*size(Y1,1));  				% the probability of flipping a volatility state in a proposal
prop_cov        = diag(abs(params))*1e-5;  			% scales proposal covariance matrix

% sets certain parameters to constants
prop_cov(10,10) = 0;  % AR(1) intercept
prop_cov(11,11) = 0;  % AR(1) intercept
prop_cov(23,23) = 0;  % unused parameter

restart 		= 0;

if restart == 1
    
    load('workspace');
	jj1_max         = 25000;
	prob_move       = 4/(tau*size(Y1,1));
    
else
    
    sT = [zeros(tau,1); zeros(T,1)];
    for t = 1:T+tau
        if log(rand) <= log(params(22))
            sT(t) = 1;
        end
    end
    
    pdraw           = params;
    p_max           = params;
    params_chain    = [];
    sT_chain        = [];
    sT_max          = sT;
    sT_cand         = sT;
    ll_pcand        = -9e+200;
    ll_pdraw        = -9e+200;
    pdraw_prior     = -9e+200;
    pcand_prior     = -9e+200;
    lpost_max       = -9e+200;
    accept_count    = 0;
    accept_count_sT = 0;
    outside_count   = 0;
    ll_fail_count 	= 0;
    count           = 0;
    fail_count      = 0;
    
    [pdraw,fail] = solve_steady(params);
    [A_draw,B_draw,D,P_draw,p_draw] = first_guess(pdraw,n_hist,n_shk,n_Y,kbar);
    [A_draw,B_draw,C1_draw,C2_draw,K_draw,gx_pc_draw,gx_pd_draw,gx_rf_draw,fail] = solve_AB(pdraw,A_draw,B_draw,D,P_draw,p_draw,kbar,n_shk,n_Y,n_hist,jlead1,jlead0,jlag,tol,fail,A_old,B_old,C1,C2,K,K_old,L,gx_pc,gx_pd,gx_rf,gx_pc_old,gx_pd_old,gx_rf_old,A1,A2,A2_big,A3,A3_big,B1,B2,B2_big,e_1,H);
    [Pz_draw] = solve_condvcov(A_draw,B_draw,C1_draw,C2_draw,D,n_shk,n_Y,kbar,n_hist,jlag);
    [g0_pc_draw,g0_pd_draw,g0_rf_draw,fail] = solve_intercepts(gx_pc_draw,gx_pd_draw,gx_rf_draw,D,Pz_draw,pdraw,n_hist,jlead1,jlead0);
    [Pi_draw] = solve_xvar(A_draw,C1_draw,C2_draw,D,K_draw,n_hist,jlag);
    
end

tic;
for jj1 = 1:jj1_max
    
    count = count+1;
    
    % randomly generate two blocks
    for pp = 1:n_p
        block_draw(pp) = rand;
    end
    [draw_sorted, params_order] = sort(block_draw);
    block{1} = params_order <= 11;
    block{2} = params_order >  11;
    
    %%% sample parameter blocks 1 and 2 %%%%
    
    for bb = 1:n_block
        
        clear mex;
        
        if count == 0
            pcand = params;
            pdraw = params;
        else
            pcand = pdraw;
        end
        
        pcand(block{bb}) = pdraw(block{bb}) + mvnrnd(zeros(sum(block{bb}),1),prop_cov(block{bb},block{bb}))';
        
        if min(pcand>params_LB)==1 && min(pcand<params_UB)==1  && abs(log_prior(pcand))<inf
            
            % attempt to solve model for parameter candidate
            
            [pcand,fail] = solve_steady(pcand);
            [A_cand,B_cand,D,P_cand,p_cand] = first_guess(params,n_hist,n_shk,n_Y,kbar);
            try
                [A_cand,B_cand,C1_cand,C2_cand,K_cand,gx_pc_cand,gx_pd_cand,gx_rf_cand,fail] = ...
                    solve_AB_mex(pcand,A_cand,B_cand,D,P_cand,p_cand,kbar,n_shk,n_Y,n_hist,jlead1,jlead0,jlag,tol,fail,A_old,B_old,C1,C2,K,K_old,L,gx_pc,gx_pd,gx_rf,gx_pc_old,gx_pd_old,gx_rf_old,A1,A2,A2_big,A3,A3_big,B1,B2,B2_big,e_1,H);
            catch
                [A_cand,B_cand,C1_cand,C2_cand,K_cand,gx_pc_cand,gx_pd_cand,gx_rf_cand,fail] = ...
                    solve_AB(pcand,A_cand,B_cand,D,P_cand,p_cand,kbar,n_shk,n_Y,n_hist,jlead1,jlead0,jlag,tol,fail,A_old,B_old,C1,C2,K,K_old,L,gx_pc,gx_pd,gx_rf,gx_pc_old,gx_pd_old,gx_rf_old,A1,A2,A2_big,A3,A3_big,B1,B2,B2_big,e_1,H);
            end
            
            if fail == 0
                [Pz_cand] = solve_condvcov_mex(A_cand,B_cand,C1_cand,C2_cand,D,n_shk,n_Y,kbar,n_hist,jlag);
                [~,g0_pd_cand,g0_rf_cand,fail] = solve_intercepts_mex(gx_pc_cand,gx_pd_cand,gx_rf_cand,D,Pz_cand,pcand,n_hist,jlead1,jlead0);
                [Pi_cand] = solve_xvar_mex(A_cand,C1_cand,C2_cand,D,K_cand,n_hist,jlag);
            end
            
			% stochastic acceptance/rejection of parameter candidate
            
            if fail == 0
                pcand_prior = log_prior(pcand);
                ll_pcand    = log_likelihood(Y1,Y2,T,sT,pcand,kbar,tau,n_Y,n_shk,not0,jlag,e_1,H,A_cand,B_cand,D,g0_pd_cand,g0_rf_cand,gx_pd_cand,gx_rf_cand,Pi_cand);
                laccratio   = ll_pcand - ll_pdraw + pcand_prior - pdraw_prior + ...
                    log(bern_pmf(T+tau,pcand(22),sum(sT==1))) - log(bern_pmf(T+tau,pdraw(22),sum(sT==1)));
                if isnan(laccratio)
                    laccratio = -inf;
                end
            else
                ll_pcand   = -inf;
                laccratio  = -inf;
                fail_count = fail_count + 1;
            end
            
        else
            laccratio       = -inf;
        end
        
        if ll_pcand + pcand_prior >= lpost_max
            lpost_max   = ll_pcand + pdraw_prior;
            p_max       = pcand;
            sT_max      = sT;
        end
        
        if log(rand) < laccratio
            ll_pdraw    = ll_pcand;
            pdraw_prior = pcand_prior;
            pdraw       = pcand;
            A_draw      = A_cand;
            B_draw      = B_cand;
            g0_pd_draw  = g0_pd_cand;
            g0_rf_draw  = g0_rf_cand;
            gx_pd_draw  = gx_pd_cand;
            gx_rf_draw  = gx_rf_cand;
            Pi_draw     = Pi_cand;
            
            accept_count = accept_count + 1/2;
        end
        
    end
    
	%%% sample block 3, the history of volatility states %%%
    
    sT_cand = sT;
    
	% sample new volatility history candidate
    for tt = 1:T+tau
        if log(rand) <= log(prob_move)
            if sT_cand(tt) == 1
                sT_cand(tt) = 0;
            else
                sT_cand(tt) = 1;
            end
        end
    end
    
    ll_sT_cand = log_likelihood(Y1,Y2,T,sT_cand,pdraw,kbar,tau,n_Y,n_shk,not0,jlag,e_1,H,A_draw,B_draw,D,g0_pd_draw,g0_rf_draw,gx_pd_draw,gx_rf_draw,Pi_draw);
    
    if ll_sT_cand > -9e200
        logaccept_prob_sT  = ll_sT_cand - ll_pdraw + log(bern_pmf(T+tau,pdraw(22),sum(sT_cand==1))) - log(bern_pmf(T+tau,pdraw(22),sum(sT==1)));
        
        if ll_sT_cand + pdraw_prior >= lpost_max
            lpost_max 	= ll_sT_cand + pdraw_prior;
            p_max  		= pdraw;
            sT_max      = sT_cand;
        end
        
        if log(rand) < logaccept_prob_sT
            accept_count_sT = accept_count_sT + 1;
            sT 				= sT_cand;
        end
    end
    
	%%% print markov chain analytics, periodically save data %%%
	
    if mod(count,100) == 0
        toc;
		
		bar(sT);
        
        disp(['Total attempts: ', num2str(count)]);
        disp(['Number params accepted: ', num2str(accept_count)]);
        disp(['Number sT accepted: ', num2str(accept_count_sT)]);
        disp(['(sT == 1) share: ', num2str(sum(sT==1)/(T+tau))]);
        disp(['Recent fails: ', num2str(fail_count)]);
        disp('Posterior max and current draw: ')
        disp(' ');
        disp([p_max, pdraw]);
        disp(' ');
        
        params_chain    = [params_chain; pdraw'];
        sT_chain        = [sT_chain; sT'];
        
        save('params_chain','params_chain');
        save('sT_chain','sT_chain');       
        save('workspace');
        
        outside_count   = 0;
        fail_count      = 0;
        
        tic;
    end
    
end

quit;