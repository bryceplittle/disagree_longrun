function [params, fail] = solve_steady(params)

options = optimset('MaxIter',1e20,'MaxFunEvals',1e20,'TolFun',1e-6,'Display','off');
fail    = 0;

try
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

end

