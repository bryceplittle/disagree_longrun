function [Params_UB, Params_LB] = param_bounds()

Params_UB       = zeros(23,1);
Params_LB       = zeros(23,1);

Params_LB(1)   = 0;
Params_LB(2)   = 0;
Params_LB(3)   = 0;
Params_LB(4)   = -inf;
Params_LB(5)   = -inf;
Params_LB(6)   = -inf;
Params_LB(7)   = -inf;
Params_LB(8)   = -inf;
Params_LB(9)   = -inf;
Params_LB(10)  = -inf;
Params_LB(11)  = -inf;
Params_LB(12)  = -inf;
Params_LB(13)  = -inf;
Params_LB(14)  = -inf;
Params_LB(15)  = 0;
Params_LB(16)  = 0;
Params_LB(17)  = 0;
Params_LB(18)  = 0;
Params_LB(19)  = 0;
Params_LB(20)  = 0;
Params_LB(21)  = 0;
Params_LB(22)  = 0;
Params_LB(23)  = -inf;

Params_UB(1)   = inf;
Params_UB(2)   = inf;
Params_UB(3)   = 1;
Params_UB(4)   = inf;
Params_UB(5)   = inf;
Params_UB(6)   = inf;
Params_UB(7)   = inf;
Params_UB(8)   = inf;
Params_UB(9)   = inf;
Params_UB(10)  = inf;
Params_UB(11)  = inf;
Params_UB(12)  = inf;
Params_UB(13)  = inf;
Params_UB(14)  = inf;
Params_UB(15)  = inf;
Params_UB(16)  = inf;
Params_UB(17)  = inf;
Params_UB(18)  = inf;
Params_UB(19)  = inf;
Params_UB(20)  = inf;
Params_UB(21)  = inf;
Params_UB(22)  = 1;
Params_UB(23)  = inf;

end

