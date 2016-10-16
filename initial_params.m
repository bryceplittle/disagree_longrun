function [params] = initial_params()

%% Quaterly Data %%

% gamma   = 10;                                   % risk aversion
% psi     = 1.90;                                 % ies
% beta    = 0.990;                                % time preference
% pc      = 7;                                    % price-cons ratio
% pd      = 7;                                    % price-div ratio
% kap_1   = exp(pc)/(1+exp(pc));                  % linearization constant
% kap_0   = log(1+exp(pc))-kap_1*pc;              % linearization constant
% kap_1m  = exp(pd)/(1+exp(pd));                  % linearization constant
% kap_0m  = log(1+exp(pd))-kap_1m*pd;             % linearization constant
% mu_c    = 0.0043;                               % avg consumption gr rate
% mu_d    = 0.0048;                               % avg dividend gr rate
% rho_x   = 0.99;                                 % lrr persistence
% rho_d   = 3.00;                                 % leverage
% rho_cd  = 1.00;                                 % cons/div cov
% phi_x   = 0.0028*0.038;                         % lrr stdev
% phi_c   = 0.0028;                               % consumption gr stdev
% phi_d   = 0.0300;                               % div gr stdev
% phi_rm  = 0.01;                                 % rm measurement error
% phi_rf  = 0.01;                                 % rf measurement error
% phi_s   = 0.001;                                % signal precision
% sig_h   = 2.00;                                 % vol multiplier
% pi      = 0.05;                                 % prob of low vol state
% phi_ce  = 0.00;                                 % xvar scale

%% Annual Data %%

gamma   = 10;                                   % risk aversion
psi     = 1.5;                                  % ies
beta    = 0.9989;                               % time preference
pc      = 7;                                    % price-cons ratio
pd      = 7;                                    % price-div ratio
kap_1   = exp(pc)/(1+exp(pc));                  % linearization constant
kap_0   = log(1+exp(pc))-kap_1*pc;              % linearization constant
kap_1m  = exp(pd)/(1+exp(pd));                  % linearization constant
kap_0m  = log(1+exp(pd))-kap_1m*pd;             % linearization constant
mu_c    = 0.0015;                               % avg consumption gr rate
mu_d    = 0.0015;                               % avg dividend gr rate
rho_x   = 0.975;                                % lrr persistence
rho_d   = 2.50;                                 % leverage
rho_cd  = 2.60;                                 % cons/div cov
phi_x   = 0.0072*0.038;                         % lrr stdev
phi_c   = 0.0072;                               % consumption gr stdev
phi_d   = 5.96*phi_c;                           % div gr stdev
phi_rm  = 0.02;                                 % rm measurement error
phi_rf  = 0.02;                                 % rf measurement error
phi_s   = 2*phi_x;                        		% signal precision
sig_h   = 2.00;                                 % vol multiplier
pi      = 0.05;                                 % prob of low vol state
phi_ce  = 0.00;                                 % xvar scale


params = [gamma,psi,beta,pc,pd,kap_1,kap_0,kap_1m,kap_0m,mu_c,mu_d,rho_x,rho_d,rho_cd,phi_x,phi_c,phi_d,phi_rm,phi_rf,phi_s,sig_h,pi,phi_ce]';

% params(1)   = gamma;
% params(2)   = psi;
% params(3)   = beta;
% params(4)   = pc;
% params(5)   = pd;
% params(6)   = kap_1;
% params(7)   = kap_0;
% params(8)   = kap_1m;
% params(9)   = kap_0m;
% params(10)  = mu_c;
% params(11)  = mu_d;
% params(12)  = rho_x;
% params(13)  = rho_d;
% params(14)  = rho_cd;
% params(15)  = phi_x;
% params(16)  = phi_c;
% params(17)  = phi_d;
% params(18)  = phi_rm;
% params(19)  = phi_rf;
% params(20)  = phi_s;
% params(21)  = sig_h;
% params(22)  = pi;
% params(23)  = phi_ce;

end

