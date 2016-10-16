function [LP] = log_prior(params)

LP      = 0;

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

%% Preferences %%

% gamma
a = 1;
b = 20;

prior_pd = makedist('Uniform','lower',a,'upper',b);

LP = LP + log(pdf(prior_pd, gamma));

% psi
a = 1e-1;
b = 2;

prior_pd = makedist('Uniform','lower',a,'upper',b);

LP = LP + log(pdf(prior_pd, psi));

% beta
a = 0.900;
b = 0.999999;

prior_pd = makedist('Uniform','lower',a,'upper',b);

LP = LP + log(pdf(prior_pd, beta));

%% Consumption Growth %%

% mu_c
a = -0.06;
b =  0.06;

prior_pd = makedist('Uniform','lower',a,'upper',b);

LP = LP + log(pdf(prior_pd, mu_c));

% phi_c
a = 0;
b = 1;

prior_pd = makedist('Uniform','lower',a,'upper',b);

LP = LP + log(pdf(prior_pd, phi_c));

% rho_x
a = -0.999999;
b =  0.999999;

prior_pd = makedist('Uniform','lower',a,'upper',b);

LP = LP + log(pdf(prior_pd, rho_x));

% phi_x
a = 1e-5;
b = 1e-1;

prior_pd = makedist('Uniform','lower',a,'upper',b);

LP = LP + log(pdf(prior_pd, phi_x));

%% Dividend Growth %%

% mu_d
a = -0.06;
b =  0.06;

prior_pd = makedist('Uniform','lower',a,'upper',b);

LP = LP + log(pdf(prior_pd, mu_d));

% rho_d
a = -10;
b =  10;

prior_pd = makedist('Uniform','lower',a,'upper',b);

LP = LP + log(pdf(prior_pd, rho_d));

% rho_cd
a = -10;
b =  10;

prior_pd = makedist('Uniform','lower',a,'upper',b);

LP = LP + log(pdf(prior_pd, rho_cd));

% phi_d
a = 0;
b = 1;

prior_pd = makedist('Uniform','lower',a,'upper',b);

LP = LP + log(pdf(prior_pd, phi_d));

%% Noise %%

% phi_s
a = 0;
b = 1;

prior_pd = makedist('Uniform','lower',a,'upper',b);

LP = LP + log(pdf(prior_pd, phi_s));

% phi_rm
a = 1e-2/2;
b = 2e-2;

prior_pd = makedist('Uniform','lower',a,'upper',b);

LP = LP + log(pdf(prior_pd, phi_rm));

% phi_rf
a = 1e-2/2;
b = 2e-2;

prior_pd = makedist('Uniform','lower',a,'upper',b);

LP = LP + log(pdf(prior_pd, phi_rf));

%% Volatility %%

% sig_h
a = 0;
b = 1;

prior_pd = makedist('Uniform','lower',a,'upper',b);

LP = LP + log(pdf(prior_pd, sig_h));

% pi
a = 0;
b = 1;

prior_pd = makedist('Uniform','lower',a,'upper',b);

LP = LP + log(pdf(prior_pd, pi));

% phi_ce
% a = 0;
% b = 1;

% prior_pd = makedist('Uniform','lower',a,'upper',b);

% LP = LP + log(pdf(prior_pd, phi_ce));

end

