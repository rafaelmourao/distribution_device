% SOVEREIGN DEFAULT AS A DISTRIBUTION DEVICE - November 2016

%% ECONOMY PARAMETERS

param = [];

param.has_default = true;
param.has_partial_default = true;     

%Consumers
param.beta = .95;      %Intertemporal discount rate
param.sigma.r = 2;    %Utility function parameter: risk aversion
param.sigma.f = 2;
param.Ar = 0;     %Resident's fixed income stream
param.Af = 0;     %Foreigner's fixed income stream

%Government
param.sigma.g = 2;    %Utility function parameter: risk aversion
param.phi = .3;    %Probability of redemption (Arellano)
param.lambda = 1/2;    %Government preference parameter: foreigners relative to residents
param.theta = 1;
param.tc = .3;   %Tax rate over CONSUMPTION
param.Ag = .3;     %Govenment's fixed income stream


%Firm
param.alpha = .3;            %Participation of capital on productio
param.rho = -1;              %Elasticity of Substitution between capital and labor is 1/2 (=1/(1-rho))


% Matrix: discretization of AR(1) process
param.n_s = 5;          % Number of states of nature
param.mu = .5;           % Average of foreigner shock
param.gamma = .5;       % Autoregressive coefficient
param.nu = 1;           % Variance of stochastic shock
param.kappa = 2;        % Number of SD that will deviate from the mean to...
                        % form the grid 
[param.e.f , param.prob] = ar1(param);
                        
                        
% Matrix: old version (still working though)
%param.e.f = [.5;10;15];
%param.prob = [.4 .5 .1;.3 .4 .3;.1 .5 .4];


%Foreigner wealth evolution
%param.e.f = [2;3;5];
%param.e.f = [.8;1.5;2.5];
%param.e.f = [1;2;5;8;10];
%param.prob = [.2 .3 .3 .1 .1;...
%               .1 .4 .3 .1 .1;...
%               .1 .3 .3 .2 .1;...
%               .1 .2 .3 .3 .1;...
%               .1 .2 .2 .3 .2];     %Construction od the Probability matrix


% GRID
%Public Bonds
param.min_b = 0;   %Minimum value for bonds
param.max_b = 3;  %Maximum value for bonds
param.n_bonds = 15;  %Quantity of points on the grid for the investors

%% ITERATIONS

iter = Economy(param); % Construct economy with given parameters

epsilon = 1e-2;                                     %Tolerance level for 
dist = 100;                                         %Distance between previous and current price and bond functions
t = 0;                                              %Number of interations

hah = tic;
while dist > epsilon && t <= 10000
    tic
    t = t+1;
    
    old_iter = iter;
    iter = iter.update(20); % Sets maximum number of parallel workers

    time = toc;
    dist = max(abs(iter.Vo(:) - old_iter.Vo(:)));
    fprintf('Iter: %d, distance: %.6f, time: %.2f seconds\n',t,dist,time) 
    disp('Prices:')
    disp(iter.q(1:10))
    fprintf('Default proportion: %.4f\n',1-mean(iter.delta(:))) 
    
end
toc(hah)

addpath('plots/')
plot_quantities(iter)
plot_default(iter)
plot_prices(iter)

function [S, P] = ar1(param)

    % basic parameters
    n_s = param.n_s;
    mu = param.mu;
    gamma = param.gamma;
    nu = param.nu;
    kappa = param.kappa;

    % grid for shock
    nu_s = nu/sqrt(1-gamma^2);
    mu_s = mu/(1-gamma);
    smin = mu_s - kappa*nu_s;
    smax = mu_s + kappa*nu_s;
    sstep = (smax - smin)/(n_s-1);
    S = exp(linspace(smin,smax,n_s))';

    % transition matrix
    normarg1 = (smin - mu - gamma*log(S))/nu + 0.5*sstep/nu;
    P(:,1) = 0.5 + 0.5*erf(normarg1/sqrt(2));
    for j = 2:(n_s-1)
        normarg1 = (log(S(j)) - mu - gamma*log(S))/nu + 0.5*sstep/nu;
        normarg2 = (log(S(j)) - mu - gamma*log(S))/nu - 0.5*sstep/nu;
        P(:,j) = 0.5*erf(normarg1/sqrt(2)) - 0.5*erf(normarg2/sqrt(2));
    end
    P(:,n_s) = 1 - sum(P(:,1:(n_s-1))')';
end