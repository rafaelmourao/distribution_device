% SOVEREIGN DEFAULT AS A DISTRIBUTION DEVICE - November 2016

%% ECONOMY PARAMETERS

param = [];

param.has_default = true;
param.has_partial_default = true;     

%Consumers
param.beta = .8;      %Intertemporal discount rate
param.sigma.r = 2;    %Utility function parameter: risk aversion
param.sigma.f = 2;
param.Ar = 0;     %Resident's fixed income stream
param.Af = 0;     %Foreigner's fixed income stream

%Government
param.sigma.g = 2;    %Utility function parameter: risk aversion
param.phi = .5;    %Probability of redemption (Arellano)
param.lambda = 1;    %Government preference parameter: foreigners relative to residents
param.tc = .2;   %Tax rate over CONSUMPTION
param.Ag = .65;     %Govenment's fixed income stream


%Firm
param.alpha = .2;            %Participation of capital on productio
param.rho = .5;              %Elasticity of Substitution between capital and labor is 1/2 (=1/(1-rho))

%Foreigner wealth evolution
param.e.f = [1;2;5];
% param.e.f = [0.1;.3;0.5;2;5];

%Transition Matrix
param.prob = [.4 .5 .1;.3 .4 .3;.1 .5 .4];
% param.prob = [.2 .3 .3 .1 .1;...
%               .1 .4 .3 .1 .1;...
%               .1 .3 .3 .2 .1;...
%               .1 .2 .2 .3 .2;...
%               .1 .2 .2 .2 .3];     %Construction od the Probability matrix

% GRID

%Public Bonds
param.min_b = 0;   %Minimum value for bonds
param.max_b = .4;  %Maximum value for bonds
param.n_bonds = 15;  %Quantity of points on the grid for the investors

%% ITERATIONS

iter = Economy(param); % Construct economy with given parameters

epsilon = 1e-5;                                     %Tolerance level for 
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
