%% SOVEREIGN DEFAULT AS A DISTRIBUTION DEVICE - November 2016

% ASSUMPTIONS:  %BOND MARKET: TELMER ARTICLE TO SOLVE EULER EQUATIONS
%DEFAULT DECISION IS FOUND BY A LINEAR PROGRAMMING METHOD
%CONCAVE UTILITY FUNCTIONS
%NO PRODUCTIVITY SHOCKS ONTO THE FIRM

%CONSTANT TAX RATE

%SYMMETRIC TRANSITION MATRIX AMONG STATES OF NATURE
%REGARDING SHOCKS

%THERE IS WEALTH GROWTH AND IT IS STOCHASTIC: THERE WILL BE
%A GOOD AND A BAD STATE OF NATURE

%IN THE GOOD STATE, WEALTH IS SPLITTED EVENLY AMONG
%CONSUMERS

%IN THE BAD STATE THE DISTRIBUTION OF WEALTH IS
%STOCHASTIC, BUT THE LIKELIHOOD OF EITHER A BAD OR GOOD
%STATES ARE THE SAME

%WITH THAT, WE CAN REWRITE ALL VARIABLES, AND MAKE THEM
%RELATIVE TO AGGREGATE WEALTH

%BOND'S MARKET WILL BE CONSIDERED CLEARED FOR ALL THE
%STATES

%% PARAMETERS

epsilon = 1e-3;                                     %Tolerance level

%Consumers
param.beta = .95;                                         %Intertemporal discount rate
param.sigma.r = 2;                                        %Utility function parameter: risk aversion
param.sigma.f = 2;

%Government
param.sigma.g = 2;                                        %Utility function parameter: risk aversion
param.phi = .282;                                         %Probability of redemption (Arellano)
param.lambda = 1;                                         %Government preference parameter: foreigners relative to residents
param.tc = .3;                                            %Tax rate over CONSUMPTION

%Firm
param.alpha = .3;                                         %Participation of capital on productio
param.rho = -1;                                           %Elasticity of Substitution between capital and labor is 1/2 (=1/(1-rho))

%Foreigner wealth evolution
param.e.f = [.1;.2;.3];

%Transition Matrix
param.prob = [.3 .6 .1;.2 .6 .2;.2 .5 .3];                %Construction od the Probability matrix
param.n_states = size(param.prob,1);                            %Numbers of States of Nature

%% GRID

%Public Bonds
min_b = 0;                                          %Minimum value for bonds
max_b = 2;                                         %Maximum value for bonds
param.n_bonds = 21;                                       %Quantity of points on the grid for the investors

param.grid.b_r = linspace(min_b,max_b,param.n_bonds);            %Grid for resident bonds:
param.grid.b_f = param.grid.b_r;                                 %Grid for foreigner bonds;
param.grid.r_aux = repmat(param.grid.b_r,1,param.n_bonds);             %Combinations of bonds
param.grid.f_aux = kron(param.grid.b_f,ones(1,param.n_bonds));
param.grid.b_g = param.grid.r_aux + param.grid.f_aux;
param.n_bonds_g = length(param.grid.b_g);                       %Number of combination of bonds


iter = Economy(param);

dist = 100;                                         %Distance between previous and current price and bond functions
t = 1;                                              %Number of interations
while dist > epsilon && t <= 1
    tic
    t = t+1;
    
    old_iter = iter;
    iter = iter.update(1);

    time = toc;
    dist = max(abs(iter.Vo(:) - old_iter.Vo(:)));
    fprintf('Iter: %d, distance: %.6f, time: %.2f seconds\n',t,dist,time) 
    disp('Prices:')
    disp(iter.q(1:10))
    
end

toc