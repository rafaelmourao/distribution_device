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

%% ECONOMY PARAMETERS

param = [];
param.p_max = 10;
%Consumers
param.beta = .8;                                         %Intertemporal discount rate
param.sigma.r = 2;                                        %Utility function parameter: risk aversion
param.sigma.f = 2;

%Government
param.sigma.g = 2;                                        %Utility function parameter: risk aversion
param.phi = .2;                                         %Probability of redemption (Arellano)
param.lambda = 1;                                         %Government preference parameter: foreigners relative to residents
param.tc = .2;                                            %Tax rate over CONSUMPTION
param.A = 1;                                              % Fixed income stream

%Firm
param.alpha = .3;                                         %Participation of capital on productio
param.rho = .5;                                           %Elasticity of Substitution between capital and labor is 1/2 (=1/(1-rho))

%Foreigner wealth evolution
% param.e.f = [1;2;3];
param.e.f = [0.1;.3;0.5;2;5];

%Transition Matrix
% param.prob = [.4 .5 .1;.3 .4 .3;.1 .5 .4];
param.prob = [.2 .3 .3 .1 .1;...
              .1 .4 .3 .1 .1;...
              .1 .3 .3 .2 .1;...
              .1 .2 .2 .3 .2;...
              .1 .2 .2 .2 .3];                %Construction od the Probability matrix
param.n_states = size(param.prob,1);                            %Numbers of States of Nature

% GRID

%Public Bonds
param.min_b = 0;                                          %Minimum value for bonds
param.max_b = .2;                                         %Maximum value for bonds
param.n_bonds = 15;                                       %Quantity of points on the grid for the investors

%%

iter = Economy(param);

epsilon = 1e-3;                                     %Tolerance level
dist = 100;                                         %Distance between previous and current price and bond functions
t = 1;                                              %Number of interations

while dist > epsilon && t <= 10000
    tic
    t = t+1;
    
    old_iter = iter;
    iter = iter.update(20);

    time = toc;
    dist = max(abs(iter.Vo(:) - old_iter.Vo(:)));
    fprintf('Iter: %d, distance: %.6f, time: %.2f seconds\n',t,dist,time) 
    disp('Prices:')
    disp(iter.q(1:10))
    fprintf('Default proportion: %.4f\n',1-mean(iter.z(:))) 
end

toc
