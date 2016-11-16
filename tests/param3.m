param.p_max = 10;
%Consumers
param.beta = .8;                                         %Intertemporal discount rate
param.sigma.r = 2;                                        %Utility function parameter: risk aversion
param.sigma.f = 2;

%Government
param.sigma.g = 2;                                        %Utility function parameter: risk aversion
param.phi = .1;                                         %Probability of redemption (Arellano)
param.lambda = 1;                                         %Government preference parameter: foreigners relative to residents
param.tc = .2;                                            %Tax rate over CONSUMPTION
param.A = 0;                                              % Fixed income stream

%Firm
param.alpha = .3;                                         %Participation of capital on productio
param.rho = .5;                                           %Elasticity of Substitution between capital and labor is 1/2 (=1/(1-rho))

%Foreigner wealth evolution
param.e.f = [1;2;5];
% param.e.f = [0.1;.3;0.5;2;5];

%Transition Matrix
param.prob = [.4 .5 .1;.3 .4 .3;.1 .5 .4];
% param.prob = [.2 .3 .3 .1 .1;...
%               .1 .4 .3 .1 .1;...
%               .1 .3 .3 .2 .1;...
%               .1 .2 .2 .3 .2;...
%               .1 .2 .2 .2 .3];                %Construction od the Probability matrix
param.n_states = size(param.prob,1);                            %Numbers of States of Nature
% GRID

%Public Bonds
param.min_b = 0;                                          %Minimum value for bonds
param.max_b = .2;                                         %Maximum value for bonds
param.n_bonds = 10;                                       %Quantity of points on the grid for the investors
