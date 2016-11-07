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
beta = .95;                                         %Intertemporal discount rate
sigma.r = 2;                                        %Utility function parameter: risk aversion
sigma.f = 2;

%Government
sigma.g = 2;                                        %Utility function parameter: risk aversion
phi = .282;                                         %Probability of redemption (Arellano)
lambda = 1;                                         %Government preference parameter: foreigners relative to residents
tc = .3;                                            %Tax rate over CONSUMPTION

%Firm
alpha = .3;                                         %Participation of capital on productio
rho = -1;                                           %Elasticity of Substitution between capital and labor is 1/2 (=1/(1-rho))

%Foreigner wealth evolution
e.f = [.1;.2;.3];

%Transition Matrix
prob = [.3 .6 .1;.2 .6 .2;.2 .5 .3];                %Construction od the Probability matrix
n_states = size(prob,1);                            %Numbers of States of Nature

s_par = struct('epsilon',epsilon,'beta',beta,'sigma',sigma,...
    'phi',phi,'lambda',lambda,'tc',tc,'alpha',alpha,'rho',rho,...
    'e',e,'prob',prob,'n_states',n_states);

%% GRID

%Public Bonds
min_b = 0;                                          %Minimum value for bonds
max_b = 1;                                         %Maximum value for bonds
n_bonds = 21;                                       %Quantity of points on the grid for the investors

grid_b_r = linspace(min_b,max_b,n_bonds);            %Grid for resident bonds:
grid_b_f = grid_b_r;                                 %Grid for foreigner bonds;
grid_r_aux = repmat(grid_b_r,1,n_bonds);             %Combinations of bonds
grid_f_aux = kron(grid_b_f,ones(1,n_bonds));
grid_b_g = grid_r_aux + grid_f_aux;
n_bonds_g = length(grid_b_g);                       %Number of combination of bonds

s_grid = struct('grid_b_r',grid_b_r,'grid_b_f',grid_b_f,...
    'grid_b_g',grid_b_g,'grid_r_aux',grid_r_aux,'grid_f_aux',grid_f_aux);

%% VALUE FUNCTIONS FOR THE GOVERNMENT

Vd = zeros(n_states,1);                         %Value function in case of DEFAULT
Vnd = zeros(n_states,n_bonds,n_bonds);          %Value function in case of NO DEFAULT
Vo = zeros(n_states,n_bonds,n_bonds);           %Value function for the OPTION of DEFAULT

%% POLICY FUNCTIONS

%Consumers

cr = zeros(n_states,n_bonds,n_bonds);          %CONSUMPTION policy funtion for RESIDENTS
cf = zeros(n_states,n_bonds,n_bonds);          %CONSUMPTION policy funtion for FOREIGNERS

kr = zeros(n_states,n_bonds,n_bonds);          %CAPITAL policy funtion for RESIDENTS
kf = zeros(n_states,n_bonds,n_bonds);          %CAPITAL policy funtion for FOREIGNERS

br = zeros(n_states,n_bonds,n_bonds);          %BONDS policy funtion for RESIDENTS
bf = zeros(n_states,n_bonds,n_bonds);          %BONDS policy funtion for FOREIGNERS

%Government
bg = 1e-10*ones(n_states,n_bonds,n_bonds);     %BONDS policy funtion for the GOVERNMENT
g = zeros(n_states,n_bonds,n_bonds);           %PUBLIC EXPENDITURE policy funtion for the GOVERNMENT
z = ones(n_states,n_bonds,n_bonds);            %DEFAULT policy funtion for the GOVERNMENT

%% PRICE FUNCTIONS

r = zeros(n_states,n_bonds,n_bonds);           %Interest Rate
w = zeros(n_states,n_bonds,n_bonds);           %Wage
q = zeros(n_states,n_bonds,n_bonds);           %Price of Public Bond

%% ITERATION

%Creation of Functions to Iteration
Vo1 = Vo;
Vd1 = Vd;
Vnd1 = Vnd;

br1 = br;
bf1 = bf;

kr1 = repmat([0;0;0],[1 n_bonds n_bonds]);    %Residents don't have any resources to lend in the very beginning
kf1 = repmat(e.f,[1 n_bonds n_bonds]);        %Foreigners only have their endowments to lend in the very beginning

r1 = r + alpha*((kr1+kf1).^(rho-1)).*((alpha*((kr1+kf1).^rho) + (1-alpha)).^(1/rho-1));
w1 = w + (1-alpha)*((alpha*((kr1+kf1).^rho) + (1-alpha)).^(1/rho-1));
q1 = q + 1;

cr1 = cr + (1/(1+tc))*((1+r1).*kr1 + w1);
cf1 = cf + (1/(1+tc))*(1+r1).*kf1;

bg1 = bg;
g1 = g + tc.*(cr1 + cf1);
z1 = z;

% Calculation of default outcome
r_d = zeros(n_states,1); % Variables in case of default
w_d = zeros(n_states,1);
cr_d = zeros(n_states,1);
cf_d = zeros(n_states,1);
g_d = zeros(n_states,1);
Wd = zeros(n_states,1);
for n = 1:n_states
    r_d(n) = ...
        alpha*(e.f(n).^(rho-1)).*((alpha*(e.f(n).^rho) + (1-alpha)).^(1/rho-1));
    w_d(n) = ...
        (1-alpha)*((alpha*(e.f(n))^(rho)) + (1-alpha))^(1/rho-1);
    cr_d(n) = (1/(1+tc))*w_d(n);
    cf_d(n) = (1/(1+tc))*(1+r_d(n))*e.f(n);
    g_d(n) = tc*(cr_d(n) + cf_d(n));
    Wd(n) =  Utility_Function(cr_d(n),sigma.r) +  lambda*Utility_Function(cf_d(n),sigma.f) ...
        +  Utility_Function(g_d(n),sigma.g);
end

dist = 100;                                         %Distance between previous and current price and bond functions
t = 1;                                              %Number of interations
while dist > epsilon && t <= 200
    
    tic
    t = t+1;
    
    %Updating process of Iteration
    %Observation: For the first iteration, the number of bonds owned by
    %each agent is ZERO, because it is considered the last period and would
    %be no future period in which the bonds would be paid back.
    Vo0 = Vo1;
    Vd0 = Vd1;
    Vnd0 = Vnd1;
    
    kr0 = kr1;
    kf0 = kf1;
    
    cr0 = cr1;
    cf0 = cf1;
    
    br0 = br1;
    bf0 = bf1;
    
    bg0 = bg1;
    g0 = g1;
    z0 = z1;
    
    r0 = r1;
    w0 = w1;
    q0 = q1;
    
    
    s_investors = struct('r0',r0,'w0',w0,'q0',q0,'z0',z0,'br0',br0,'bf0',bf0);

    Vd1 = Wd + beta * prob * ( phi * Vo0(:,1,1) + (1-phi) * Vd0 );
    
    for id_br = 1:n_bonds                               %RESIDENTS bonds from previous period
        brt_1 = grid_b_r(id_br);
        
        for id_bf = 1:n_bonds                      %FOREIGNERS bonds from previous period
            bft_1 = grid_b_f(id_bf);
            
            for n = 1:n_states                          %State of Nature
                
                bgt_1 = brt_1 + bft_1;                  %Amount of bonds issued by the Government in the previous period
                probt = prob(n,:);                      %Probability measure of next period's states of nature (Row Vector)
                rt = r0(n,id_br,id_bf);                 %Current Interest Rate
                wt = w0(n,id_br,id_bf);                 %Current Interest Rate
                zt = z0(n,id_br,id_bf);                 %Current Default Decision
                
                s_state = struct('rt',rt,'wt',wt,'zt',zt,'n',n,...
                    'brt_1',brt_1,'bft_1',bft_1,'bgt_1',bgt_1);
                
                %Variables that agents will use as known values when making their decisions
                
                %Government
                crt = cr0(n,id_br,id_bf);               %Current Consumption for RESIDENT investor
                cft = cf0(n,id_br,id_bf);               %Current Consumption for FOREIGN investor
                
                s_gov = struct('crt',crt,'cft',cft,'bg0',bg0,...
                    'cr0', cr0, 'cf0', cf0);
                
                %**********%
                [p, br_s, bf_s, bg_s] = ...
                    Solution(s_par,s_grid,s_state,s_gov,s_investors);
                id_br_s = find(grid_b_r == br_s);
                id_bf_s = find(grid_b_f == bf_s);
                %**********%
                
                q1(n,id_br,id_bf) = p;                              %Equilibrium price for the Bonds' market
                
                br1(n,id_br,id_bf) = br_s;                          %RESIDENTS demand for Bonds
                bf1(n,id_br,id_bf) = bf_s;                          %FOREIGNERS demand for Bonds
                bg1(n,id_br,id_bf) = bg_s;                          %GOVERNEMNT supply of Bonds
                
                kr1(n,id_br,id_bf) = brt_1;                         %RESIDENTS capital supply
                kf1(n,id_br,id_bf) = e.f(n) + bft_1;                %FOREIGNERS capital supply
                
                r1(n,id_br,id_bf) = ...
                    alpha*((kr1(n,id_br,id_bf)+kf1(n,id_br,id_bf))^(rho-1))*...
                    ((alpha*((kr1(n,id_br,id_bf)+kf1(n,id_br,id_bf))^(rho)) + (1-alpha))^(1/rho-1));            %Interest rate
                w1(n,id_br,id_bf) = ...
                    (1-alpha)*(alpha*((kr1(n,id_br,id_bf)+kf1(n,id_br,id_bf))^(rho)) + (1-alpha))^(1/rho-1);    %Wage
                
                cr1(n,id_br,id_bf) = (1/(1+tc))*((1+r1(n,id_br,id_bf))*kr1(n,id_br,id_bf)...
                    + w1(n,id_br,id_bf) - p*br_s);                                                              %RESIDENT consumption
                cf1(n,id_br,id_bf) = (1/(1+tc))*((1+r1(n,id_br,id_bf))*kf1(n,id_br,id_bf)...
                    - p*bf_s);                                                                                  %FOREIGN investors' consumption
                
                g1(n,id_br,id_bf) = tc*(cr1(n,id_br,id_bf) + cf1(n,id_br,id_bf)) + p*bg_s - bgt_1;
                
                Wnd = Utility_Function(cr1(n,id_br,id_bf),sigma.r) + ...
                    lambda*Utility_Function(cf1(n,id_br,id_bf),sigma.f) + ...
                    Utility_Function(g1(n,id_br,id_bf),sigma.g);
                
                Vnd1(n,id_br,id_bf) = Wnd + beta*(probt * Vo0(:,id_br_s,id_bf_s));
                
                if ( Vnd1(n,id_br,id_bf) > Vd1(n) )
                    z1(n,id_br,id_bf) = 1;
                    Vo1(n,id_br,id_bf) = Vnd1(n,id_br,id_bf);
                else
                    z1(n,id_br,id_bf) = 0;
                    Vo1(n,id_br,id_bf) = Vd1(n);
                end
                
            end
            
        end
        
    end
    
    time = toc;
    fprintf('Iter: %d, distance: %.6f, time: %.2f seconds\n',t,dist,time) 
    disp('Prices in the 1 state:')
    disp(q1(1:10))
    
    
end

toc