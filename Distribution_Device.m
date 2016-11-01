%% SOVEREIGN DEFAULT AS A DISTRIBUTION DEVICE - First version - 2014

% ASSUMPTIONS:  %BOND MARKET: TELMER ARTICLE TO SOLVE EULER EQUATIONS

%DEFAULT DECISION IS FOUND BY A LINEAR PROGRAMMING METHOD

%CONCAVE UTILITY FUNCTIONS

%NO PRODUCTIVITY SHOCKS ONTO THE FIRM

%CONSTANT TAX RATE
%THE WEALTH SHOCKS ARE MODELLED THE SAME WAY AS IN
%MANKIW(1986)

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
gamma = .15;                                        %Degree of wealth heterogeneity
Qr = [1;1+gamma;1-gamma];                        %(Fraction of)Wealth vector for the RESIDENT consumer
Qf = [1;1-gamma;1+gamma];                        %(Fraction of)Wealth vector for the FOREIGN consumer                                                                                        
sigma.r = 2;                                     %Utility function parameter: risk aversion
sigma.f = 2;

%Firm
alpha = .3;                                         %Participation of capital on production (Sosa)

%Government
sigma.g = 2;                                        %Utility function parameter: risk aversion
phi = .282;                                         %Probability of redemption (Arellano)
theta = .2;                                         %Government preference parameter: foreigners relative to residents
tc = .3;                                            %Tax rate over CONSUMPTION

%Work Sector
eta = .4;                                          %Proportion of labor sector on total initial wealth

%Wealth evolution
lambda = [1.2;.8;.8];

%Transition Matrix
psi = .8;                                           %Parameter of the transition matrix
psi2 = .5;                                          %Parameter for the probability of a bigger wealth in case of a bad shock
w_matrix = [psi 1-psi;1-psi psi];                   %Transition matrix for wealth growth rate
prob = [1 0;0 1;0 1]*w_matrix*[1 0 0;0 1 1];        %Construction od the Probability matrix
prob(:,2) = prob(:,2)*psi2;
prob(:,3) = prob(:,3)*(1-psi2);
n_states = size(prob,1);                            %Numbers of States of Nature

s_par = struct('epsilon',epsilon,'beta',beta,'Qr',Qr,'Qf',Qf,'sigma',...
    sigma,'lambda',lambda,'n_states',n_states,'prob',prob,'tc',tc);

%% GRID

%Public Bonds
min_b = 0;                                          %Minimum value for bonds
max_b = .8;                                         %Maximum value for bonds
grid_b_d_l = 30;                                    %Quantity of points on the grid for the investors
width_b = (max_b-min_b)/(grid_b_d_l-1);             %Grid's width
grid_b_r = min_b:width_b:max_b;                     %Grid for public bonds: RESIDENTS
grid_b_f = grid_b_r + (1:grid_b_d_l)*1e-10;         %Grid for public bonds:  #### WHY???
[X, Y] = meshgrid(grid_b_r,grid_b_f);
b_g_finder = X+Y;                                   %This is an auxiliary object used to identify a value from the GOVERNMENT's grid into it's policy function
grid_b_g = unique(b_g_finder)';                     %Grid for public bonds: GOVERNMENT
grid_b_g_l = length(grid_b_g);                      %Number of points on the government grid
r_func_aux = zeros(1,grid_b_g_l);                   %This vector is used to evaluate the functions for the Resident's Euler Equation
f_func_aux = zeros(1,grid_b_g_l);                   %This vector is used to evaluate the functions for the Foregin's Euler Equation
g_func_aux = zeros(1,grid_b_g_l);                   %This vector is used to evaluate the functions for the Government's Euler Equation
for i = 1:grid_b_g_l
    g_func_aux(i) = find(b_g_finder == grid_b_g(i));
    [f_func_aux(i) r_func_aux(i)] = find(b_g_finder == grid_b_g(i));
end
grid_r_aux = grid_b_r(r_func_aux);
grid_f_aux = grid_b_f(f_func_aux);

s_grid = struct('grid_b_r',grid_b_r,'grid_b_f',grid_b_f,'grid_b_g',grid_b_g,'g_func_aux',g_func_aux,'grid_r_aux',grid_r_aux,'grid_f_aux',grid_f_aux,...
    'b_g_finder',b_g_finder);

%% VALUE FUNCTIONS FOR THE GOVERNMENT

Vd = zeros(n_states,1);                             %Value function in case of DEFAULT
Vnd = zeros(n_states,grid_b_d_l,grid_b_d_l);        %Value function in case of NO DEFAULT
Vo = zeros(n_states,grid_b_d_l,grid_b_d_l);         %Value function for the OPTION of DEFAULT

%% POLICY FUNCTIONS

%Consumers
kf = zeros(n_states,grid_b_d_l,grid_b_d_l);          %CAPITAL policy funtion for FOREIGNERS

cr = zeros(n_states,grid_b_d_l,grid_b_d_l);          %CONSUMPTION policy funtion for RESIDENTS
cf = zeros(n_states,grid_b_d_l,grid_b_d_l);          %CONSUMPTION policy funtion for FOREIGNERS
cw = zeros(n_states,grid_b_d_l,grid_b_d_l);          %CONSUMPTION policy funtion for WORKERS
kr = zeros(n_states,grid_b_d_l,grid_b_d_l);          %CAPITAL policy funtion for RESIDENTS

br = zeros(n_states,grid_b_d_l,grid_b_d_l);          %BONDS policy funtion for RESIDENTS
bf = zeros(n_states,grid_b_d_l,grid_b_d_l);          %BONDS policy funtion for FOREIGNERS

%Government
bg = 1e-10*ones(n_states,grid_b_d_l,grid_b_d_l);     %BONDS policy funtion for the GOVERNMENT
g = zeros(n_states,grid_b_d_l,grid_b_d_l);           %PUBLIC EXPENDITURE policy funtion for the GOVERNMENT
z = ones(n_states,grid_b_d_l,grid_b_d_l);            %DEFAULT policy funtion for the GOVERNMENT

%% PRICE FUNCTIONS

r = zeros(n_states,grid_b_d_l,grid_b_d_l);           %Interest Rate
w = zeros(n_states,grid_b_d_l,grid_b_d_l);           %Wage
q = zeros(n_states,grid_b_d_l,grid_b_d_l);           %Price of Public Bond

%% ITERATION

%Creation of Functions to Iteration
br1 = br;
bf1 = bf;

kr1 = repmat(Qr,[1 grid_b_d_l grid_b_d_l]);
kf1 = repmat(Qf,[1 grid_b_d_l grid_b_d_l]);

r1 = r + alpha*((kr1+kf1)/eta).^(alpha-1);
w1 = w + (1-alpha)*((kr1+kf1)/eta).^alpha;
q1 = q;

cr1 = cr + (1-tc)*(1+r1).*kr1;
cf1 = cf + (1+r1).*kf1;
cw1 = cw + (1-tc)*eta.*w1;

bg1 = bg;
g1 = g + tc.*(cr1 + cw1);
z1 = z;

dist = 100;                                         %Distance between previous and current price and bond functions
t = 1;                                              %Number of interations

while dist > epsilon && t <= 200
    
    tic
    
    %Updating process of Iteration
    %Observation: For the first iteration, the number of bonds owned by
    %each agent is ZERO, because it is considered the last period and would
    %be no future period in which the bonds would be paid.
    kr0 = kr1;
    kf0 = kf1;
    
    cr0 = cr1;
    cf0 = cf1;
    cw0 = cw1;
    
    br0 = br1;
    bf0 = bf1;
    
    bg0 = bg1;
    g0 = g1;
    z0 = z1;
    
    r0 = r1;
    w0 = w1;
    q0 = q1;
    
    future = @(x) [x(1,g_func_aux);x(2,g_func_aux);x(3,g_func_aux)];
    rt1 = future(r0);         %Future interest rate
    qt1 = future(q0);         %Future bond price
    zt1 = future(z0);         %Future Default
    brt1 = future(br0);     %RESIDENT's future bond supply
    bft1 = future(bf0);     %FOREIGN's future bond supply
    bgt1 = future(bg0);     %GOVERNMENT's future bond supply
    crt1 = future(cr0);     %Future consumption for RESIDENT investor
    cwt1 = future(cw0);     %Future consumption for RESIDENT worker
    
    s_investors = struct('r0',r0,'q0',q0,'z0',z0,'br0',br0,'bf0',bf0,'rt1',rt1,'qt1',qt1,'zt1',zt1,'brt1',brt1,'bft1',bft1);
    
    parfor id_br = 1:grid_b_d_l                               %RESIDENTS bonds from previous period
        brt_1 = grid_b_r(id_br);
        
        for id_bf = 1:grid_b_d_l                      %FOREIGNERS bonds from previous period
            bft_1 = grid_b_f(id_bf);
            
            for n = 1:n_states                          %State of Nature
                
                bgt_1 = brt_1 + bft_1;                  %Amount of bonds issued by the Government in the previous period
                lambdat = lambda(n);                    %Wealth growth rate: states 2 and 3 have the same lambda
                Qrt = Qr(n);                            %Share of total income in hands of RESIDENTS
                Qft = Qf(n);                            %Share of total income in hands of FOREIGNERS
                probt = prob(n,:);                      %Probability measure of next period's states of nature (Row Vector)
                rt = r0(n,id_br,id_bf);                 %Current Interest Rate
                
                s_state = struct('rt',rt,'n',n,'brt_1',brt_1,'bft_1',bft_1,'bgt_1',bgt_1);
                
                %Variables that agents will use as known values when making their decisions
                
                %Government
                crt = cr0(n,id_br,id_bf);               %Current Consumption for RESIDENT investor
                cwt = cw0(n,id_br,id_bf);               %Current Consumption for RESIDENT worker
                
                s_gov = struct('crt',crt,'cwt',cwt,'bgt1',bgt1,'bg0',bg0,'crt1',crt1,'cr0',cr0,'cwt1',cwt1,'cw0',cw0);
                
                
                %**********%
                
                if dist < 10
                    
                    vamo_ver = 1;
                    
                end
                
                [p br_s bf_s bg_s] = Solution(s_par,s_grid,s_state,s_gov,s_investors);
                
                %**********%
                
                q1(n,id_br,id_bf) = p;     %Equilibrium price for the Bonds' market
                
                br1(n,id_br,id_bf) = br_s;	%RESIDENTS demand for Bonds
                bf1(n,id_br,id_bf) = bf_s;	%FOREIGNERS demand for Bonds
                bg1(n,id_br,id_bf) = bg_s; 	%GOVERNEMNT supply of Bonds
                
                kr1(n,id_br,id_bf) = Qrt + brt_1/lambdat - p*br_s;      %RESIDENTS capital supply
                kf1(n,id_br,id_bf) = Qft + bft_1/lambdat - p*bf_s;       %FOREIGNERS capital supply
                
                r1(n,id_br,id_bf) = alpha*((kr1(n,id_br,id_bf)+kf1(n,id_br,id_bf))/eta)^(alpha-1);  %Interest rate
                w1(n,id_br,id_bf) = (1-alpha)*((kr1(n,id_br,id_bf)+kf1(n,id_br,id_bf))/eta)^alpha;  %Wage
                
                cr1(n,id_br,id_bf) = (1-tc)*(1+r1(n,id_br,id_bf))*kr1(n,id_br,id_bf);                  %RESIDENT investors' consumption
                cw1(n,id_br,id_bf) = (1-tc)*w1(n,id_br,id_bf)*eta;                              %RESIDENT workers' consumption
                cf1(n,id_br,id_bf) = (1+r1(n,id_br,id_bf))*kf1(n,id_br,id_bf);                  %FOREIGN investors' consumption
                
            end
            
        end
        
    end
    
    t = t + 1;
    
    dist = sum(q1(:) - q0(:))^2
    
    toc
    
    
end

toc