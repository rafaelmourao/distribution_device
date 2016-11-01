function [varargout] = Euler_sum_int(p,s_par,s_grid,s_state,s_gov,s_investors)

%% VARIABLES NEEDED

%Parameters
sigma = s_par.sigma;
lambda = s_par.lambda;
Qr = s_par.Qr;
Qf = s_par.Qf;
tc = s_par.tc;
n_states = s_par.n_states;

%Grid
grid_r = s_grid.grid_r_aux;
grid_f = s_grid.grid_f_aux;
grid_g = s_grid.grid_b_g;
l_grid_g = length(grid_g);

%State
n = s_state.n;
brt_1 = s_state.brt_1;
bft_1 = s_state.bft_1;
bgt_1 = s_state.bgt_1;
rt = s_state.rt;
lambdat = lambda(n);
Qrt = Qr(n);
Qft = Qf(n);
Qg = lambda/2;
Qgt = lambdat/2;

%Government
crt = s_gov.crt;
cwt = s_gov.cwt;
bgt1 = s_gov.bgt1;
crt1 = s_gov.crt1;
cwt1 = s_gov.cwt1;

%Investors
rt1 = s_investors.rt1;
qt1 = s_investors.qt1;
zt1 = s_investors.zt1;
brt1 = s_investors.brt1;
bft1 = s_investors.bft1;

%% ALGORITHM

%JUST FOR YOU TO KNOW -> YOU'RE GOING TO HAVE PROBLEMS WITH
%THE 'zt1' IN ALL THE NUMERATORS: IT CAN BECOME ZERO AND
%IT'S RAISED TO THE POWER OF A NEGATIVE NUMBER. YOU MUST
%FIX THAT LATER

denom_g = Qgt + tc*(crt+cwt) - bgt_1/lambdat + p*grid_g';
num_g = (zt1.^(-1/sigma)).*(Qg*ones(1,l_grid_g) + (tc*lambda)*ones(1,l_grid_g).*(crt1+cwt1) - ones(n_states,1)*grid_g +...
                (lambda*ones(1,l_grid_g)).*qt1.*zt1.*bgt1);
ratio_g = Euler_ratio(s_par,s_state,num_g,denom_g);
euler_g = abs(ratio_g - p);

denom_r = ((1+rt).^((sigma-1)/sigma)).*(Qrt + brt_1/lambdat - p*grid_r');
num_r = ((zt1.^(-1/sigma)).*((1+rt1).^((sigma-1)/sigma))).*...
        ((lambda.*Qr)*ones(1,l_grid_g) + ones(n_states,1)*grid_r - (lambda*ones(1,l_grid_g)).*qt1.*zt1.*brt1);
ratio_r = Euler_ratio(s_par,s_state,num_r,denom_r);
euler_r = abs(p - ratio_r);


denom_f = ((1+rt).^((sigma-1)/sigma)).*(Qft + bft_1/lambdat - p*grid_f');
num_f = ((zt1.^(-1/sigma)).*((1+rt1).^((sigma-1)/sigma))).*...
        ((lambda.*Qf)*ones(1,l_grid_g) + ones(n_states,1)*grid_f - (lambda*ones(1,l_grid_g)).*qt1.*zt1.*bft1);
ratio_f = Euler_ratio(s_par,s_state,num_f,denom_f);
euler_f = abs(p - ratio_f);


[varargout{1} , b_star] = min(euler_g + euler_r + euler_f);

Bt = grid_g(b_star);

brt = grid_r(b_star);

bft = grid_f(b_star);

if nargout == 3
    
    varargout{1} = Bt;
    varargout{2} = brt;
    varargout{3} = bft;
    
    return
    
end
