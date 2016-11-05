function [fval, Bt, brt, bft ] = Euler_sum_int(p,s_par,s_grid,s_state,s_gov,s_investors)

%% VARIABLES NEEDED

%Parameters
sigma = s_par.sigma;
e.f = s_par.e.f;
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
wt = s_state.wt;
eft = e.f(n);

%Government
crt = s_gov.crt;
cft = s_gov.cft;
bgt1 = s_gov.bg0(:,:);
crt1 = s_gov.cr0(:,:);
cft1 = s_gov.cf0(:,:);

%Investors
rt1 = s_investors.r0(:,:);
wt1 = s_investors.w0(:,:);
qt1 = s_investors.q0(:,:);
zt1 = s_investors.z0(:,:);
brt1 = s_investors.br0(:,:);
bft1 = s_investors.bf0(:,:);

%% ALGORITHM

%YOU DON'T NEED TO CONSIDER 'zt' IN THE DENOMINADOR, SINCE 
%THE BOND MARKET MUST BE OPEN FOR THIS CALCULATION DO BE DONE
%YOU DO NEED 'zt1'

denom_g = tc*(crt+cft) + (p*grid_g' - bgt_1);
num_g = (tc*(crt1+cft1) + zt1.*(qt1.*bgt1...
                - repmat(grid_g,n_states,1)));
ratio_g = Euler_ratio(s_par,s_state,zt1,num_g,denom_g,'g');
euler_g = abs(ratio_g - p);

denom_r = ((1+rt)*brt_1 + wt - p*grid_r');
num_r = ((1+rt1).^(-1/sigma.r)).*...
        ((1+rt1).*zt1.*repmat(grid_r,n_states,1) + wt1 - zt1.*qt1.*brt1);
ratio_r = Euler_ratio(s_par,s_state,zt1,num_r,denom_r,'r');
euler_r = abs(p - ratio_r);


denom_f = ((1+rt)*(eft + bft_1) - p*grid_f');
num_f = ((1+rt1).^(-1/sigma.f)).*...
        ((1+rt1).*(repmat(e.f,1,l_grid_g) + zt1.*repmat(grid_f,n_states,1)) - zt1.*qt1.*bft1);
ratio_f = Euler_ratio(s_par,s_state,zt1,num_f,denom_f,'f');
euler_f = abs(p - ratio_f);


[fval , b_star] = min(euler_g + euler_r + euler_f);
Bt = grid_g(b_star);
brt = grid_r(b_star);
bft = grid_f(b_star);

end
