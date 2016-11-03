function [fval, Bt, bft] = Euler_sum_cor2(p,s_par,s_grid,s_state,s_gov,s_investors)

%This is the case where the RESISENT investor chooses zero bonds:'br = 0'

%% VARIABLES NEEDED

%Parameters
sigma = s_par.sigma;
e.f = s_par.e.f;
tc = s_par.tc;
n_states = s_par.n_states;

%Grid
grid_f = s_grid.grid_b_f;
l_grid_f = length(grid_f);

%State
n = s_state.n;
bft_1 = s_state.bft_1;
bgt_1 = s_state.bgt_1;
rt = s_state.rt;
eft = e.f(n);

%Investors
r0 = s_investors.r0;
q0 = s_investors.q0;
z0 = s_investors.z0;
bf0 = s_investors.bf0;
bft1 = squeeze(bf0(:,1,:));

%Government
crt = s_gov.crt;
rt1 = squeeze(r0(:,1,:));
qt1 = squeeze(q0(:,1,:));
zt1 = squeeze(z0(:,1,:));
bgt1 = squeeze(s_gov.bg0(:,1,:));
crt1 = squeeze(s_gov.cr0(:,1,:));
cft1 = squeeze(s_gov.cf0(:,1,:));

%% ALGORITHM

%JUST FOR YOU TO KNOW -> YOU'RE GOING TO HAVE PROBLEMS WITH
%THE 'zt1' IN ALL THE NUMERATORS: IT CAN BECOME ZERO AND
%IT'S RAISED TO THE POWER OF A NEGATIVE NUMBER. YOU MUST FIX THAT LATER

denom_g = tc*(crt) - bgt_1 + p*grid_f';
num_g = (tc*(crt1+cft1) - ones(n_states,1)*grid_f +...
                qt1.*zt1.*bgt1);
ratio_g = Euler_ratio(s_par,s_state,zt1,num_g,denom_g,'g');
euler_g = abs(ratio_g - p);
                

denom_f = ((1+rt).^(-1/sigma.f)).*((1+rt)*(eft + bft_1) - p*grid_f');
num_f = ((1+rt1).^(-1/sigma.f)).*...
        ((1+rt1).*(repmat(e.f,1,l_grid_f) + repmat(grid_f,n_states,1)) - zt1.*qt1.*bft1);
ratio_f = Euler_ratio(s_par,s_state,zt1,num_f,denom_f,'f');
euler_f = abs(p - ratio_f);

[fval, b_star] = min(euler_g + euler_f);


Bt = grid_f(b_star);
bft = grid_f(b_star);
    
end