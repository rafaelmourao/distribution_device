function [fval, Bt, brt] = Euler_sum_cor3(p,s_par,s_grid,s_state,s_gov,s_investors)

%This is the case where the FOREIGN investor chooses zero bonds:'bf = 0'

%% VARIABLES NEEDED

%Parameters
sigma = s_par.sigma;
tc = s_par.tc;
n_states = s_par.n_states;

%Grid
grid_r = s_grid.grid_b_r;

%State
brt_1 = s_state.brt_1;
bgt_1 = s_state.bgt_1;
rt = s_state.rt;
wt = s_state.wt;

%Investors
r0 = s_investors.r0;
w0 = s_investors.w0;
q0 = s_investors.q0;
z0 = s_investors.z0;
br0 = s_investors.br0;
brt1 = squeeze(br0(:,1,:));

%Government
crt = s_gov.crt;
cft = s_gov.cft;
rt1 = r0(:,:,1);
wt1 = w0(:,:,1);
qt1 = q0(:,:,1);
zt1 = z0(:,:,1);
bgt1 = s_gov.bg0(:,:,1);
crt1 = s_gov.cr0(:,:,1);
cft1 = s_gov.cf0(:,:,1);

%% ALGORITHM

denom_g = tc*(crt+cft) - bgt_1 + p*grid_r';
num_g = ((tc*(crt1+cft1) - repmat(grid_r,n_states,1) +...
                zt1.*qt1.*bgt1));
ratio_g = Euler_ratio(s_par,s_state,zt1,num_g,denom_g,'g');
euler_g = abs(ratio_g - p);


denom_r = ((1+rt)*brt_1 + wt - p*grid_r');
num_r = ((1+rt1).^(-1/sigma.r)).*...
        (bsxfun(@times,(1+rt1),grid_r) + wt1 - zt1.*qt1.*brt1);
ratio_r = Euler_ratio(s_par,s_state,zt1,num_r,denom_r,'r');
euler_r = abs(p - ratio_r);

[fval, b_star] = min(euler_g + euler_r);


Bt = grid_r(b_star);
brt = grid_r(b_star);

end