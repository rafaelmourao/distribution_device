function [varargout] = Euler_sum_cor2(p,s_par,s_grid,s_state,s_gov,s_investors)

%This is the case where the RESISENT investor chooses zero bonds:'br = 0'

%% VARIABLES NEEDED

%Parameters
sigma = s_par.sigma;
lambda = s_par.lambda;
Qf = s_par.Qf;
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
lambdat = lambda(n);
Qft = Qf(n);

%Investors
r0 = s_investors.r0;
q0 = s_investors.q0;
z0 = s_investors.z0;
bf0 = s_investors.bf0;
bft1 = squeeze(bf0(:,1,:));

%Government
crt = s_gov.crt;
cwt = s_gov.cwt;
rt1 = squeeze(r0(:,1,:));
qt1 = squeeze(q0(:,1,:));
zt1 = squeeze(z0(:,1,:));
bg0 = s_gov.bg0;
bgt1 = squeeze(bg0(:,1,:));
cr0 = s_gov.cr0;
crt1 = squeeze(cr0(:,1,:));
cw0 = s_gov.cw0;
cwt1 = squeeze(cw0(:,1,:));

%% ALGORITHM

%JUST FOR YOU TO KNOW -> YOU'RE GOING TO HAVE PROBLEMS WITH
%THE 'zt1' IN ALL THE NUMERATORS: IT CAN BECOME ZERO AND
%IT'S RAISED TO THE POWER OF A NEGATIVE NUMBER. YOU MUST FIX THAT LATER

denom_g = tc*(crt+cwt) - bgt_1/lambdat + p*grid_f';

num_g = (zt1.^(-1/sigma)).*(((tc*lambda)*ones(1,l_grid_f)).*(crt1+cwt1) - ones(n_states,1)*grid_f +...
                (lambda*ones(1,l_grid_f)).*qt1.*zt1.*bgt1);

ratio_g = Euler_ratio(s_par,s_state,num_g,denom_g);

euler_g = abs(ratio_g - p);
                

denom_f = ((1+rt).^((sigma-1)/sigma)).*(Qft + bft_1/lambdat - p*grid_f');

num_f = ((zt1.^(-1/sigma)).*((1+rt1).^((sigma-1)/sigma))).*...
        ((lambda.*Qf)*ones(1,l_grid_f) + ones(n_states,1)*grid_f - (lambda*ones(1,l_grid_f)).*qt1.*zt1.*bft1);
    
ratio_f = Euler_ratio(s_par,s_state,num_f,denom_f);

euler_f = abs(p - ratio_f);

[varargout{1}, b_star] = min(euler_g + euler_f);


Bt = grid_f(b_star);

bft = grid_f(b_star);

if nargout == 2
    
    varargout{1} = Bt;
    varargout{2} = bft;
    
    return
    
end

end