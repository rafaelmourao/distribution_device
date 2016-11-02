function [varargout] = Euler_sum_cor3(p,s_par,s_grid,s_state,s_gov,s_investors)

%This is the case where the FOREIGN investor chooses zero bonds:'bf = 0'

%% VARIABLES NEEDED

%Parameters
sigma = s_par.sigma;
lambda = s_par.lambda;
Qr = s_par.Qr;
tc = s_par.tc;
n_states = s_par.n_states;

%Grid
grid_r = s_grid.grid_b_r;
grid_r_g = grid_r + 1e-10;
l_grid_r = length(grid_r);


%State
n = s_state.n;
brt_1 = s_state.brt_1;
bgt_1 = s_state.bgt_1;
rt = s_state.rt;
lambdat = lambda(n);
Qrt = Qr(n);

%Investors
r0 = s_investors.r0;
q0 = s_investors.q0;
z0 = s_investors.z0;
br0 = s_investors.br0;
brt1 = squeeze(br0(:,1,:));

%Government
crt = s_gov.crt;
cwt = s_gov.cwt;
rt1 = squeeze(r0(:,:,1));
qt1 = squeeze(q0(:,:,1));
zt1 = squeeze(z0(:,:,1));
bg0 = s_gov.bg0;
bgt1 = squeeze(bg0(:,:,1));
cr0 = s_gov.cr0;
crt1 = squeeze(cr0(:,:,1));
cw0 = s_gov.cw0;
cwt1 = squeeze(cw0(:,:,1));

%% ALGORITHM

%JUST FOR YOU TO KNOW -> YOU'RE GOING TO HAVE PROBLEMS WITH
%THE 'zt1' IN ALL THE NUMERATORS: IT CAN BECOME ZERO AND
%IT'S RAISED TO THE POWER OF A NEGATIVE NUMBER. YOU MUST
%FIX THAT LATER

denom_g = tc*(crt+cwt) - bgt_1/lambdat + p*grid_r_g';

num_g = (zt1.^(-1/sigma.g)).*(((tc*lambda)*ones(1,l_grid_r)).*(crt1+cwt1) - ones(n_states,1)*grid_r_g +...
                (lambda*ones(1,l_grid_r)).*qt1.*zt1.*bgt1);

ratio_g = Euler_ratio(s_par,s_state,num_g,denom_g,'g');

euler_g = abs(ratio_g - p);


denom_r = ((1+rt).^((sigma.r-1)/sigma.r)).*(Qrt + brt_1/lambdat - p*grid_r');

num_r = ((zt1.^(-1/sigma.r)).*((1+rt1).^((sigma.r-1)/sigma.r))).*...
        ((lambda.*Qr)*ones(1,l_grid_r) + ones(n_states,1)*grid_r - (lambda*ones(1,l_grid_r)).*qt1.*zt1.*brt1);
            
ratio_r = Euler_ratio(s_par,s_state,num_r,denom_r,'r');

euler_r = abs(p - ratio_r);

[varargout{1}, b_star] = min(euler_g + euler_r);

Bt = grid_r_g(b_star);

brt = grid_r(b_star);

if nargout == 2
    
    varargout{1} = Bt;
    varargout{2} = brt;
    
    return
    
end

end