function [status, p, br, bf, B] = Corner_Solution(s_par,s_grid,s_state,s_gov,s_investors)

%Each output variable has '6' lines, each refering to one of the '6'
%possible cases of corner solutions.

%% VARIABLES NEEDED

%Parameters
epsilon = .1;
sigma = s_par.sigma;
tc = s_par.tc;
e = s_par.e;

%Grid
grid_r = s_grid.grid_b_r;
grid_f = s_grid.grid_b_f;

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
bgt1 = s_gov.bg0;
crt1 = s_gov.cr0;
cft1 = s_gov.cf0;

%Investors
r0 = s_investors.r0;
w0 = s_investors.w0;
q0 = s_investors.q0;
z0 = s_investors.z0;
br0 = s_investors.br0;
bf0 = s_investors.bf0;


%% CORNER SOLUTIONS: CASES

%%%%% CASE 1: B = 0 %%%%%
B_c1 = 0;
br_c1 = 0;
bf_c1 = 1e-10;
rt1_c1 = r0(:,1,1);
wt1_c1 = w0(:,1,1);
qt1_c1 = q0(:,1,1);
zt1_c1 = z0(:,1,1);
brt1_c1 = br0(:,1,1);
bft1_c1 = bf0(:,1,1);

denom_r_c1 = ((1+rt)*brt_1 + wt);
num_r_c1 = ((1+rt1_c1).^(-1/sigma.r)).*...
                (wt1_c1 - zt1_c1.*qt1_c1.*brt1_c1);
            
pr_0_c1 = Euler_ratio(s_par,s_state,zt1_c1,num_r_c1,denom_r_c1,'r');

denom_f_c1 = (1+rt)*(bft_1 + eft);
num_f_c1 = ((1+rt1_c1).^(-1/sigma.f)).*...
                (bsxfun(@times,(1+rt1_c1),e.f) - zt1_c1.*qt1_c1.*bft1_c1);

pf_0_c1 = Euler_ratio(s_par,s_state,zt1_c1,num_f_c1,denom_f_c1,'f');


crt1_c1 = crt1(:,1,1);
cft1_c1 = cft1(:,1,1);
bgt1_c1 = bgt1(:,1,1);

denom_g_c1 = tc*(crt+cft) - bgt_1;
num_g_c1 = (tc*(crt1_c1 + cft1_c1) +...
                zt1_c1.*qt1_c1.*bgt1_c1);
            
pg_0_c1 = Euler_ratio(s_par,s_state,zt1_c1,num_g_c1,denom_g_c1,'g');

p_c1 = max(pr_0_c1,pf_0_c1);

status1 = p_c1 < pg_0_c1;


%%%%% CASE 2: br = 0 %%%%%
br_c2 = 0;

Obj_func2 = @(p) Euler_sum_cor2(p,s_par,s_grid,s_state,s_gov,s_investors);

status2 = 1;
p_max = 10;

[p_c2, sum_euler2] = fminbnd(Obj_func2,0,p_max);

[~, B_c2, bf_c2] = Euler_sum_cor2(p_c2,s_par,s_grid,s_state,s_gov,s_investors);
id_bf_c2 = find(grid_f == bf_c2);

rt1_c2 = r0(:,1,id_bf_c2);
wt1_c2 = w0(:,1,id_bf_c2);
qt1_c2 = q0(:,1,id_bf_c2);
zt1_c2 = z0(:,1,id_bf_c2);
brt1_c2 = br0(:,1,id_bf_c2);

denom_r_c2 = (1+rt)*brt_1 + wt;
num_r_c2 = ((1+rt1_c2).^(-1/sigma.r)).*...
                (wt1_c2 - zt1_c2.*qt1_c2.*brt1_c2);
            
pr_0_c2 = Euler_ratio(s_par,s_state,zt1_c2,num_r_c2,denom_r_c2,'r');

if sum_euler2 > epsilon || pr_0_c2 > p_c2
    
    status2 = 0;
    p_c2 = 0;
    B_c2 = 1e-10;
    bf_c2 = 1e-10;
    
end

%%%%% CASE 3: bf = 0 %%%%%
bf_c3 = 0;

Obj_func3 = @(p) Euler_sum_cor3(p,s_par,s_grid,s_state,s_gov,s_investors);

status3 = 1;
p_max = 10;

[p_c3, sum_euler3] = fminbnd(Obj_func3,0,p_max);

[~, B_c3, br_c3] = Euler_sum_cor3(p_c3,s_par,s_grid,s_state,s_gov,s_investors);

id_br_c3 = find(grid_r == br_c3);
rt1_c3 = r0(:,id_br_c3,1);
qt1_c3 = q0(:,id_br_c3,1);
zt1_c3 = z0(:,id_br_c3,1);
bft1_c3 = bf0(:,id_br_c3,1);

denom_f_c3 = (1+rt)*(eft + bft_1);
num_f_c3 = ((1+rt1_c3).^(-1/sigma.f)).*...
                (bsxfun(@times,(1+rt1_c3),e.f) - zt1_c3.*qt1_c3.*bft1_c3);
            
pf_0_c3 = Euler_ratio(s_par,s_state,zt1_c3,num_f_c3,denom_f_c3,'f');

if sum_euler3 > epsilon || pf_0_c3 > p_c3    
    status3 = 0;
    p_c3 = 0;
    B_c3 = 1e-10;
    br_c3 = 0;
end


%CASE 4: B = 0 and br = 0
B_c4 = 0;
br_c4 = 0;
bf_c4 = 0;  
pr_0_c4 = pr_0_c1;
pf_0_c4 = pf_0_c1;                        
pg_0_c4 = pg_0_c1;
p_c4 = pf_0_c4;
status4 = (p_c4 < pg_0_c4) && (p_c4 > pr_0_c4);


%CASE 5: B = 0 and bf = 0
B_c5 = 0;
br_c5 = 0;
bf_c5 = 0;
pr_0_c5 = pr_0_c1;
pf_0_c5 = pf_0_c1;                        
pg_0_c5 = pg_0_c1;
p_c5 = pr_0_c5;
status5 = (p_c5 < pg_0_c5) && (p_c5 > pf_0_c5);


%CASE 6: br = 0 and bf = 0
B_c6 = 0;
br_c6 = 0;
bf_c6 = 0;
pr_0_c6 = pr_0_c1;
pf_0_c6 = pf_0_c1;                        
pg_0_c6 = pg_0_c1;
p_c6 = pg_0_c6;
status6 = (p_c6 > pr_0_c6) && (p_c6 > pf_0_c6);

%% FINALLY

status = [status1;status2;status3;status4;status5;status6];
p = [p_c1;p_c2;p_c3;p_c4;p_c5;p_c6];
B = [B_c1;B_c2;B_c3;B_c4;B_c5;B_c6];
br = [br_c1;br_c2;br_c3;br_c4;br_c5;br_c6];
bf = [bf_c1;bf_c2;bf_c3;bf_c4;bf_c5;bf_c6];

end