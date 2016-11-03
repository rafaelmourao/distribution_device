function [status, p, br, bf, B] = Interior_Solution(s_par,s_grid,s_state,s_gov,s_investors)
%% VARIABLES NEEDED

%Parameters
epsilon = .1;

%% SOLUTION

Obj_func = @(p) Euler_sum_int(p,s_par,s_grid,s_state,s_gov,s_investors);

status = 1;

p_max = 4;

%grid_p = 0:.001:p_max;

%sum_euler_func = arrayfun(Obj_func,grid_p);

%[sum_euler, i] = min(sum_euler_func);

%p = grid_p(i);
[p sum_euler] = fminbnd(Obj_func,0,p_max);
%if p > 2.5
    
%    grid_p = 0:.01:p_max;
%    sum_euler_func = arrayfun(Obj_func,grid_p);
%    plot(grid_p,sum_euler_func);

%end

%bft_1 = s_state.bft_1;

%if bft_1 > .7
    
%    problem_here = 1;
%    plot(0:.001:p_max,helpgod);
    
%end

[B br bf] = Euler_sum_int(p,s_par,s_grid,s_state,s_gov,s_investors);

if sum_euler > epsilon
    
    status = 0;
    %p = 0;
    %B = 1e-10;
    %br = 0;
    %bf = 1e-10;
    
end

end
