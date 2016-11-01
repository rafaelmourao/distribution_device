function [p br_s bf_s bg_s] = Solution(s_par,s_grid,s_state,s_gov,s_investors)

[status_int p_int br_int bf_int bg_int] = Interior_Solution(s_par,s_grid,s_state,s_gov,s_investors);

[status_cor p_cor br_cor bf_cor bg_cor] = Corner_Solution(s_par,s_grid,s_state,s_gov,s_investors);

status = [status_int;status_cor];
price = [p_int;p_cor];
br = [br_int;br_cor];
bf = [bf_int;bf_cor];
bg = [bg_int;bg_cor];

status
[~, i] = find(status == 1,1);

p = price(i);
br_s = br(i);
bf_s = bf(i);
bg_s = bg(i);

if isempty(i)
   
    p = p_int;
    br_s = br_int;
    bf_s = bf_int;
    bg_s = bg_int;
    
end

end