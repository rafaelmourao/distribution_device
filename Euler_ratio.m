function ratio = Euler_ratio(s_par,s_state,zt1,x,y, type)

% Important remark: the 'x' refers to the NUMERATOR while the 'y' refers to
% the DENOMINATOR. The former is inouted as a 'n_states x length_grid'
% matriz and the later as a '1 x length_grid' vector.

%% VARIABLES NEEDED

%Parameters
sigma = s_par.sigma.(type);
beta = s_par.beta;
prob = s_par.prob;
n_states = s_par.n_states;

%State
n = s_state.n;
probt = prob(n,:);

%% PROGRAM

%ERROR FUNTION
%This function identifies where there is a problem with consumption (it
%can never be negative), for both vectors: 'x' and 'y'. They will be part 
%of the numerator and denominator, respectively, of the euler equations. 
%The first part identifies the coordinates where the denominator has 
%negative values. The second part finds if there is any entry, for every
%state of nature in the future, where consumption might be negative; the
%whole column will be shut down.

error_c = ((y < 0) + max((x < 0))') > 0;

%FIX NUMERATOR FUNCTION
%This function corrects the numerator, placing 'Inf' whenever there is a
%negative consumption on either the numerator or denominator.

fix_num = (1./(ones(n_states,1)*(1 - error_c'))).*abs(zt1.*x);

%CALCULATING THE NUMERATOR
%This functions computes the numerator itself, but for the coordinates
%where there was negative consumption for either numerator or denominator,
%there will be '0' instead.

numerator = beta*(probt*(fix_num.^(-sigma)))';

%FIX DENOMINATOR FUNTION
%This function corrects the denominator, placing 'Inf' whenever there is
%a negative consumption on either the numerator or denominator.

fix_denom = (1./(1 - error_c)).*abs(y);

%CALCULATING THE DENOMINATOR
%This functions computes the denomitor itself, but for the coordinates 
%where there was negative consumption for either numerator or denominator,
%there will be '0' instead.

denominator = fix_denom.^(-sigma);

%RATIO FUNCTION
%This function calculates the MRS. But if there is any negative consumtion
%on any of the parts used to calculate this ratio, will show up a 'NaN' 
%(not a number) on that coordinate, for the '0/0' that would be calculated.
%The Matlab won't consider this coordinates with 'NaN' when trying to find 
%the policy functions.

ratio = numerator./denominator; 

end