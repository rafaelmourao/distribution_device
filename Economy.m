classdef Economy
    properties
        % Necessary economy parameters
        beta %Intertemporal discount rate
        sigma % Utility function parameter: risk aversion (structure)
        phi %Probability of redemption (Arellano)
        lambda %Government preference parameter: foreigners relative to residents
        tc %Tax rate over CONSUMPTION
        A % Fixed income stream for the government
        alpha %Participation of capital on production
        rho % Elasticity of Substitution between capital and labor (=1/(1-rho))
        e % shocks to the consumers endowment
        prob % transition matrix
        n_bonds % number of bonds for sale
        
        % Dependent variables
        n_states % number of states
        grid % structure that defines the grid of bonds, optional
        default % default outcomes (as below)
        
        % Iteration variables
        Vd % Value function in case of default
        Vnd % Value function in case of no default
        Vo % Value function
        Wnd % Welfare in case of no default
        cr % Resident consumption policy function in case of no default
        cf % Foreign consumption policy function in case of no default
        kr % Resident capital policy function in case of no default
        kf % Foreign capital policy function in case of no default
        br % Resident bond policy function holdings in case of no default
        bf % Foreign bond policy function holdings in case of no default
        bg % Government bond policy function holdings in case of no default
        z % Government default policy function
        r % Interest rates without default
        w % Wages without default
        q % Prices without default
        g % Government spending without default
    end
    
    properties (Hidden = true)
        % Extended arrays to conform with policy function dimensions
        extended_grid
        extended_default_r
        extended_default_w
    end
    
    methods
        
        function obj = Economy(param)
            % Setting parameters
            obj.beta = param.beta;
            obj.sigma = param.sigma;
            obj.phi = param.phi;
            obj.lambda = param.lambda;
            obj.tc = param.tc;
            obj.A = param.A;
            obj.alpha = param.alpha;
            obj.rho = param.rho;
            obj.e = param.e;
            obj.prob = param.prob;
            obj.n_bonds = param.n_bonds;
            
            % Setting dependent properties
            
            obj.n_states = size(obj.prob,1);        %Numbers of States of Nature
            
            % Constructing grid
            
            if isfield(param,'grid') % if grid is supplied
                obj.grid = param.grid;
            else
                obj.grid.b_r = linspace(param.min_b,param.max_b,obj.n_bonds);            %Grid for resident bonds:
                obj.grid.b_f = obj.grid.b_r;                                 %Grid for foreigner bonds;
                obj.grid.r_aux = repmat(obj.grid.b_r,1,obj.n_bonds);             %Combinations of bonds
                obj.grid.f_aux = kron(obj.grid.b_f,ones(1,obj.n_bonds));
                obj.grid.b_g = obj.grid.r_aux + obj.grid.f_aux;
            end
            
            % Generating auxiliary grid
            
            obj.extended_grid.b_r = repmat(obj.grid.b_r,[obj.n_states,1,obj.n_bonds]);
            obj.extended_grid.b_f = repmat(reshape(obj.grid.b_f,[1,1,obj.n_bonds]),...
                obj.n_states,obj.n_bonds);
            obj.extended_grid.e_f = repmat(obj.e.f,[1,obj.n_bonds,obj.n_bonds]);
            
            % Value functions
            
            obj.Vd = zeros(obj.n_states,1);                                 %Value function in case of DEFAULT
            obj.Vnd = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);          %Value function in case of NO DEFAULT
            obj.Vo = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);           %Value function for the OPTION of DEFAULT
            
            %% POLICY FUNCTIONS
            
            %Consumers
            
            obj.cr = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);          %CONSUMPTION policy funtion for RESIDENTS
            obj.cf = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);          %CONSUMPTION policy funtion for FOREIGNERS
            
            obj.kr = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);          %CAPITAL policy funtion for RESIDENTS
            obj.kf = repmat(obj.e.f,[1,obj.n_bonds,obj.n_bonds]);           %CAPITAL policy funtion for FOREIGNERS
            
            %obj.br = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);          %BONDS policy funtion for RESIDENTS
            obj.br = reshape(kron(obj.grid.r_aux,ones(1,obj.n_states)),...
                [obj.n_states,obj.n_bonds,obj.n_bonds]);
            %             obj.bf = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);          %BONDS policy funtion for FOREIGNERS
            obj.bf = reshape(kron(obj.grid.f_aux,ones(1,obj.n_states)),...
                [obj.n_states,obj.n_bonds,obj.n_bonds]);
            
            %Government
            %             obj.bg = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);     %     %BONDS policy funtion for the GOVERNMENT
            obj.bg = obj.br + obj.bf;
            obj.g = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);            %PUBLIC EXPENDITURE policy funtion for the GOVERNMENT
            obj.z = ones(obj.n_states,obj.n_bonds,obj.n_bonds);       %DEFAULT policy funtion for the GOVERNMENT
            
            %% PRICE FUNCTIONS
            
            obj.r = obj.alpha*((obj.kr+obj.kf).^(obj.rho-1)).*...
                ((obj.alpha*((obj.kr+obj.kf).^obj.rho) +...
                (1-obj.alpha)).^(1/obj.rho-1));     %Interest Rate
            
            obj.w = (1-obj.alpha)*((obj.alpha*((obj.kr+obj.kf).^obj.rho)...
                + (1-obj.alpha)).^(1/obj.rho-1));           %Wage
            
            obj.q = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);           %Price of Public Bond
            
            %% Default Outcomes
            
            obj.default.r = zeros(obj.n_states,1); % Variables in case of default
            obj.default.w = zeros(obj.n_states,1);
            obj.default.cr = zeros(obj.n_states,1);
            obj.default.cf = zeros(obj.n_states,1);
            obj.default.g = zeros(obj.n_states,1);
            obj.default.W = zeros(obj.n_states,1);
            for n = 1:obj.n_states
                obj.default.r(n) = obj.alpha*(obj.e.f(n).^(obj.rho-1)).*...
                    ((obj.alpha*(obj.e.f(n).^obj.rho) + (1-obj.alpha)).^(1/obj.rho-1));
                obj.default.w(n) = (1-obj.alpha)*((obj.alpha*(obj.e.f(n))^(obj.rho)) + (1-obj.alpha))^(1/obj.rho-1);
                obj.default.cr(n) = (1/(1+obj.tc))*obj.default.w(n);
                obj.default.cf(n) = (1/(1+obj.tc))*(1+obj.default.r(n))*obj.e.f(n);
                obj.default.g(n) = obj.A + obj.tc*(obj.default.cr(n) + obj.default.cf(n));
                obj.default.W(n) =  utility_function(obj.default.cr(n),obj.sigma.r) +...
                    obj.lambda*utility_function(obj.default.cf(n),obj.sigma.f) ...
                    +  utility_function(obj.default.g(n),obj.sigma.g);
                obj.extended_default_r = repmat(obj.default.r,...
                    1,obj.n_bonds*obj.n_bonds);
                obj.extended_default_w = repmat(obj.default.w,...
                    1,obj.n_bonds*obj.n_bonds);
            end
        end
        
        function obj = update(obj, nworkers)
            if ( ~exist('nworkers','var') || isempty(nworkers) )
                nworkers = 12;
            end
            
            n_b_states = obj.n_states*obj.n_bonds^2;
            n_b_states_size = [obj.n_states, obj.n_bonds, obj.n_bonds];
            
            p = zeros(n_b_states_size);
            br_s = zeros(n_b_states_size);
            bf_s = zeros(n_b_states_size);
            bg_s = zeros(n_b_states_size);
            next_Vo = zeros(n_b_states_size);
            
            
            parfor(i=1:n_b_states,nworkers)
                
                [n, id_br, id_bf] = ind2sub(n_b_states_size,i);
                
                [p(i), br_s(i),...
                    bf_s(i), bg_s(i)] = ...
                    Solution(obj, n, id_br, id_bf);
                
                loc_br_s = ( obj.grid.b_r == br_s(i) );
                loc_bf_s = ( obj.grid.b_f == bf_s(i) );
                % Expected value for next period
                next_Vo(i) = obj.prob(n,:)*obj.Vo(:,loc_br_s,loc_bf_s);
            end
            
            obj.q = p;               %Equilibrium price for the Bonds' market
            obj.br = br_s;                          %RESIDENTS demand for Bonds
            obj.bf = bf_s;                          %FOREIGNERS demand for Bonds
            obj.bg = bg_s;                          %GOVERNEMNT supply of Bonds
            
            obj.kr = obj.extended_grid.b_r;                       %RESIDENTS capital supply
            obj.kf = obj.extended_grid.e_f +...
                obj.extended_grid.b_f;                           %FOREIGNERS capital supply
            
            obj.r = ...
                obj.alpha*((obj.kr+obj.kf).^(obj.rho-1)).*...
                ((obj.alpha*((obj.kr+obj.kf).^obj.rho) +...
                (1-obj.alpha)).^(1/obj.rho-1));
            
            obj.w = ...
                (1-obj.alpha)*(obj.alpha*((obj.kr+obj.kf).^(obj.rho))...
                + (1-obj.alpha)).^(1/obj.rho-1);
            
            obj.cr = (1/(1+obj.tc))*((1+obj.r).*obj.kr...
                + obj.w - obj.q.*obj.br);
            obj.cr(obj.cr<0) = 0;
            
            obj.cf = (1/(1+obj.tc))*((1+obj.r).*obj.kf...
                - obj.q.*obj.bf);
            obj.cf(obj.cf<0) = 0;
            
            obj.g = obj.A + obj.tc*(obj.cr + obj.cf) + obj.q.*obj.bg - ...
                (obj.extended_grid.b_r + obj.extended_grid.b_f);
            obj.g(obj.g<0) = 0;
            
            obj.Wnd = utility_function(obj.cr,obj.sigma.r) + ...
                obj.lambda*utility_function(obj.cf,obj.sigma.f) + ...
                utility_function(obj.g,obj.sigma.g);
            
            obj.Vnd = obj.Wnd + obj.beta*next_Vo;
            
            obj.Vd = obj.default.W + obj.beta*obj.prob*...
                (obj.phi * obj.Vo(:,1,1) + (1-obj.phi) * obj.Vd);
            
            % in case of no default
            obj.z = ones(obj.n_states,obj.n_bonds,obj.n_bonds);
            obj.Vo = obj.Vnd;
            
            % check where there is default
            def = bsxfun(@lt,obj.Vnd,obj.Vd);
            [def_states, ~] = find(def); % retrieving the states of all default occ.
            obj.z(def) = 0; % updating default decision
            obj.Vo(def) = obj.Vd(def_states); % updating policy function
        end
        
        function [p, br_s, bf_s, bg_s] = Solution(obj, n, id_br, id_bf)
            % This function solves the government and consumers
            % optimization problem. It selects the highest total quantity
            % bond portfolio where the highest price both consumers accept
            % to buy is lower than the lowest price the government accepts.
            
            % Default (corner) solution
            p = 0;
            br_s = 0;
            bf_s = 0;
            bg_s = 0;
            
            eq_price_g = eq_prices_government(obj, n, id_br, id_bf);
            eq_price_r = eq_prices_residents(obj, n, id_br, id_bf);
            eq_price_f = eq_prices_foreigners(obj, n, id_br, id_bf);
            
            prices = min(eq_price_f,eq_price_r);
            prices(eq_price_g > prices) = 0;
            
            [~, sorted_bonds] = sort(obj.grid.b_g,2,'descend');
            for i = sorted_bonds
                % Accept solution if price is positive and finite
                % If none found, solution is the corner (0,0) solution with
                % zero price
                if ( prices(i) > 0 && isfinite(prices(i)))
                    p = prices(i);
                    br_s = obj.grid.r_aux(i);
                    bf_s = obj.grid.f_aux(i);
                    bg_s = obj.grid.b_g(i);
                    break
                end
            end
        end
        
        function eq_price_g = eq_prices_government(obj, n, id_br, id_bf)
            % Calculation of the minimal positive price the government may accept
            % to sell the combination of bonds through a bisection
            % procedure. Equals infinity for the cases where the government
            % wouldn't accept any price (and 0,0)
            
            % Variables
            
            % Past and present
            probt = obj.prob(n,:);
            brt_1 = obj.grid.b_r(id_br);
            bft_1 = obj.grid.b_f(id_bf);
            bgt_1 = brt_1 + bft_1;
            rt = obj.r(n, id_br, id_bf);
            wt = obj.w(n, id_br, id_bf);
            eft = obj.e.f(n);
            
            % Future
            rt1 = obj.r(:,:);
            wt1 = obj.w(:,:);
            qt1 = obj.q(:,:);
            zt1 = obj.z(:,:);
            brt1 = obj.br(:,:);
            bft1 = obj.bf(:,:);
            bgt1 = obj.bg(:,:);
            
            % In case of default, future interest rate and wages are the
            % default ones
            rt1(~zt1) = obj.extended_default_r(~zt1);
            wt1(~zt1) = obj.extended_default_w(~zt1);
            
            grid_r = obj.grid.r_aux;
            grid_f = obj.grid.f_aux;
            grid_g = obj.grid.b_g;
            l_grid_g = length(grid_g);
            
            
            num_g = obj.A + (obj.tc/(1+obj.tc))*...
                ((zt1.*bsxfun(@times,(1+rt1),grid_r) + ...
                wt1 - zt1.*qt1.*brt1) + ...
                (1+rt1).*(repmat(obj.e.f,1,l_grid_g) + ...
                bsxfun(@times,zt1,grid_f)) - zt1.*qt1.*bft1) + ...
                zt1.*(qt1.*bgt1 - repmat(grid_g,obj.n_states,1));
            
            % Disregard cases where the numerator is negative independently
            % of prices, and the (0,0) portfolio, which is the default
            % solution
            
            valid_g = all(num_g > 0);
            valid_g(1) = 0;
            
            denom_g_0 = obj.A + (obj.tc/(1+obj.tc))*...
                (((1+rt)*brt_1 + wt) + ...
                ((1+rt)*(eft + bft_1))) - ...
                bgt_1;
            
            grid_g_valid = grid_g(valid_g);
            euler_denom_g_valid = @(p) (denom_g_0 + ...
                (1/(1+obj.tc))*p.*grid_g_valid).^-obj.sigma.g;
            
            euler_num_g = obj.beta*(probt*(zt1.*(num_g.^-obj.sigma.g)));
            euler_num_g(any(num_g)<0) = NaN; % just in case
            euler_num_g_valid = euler_num_g(valid_g);
            
            ratio_g = @(p) euler_num_g_valid ./ euler_denom_g_valid(p);
            
            
            % Find price where euler denominator approaches zero
            min_feasible_price_g = ( -obj.A - obj.tc*...
                ((1+rt)*brt_1 + wt + ...
                (1+rt)*(eft + bft_1)) + ...
                (1+obj.tc)*bgt_1 ) ./ grid_g;
            
            % Find a price where Euler ratio is below price or set a
            % maximum of 1e4;
            
            pmax = max(min_feasible_price_g(valid_g)+.5,.5);
            while any((ratio_g(pmax) > pmax) & pmax < 100)
                pmax = 5*pmax;
            end
            
            % For invalid cases, set government equilibrium prices to Inf,
            % otherwise find prices which equate euler ratio to the price
            
            eq_price_g = Inf*ones(size(valid_g));
            if any(valid_g)
                eq_price_g(valid_g) = bisection(@(p) ratio_g(p) - p,...
                    min_feasible_price_g(valid_g), pmax);
                eq_price_g(isnan(eq_price_g)) = Inf;
            end
            
        end
        
        function eq_price_r = eq_prices_residents(obj, n, id_br, id_bf)
            % Calculation of the maximal prices the residents may accept to
            % buy the bond portfolio. Equals infinity if resident doesn't
            % buy anything, and zero is she accepts any positive price.
            
            % Variables
            grid_r = obj.grid.r_aux;
            
            % Past and present
            probt = obj.prob(n,:);
            brt_1 = obj.grid.b_r(id_br);
            rt = obj.r(n, id_br, id_bf);
            wt = obj.w(n, id_br, id_bf);
            
            % Future
            rt1 = obj.r(:,:);
            wt1 = obj.w(:,:);
            qt1 = obj.q(:,:);
            zt1 = obj.z(:,:);
            brt1 = obj.br(:,:);
            
            % In case of default, future interest rate and wages are the
            % default ones
            rt1(~zt1) = obj.extended_default_r(~zt1);
            wt1(~zt1) = obj.extended_default_w(~zt1);
            
            num_r = ((1+rt1).^(-1/obj.sigma.r)).*...
                (zt1.*bsxfun(@times,(1+rt1),grid_r) + wt1 - zt1.*qt1.*brt1);
            
            euler_num_r = obj.beta*(probt*(zt1.*(num_r.^-obj.sigma.r)));
            euler_num_r(any(num_r)<0) = NaN;
            
            denom_r_0 = (1+rt)*brt_1 + wt;
            ratio_r_0 = euler_num_r ./ (denom_r_0.^-obj.sigma.r);
            
            % Disregard cases where the ratio is not well defined
            % independent of prices, and all portfolios where the resident
            % buys no bonds
            
            valid_r = ~(isnan(ratio_r_0) | ratio_r_0 < 1e-6 | grid_r == 0);
            valid_r(1) = 0;
            
            grid_r_valid = grid_r(valid_r);
            euler_num_r_valid = euler_num_r(valid_r);
            
            euler_denom_r = @(p) (denom_r_0 - p.*grid_r_valid).^-obj.sigma.r;
            ratio_r = @(p) euler_num_r_valid ./ euler_denom_r(p);
            
            % max price where the denominator is still positive (ratio will
            % be always lower than price in these cases)
            max_feasible_price_r = ( (1+rt)*brt_1 + wt ) ./ grid_r ;
            
            eq_price_r = -Inf*ones(size(valid_r));
            if any(valid_r)
                eq_price_r(valid_r) = bisection(@(p) ratio_r(p) - abs(p),...
                    0,max_feasible_price_r(valid_r));
            end
            eq_price_r(grid_r == 0) = Inf;
        end
        
        function eq_price_f = eq_prices_foreigners(obj, n, id_br, id_bf)
            % Calculation of the maximal prices the foreigner may accept to
            % buy the bond portfolio. Equals infinity if foreigner doesn't
            % buy anything, and zero is she accepts any positive price.
            
            % Variables
            grid_f = obj.grid.f_aux;
            l_grid_g = length(grid_f);
            
            % Past and present
            probt = obj.prob(n,:);
            bft_1 = obj.grid.b_f(id_bf);
            rt = obj.r(n, id_br, id_bf);
            
            % Future
            rt1 = obj.r(:,:);
            qt1 = obj.q(:,:);
            zt1 = obj.z(:,:);
            bft1 = obj.br(:,:);
            eft = obj.e.f(n);
            
            % In case of default, future interest rate and wages are the
            % default ones
            rt1(~zt1) = obj.extended_default_r(~zt1);
            
            num_f = (1+rt1).^(-1/obj.sigma.f).*...
                ((1+rt1).*(repmat(obj.e.f,1,l_grid_g) +...
                bsxfun(@times,zt1,grid_f)) - zt1.*qt1.*bft1);
            
            euler_num_f = obj.beta*(probt*(zt1.*(num_f.^-obj.sigma.f)));
            euler_num_f(any(num_f)<0) = NaN;
            
            denom_f_0 = ((1+rt)*(eft + bft_1));
            ratio_f_0 = euler_num_f ./ (denom_f_0.^-obj.sigma.f);
            
            % Disregard cases where the ratio is not well defined
            % independent of prices, and all portfolios where the foreigner
            % buys no bonds
            valid_f = ~(isnan(ratio_f_0) | ratio_f_0 < 1e-6 | grid_f == 0);
            valid_f(1) = 0;
            
            grid_f_valid = grid_f(valid_f);
            euler_num_f_valid = euler_num_f(valid_f);
            
            euler_denom_f = @(p) (denom_f_0 - p.*grid_f_valid).^-obj.sigma.f;
            ratio_f = @(p) euler_num_f_valid./euler_denom_f(p);
            
            % max price where the denominator is still positive (ratio will
            % be always lower than price in these cases)
            max_feasible_price_f = denom_f_0 ./ grid_f;
            
            eq_price_f = -Inf*ones(size(valid_f));
            if any(valid_f)
                eq_price_f(valid_f) = bisection(@(p) ratio_f(p) - abs(p),...
                    0, max_feasible_price_f(valid_f));
            end
            eq_price_f(grid_f == 0) = Inf;
            
        end
        
    end
end

function f = utility_function(x, sigma)
f = (x.^(1-sigma) - 1)/(1 - sigma);
end