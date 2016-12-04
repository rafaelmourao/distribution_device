classdef Economy
    properties
        % Necessary economy parameters
        has_default = true; % By default the model allows for default
        has_partial_default = true; % and partial default
        
        beta %Intertemporal discount rate
        sigma % Utility function parameter: risk aversion (structure)
        phi %Probability of redemption (Arellano)
        lambda %Government preference parameter: foreigners relative to residents
        theta %Discount factor over utility of public good
        tc %Tax rate over CONSUMPTION
        Ag = 0 % Fixed income stream for the government
        Ar = 0% Fixed income stream for the residents
        Af = 0% Fixed income stream for the foreginers
        discount % Discount factor summed on A's
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
        delta % Government open market policy function
        z % Government default policy function
        r % Interest rates without default
        w % Wages without default
        q % Prices without default
        g % Government spending without default
    end
    
    properties (Hidden = true)
        % Extended arrays to conform with policy function dimensions
        extended_grid
    end
    
    methods
        
        function obj = Economy(param)
            
            if isfield(param,'has_default') % if grid is supplied
                obj.has_default = param.has_default;
            end
            
            if isfield(param,'has_partial_default') % if grid is supplied
                obj.has_partial_default = param.has_partial_default;
            end

            % Setting parameters
            obj.beta = param.beta;
            obj.sigma = param.sigma;
            obj.phi = param.phi;
            obj.lambda = param.lambda;
            obj.theta = param.theta;
            obj.tc = param.tc;
            obj.Ag = param.Ag;
            obj.Ar = param.Ar;
            obj.Af = param.Af;
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
            obj.bg = obj.br + obj.bf; %BONDS policy funtion for the GOVERNMENT
            obj.g = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);            %PUBLIC EXPENDITURE policy funtion for the GOVERNMENT
            obj.z = ones(obj.n_states,obj.n_bonds,obj.n_bonds);       %DEFAULT policy funtion for the GOVERNMENT
            obj.delta = ones(obj.n_states,obj.n_bonds,obj.n_bonds); %MARKET openess policy function
            
            %% PRICE FUNCTIONS
            
            obj.r = obj.alpha*((obj.kr+obj.kf).^(obj.rho-1)).*...
                ((obj.alpha*((obj.kr+obj.kf).^obj.rho) +...
                (1-obj.alpha)).^(1/obj.rho-1));     %Interest Rate
            
            obj.w = (1-obj.alpha)*((obj.alpha*((obj.kr+obj.kf).^obj.rho)...
                + (1-obj.alpha)).^(1/obj.rho-1));           %Wage
            
            obj.q = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);           %Price of Public Bond
            
            obj = set_default(obj);
            
        end
        
        function obj = set_default(obj)
            %% Default Outcomes
            
            obj.default.z = ones(obj.n_states,obj.n_bonds,obj.n_bonds);
            obj.default.kr = ones(obj.n_states,obj.n_bonds,obj.n_bonds);
            obj.default.kf = ones(obj.n_states,obj.n_bonds,obj.n_bonds);
            obj.default.r = ones(obj.n_states,obj.n_bonds,obj.n_bonds);
            obj.default.w = ones(obj.n_states,obj.n_bonds,obj.n_bonds);
            obj.default.cr = ones(obj.n_states,obj.n_bonds,obj.n_bonds);
            obj.default.cr = ones(obj.n_states,obj.n_bonds,obj.n_bonds);
            obj.default.g = ones(obj.n_states,obj.n_bonds,obj.n_bonds);
            obj.default.W = ones(obj.n_states,obj.n_bonds,obj.n_bonds);
            
            for i = 1:numel(obj.z)
                if (i <= obj.n_states || ~obj.has_partial_default)
                    % "Full default" when Bg = 0
                    obj.default.z(i) = 0;
                else
                    options = optimset('TolX',1e-6);
                    ub = 1;
                    while ~isfinite(default_calc(.8*ub,i))
                        ub = .8*ub;
                    end
                    obj.default.z(i) = fminbnd(@(x) ...
                        -default_calc(x,i),0,ub,options);
                    if (obj.default.z(i) < 1e-5)
                        obj.default.z(i) = 0;
                    elseif (obj.default.z(i) > 1 - 1e-5)
                        obj.default.z(i) = 1;
                    end
                end
                [~, def] = default_calc(obj.default.z(i),i);
                obj.default.kr(i) = def.kr;
                obj.default.kf(i) = def.kf;
                obj.default.r(i) = def.r;
                obj.default.w(i) = def.w;
                obj.default.cr(i) = def.cr;
                obj.default.cf(i) = def.cf;
                obj.default.g(i) = def.g;
                obj.default.W(i) = def.W;
            end
            
            
            function [fval, def] = default_calc(z,i)
                [n, ~, ~] = ind2sub(size(obj.z),i);
                def.kr = z.*obj.extended_grid.b_r(i);
                def.kf = obj.e.f(n) + ...
                    z.*obj.extended_grid.b_f(i);
                def.r = obj.alpha*((def.kf+def.kr)^(obj.rho-1)).*...
                    ((obj.alpha*((def.kf+def.kr)^obj.rho) +...
                    (1-obj.alpha)).^(1/obj.rho-1));
                def.w = (1-obj.alpha)*((obj.alpha*((def.kf+def.kr)).^(obj.rho)) ...
                    + (1-obj.alpha))^(1/obj.rho-1);
                def.cr = obj.Ar + (1/(1+obj.tc))*(def.w +(1+def.r).*def.kr);
                def.cf = obj.Af + (1/(1+obj.tc))*(1+def.r).*def.kf;
                def.g = max(obj.Ag + obj.tc*(def.cr + def.cf) - z.*...
                    (obj.extended_grid.b_r(i) + obj.extended_grid.b_f(i)),0);
                def.W =  obj.lambda*utility_function(def.cr,obj.sigma.r) +...
                    (1-obj.lambda)*utility_function(def.cf,obj.sigma.f) ...
                    +  obj.theta*utility_function(def.g,obj.sigma.g);
                fval = def.W;
            end
            
        end
        
        function obj = update(obj, nworkers)
            if ( ~exist('nworkers','var') || isempty(nworkers) )
                nworkers = 12;
            end
            
            n_tot_bonds = obj.n_bonds^2;
            tot_bonds_size = [obj.n_bonds, obj.n_bonds];
            n_tot_bonds_size = [obj.n_states, obj.n_bonds, obj.n_bonds];
            
            p_s = zeros(n_tot_bonds_size);
            br_s = zeros(n_tot_bonds_size);
            bf_s = zeros(n_tot_bonds_size);
            bg_s = zeros(n_tot_bonds_size);
            next_Vo = zeros(n_tot_bonds_size);
            
            
            parfor(i=1:n_tot_bonds,nworkers)
                
                [id_br, id_bf] = ind2sub(tot_bonds_size,i);
                
                [p_s(:,i), id_s, br_s(:,i),...
                    bf_s(:,i), bg_s(:,i)] = ...
                    Solution(obj, id_br, id_bf);
                
                next_Vo(:,i) = diag(obj.prob*obj.Vo(:,id_s));
            end
            
            obj.q = p_s;               %Equilibrium price for the Bonds' market
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
            
            obj.cr = obj.Ar + (1/(1+obj.tc))*((1+obj.r).*obj.kr...
                + obj.w - obj.q.*obj.br);
            obj.cr(obj.cr<0) = 0;
            
            obj.cf = obj.Af + (1/(1+obj.tc))*((1+obj.r).*obj.kf...
                - obj.q.*obj.bf);
            obj.cf(obj.cf<0) = 0;
            
            obj.g = obj.Ag + obj.tc*(obj.cr + obj.cf) + obj.q.*obj.bg - ...
                (obj.extended_grid.b_r + obj.extended_grid.b_f);
            obj.g(obj.g<0) = 0;
            
            obj.Wnd = obj.lambda*utility_function(obj.cr,obj.sigma.r) + ...
                (1-obj.lambda)*utility_function(obj.cf,obj.sigma.f) + ...
                obj.theta*utility_function(obj.g,obj.sigma.g);
            
            obj.Vnd = obj.Wnd + obj.beta*next_Vo;
            
            obj.Vd = obj.default.W + repmat(obj.beta*obj.prob*...
                (obj.phi * obj.Vo(:,1,1) + (1-obj.phi) * obj.Vd(:,1,1)),...
                [1,obj.n_bonds,obj.n_bonds]);
            
            % in case of no default
            obj.z = ones(obj.n_states,obj.n_bonds,obj.n_bonds);
            obj.delta = obj.z;
            obj.Vo = obj.Vnd;
            
            % check where there is default
            if obj.has_default
                def = ( obj.Vnd < obj.Vd );
                obj.z(def) = obj.default.z(def) ; % updating default decision
                obj.delta(def) = 0;
                obj.Vo(def) = obj.Vd(def); % updating policy function
            end
        end
        
        function [p_s, id_s, br_s, bf_s, bg_s] = Solution(obj, id_br, id_bf)
            % This function solves the government and consumers
            % optimization problem. It selects the highest total quantity
            % bond portfolio where the highest price both consumers accept
            % to buy is lower than the lowest price the government accepts.
            
            % Default (corner) solution
            p_s = zeros(obj.n_states,1);
            id_s = ones(obj.n_states,1);
            br_s = zeros(obj.n_states,1);
            bf_s = zeros(obj.n_states,1);
            bg_s = zeros(obj.n_states,1);
            
            eq_price_g = eq_prices_government(obj, id_br, id_bf);
            eq_price_r = eq_prices_residents(obj, id_br, id_bf);
            eq_price_f = eq_prices_foreigners(obj, id_br, id_bf);
            
            prices = min(eq_price_f,eq_price_r);
            prices(eq_price_g > prices) = 0;
            
            [~, sorted_bonds] = sort(obj.grid.b_g,2,'descend');
            for n = 1:obj.n_states
                for i = sorted_bonds
                    % Accept solution if price is positive and finite
                    % If none found, solution is the corner (0,0) solution with
                    % zero price
                    if ( prices(n,i) > 0 && isfinite(prices(n,i)))
                        id_s(n) = i;
                        p_s(n) = prices(n,i);
                        br_s(n) = obj.grid.r_aux(i);
                        bf_s(n) = obj.grid.f_aux(i);
                        bg_s(n) = obj.grid.b_g(i);
                        break
                    end
                end
            end
        end
        
        function eq_price_g = eq_prices_government(obj, id_br, id_bf)
            % Calculation of the minimal positive price the government may accept
            % to sell the combination of bonds through a bisection
            % procedure. Equals infinity for the cases where the government
            % wouldn't accept any price (and 0,0)
            
            % Variables
            
            % Past and present
            brt_1 = obj.grid.b_r(id_br);
            bft_1 = obj.grid.b_f(id_bf);
            bgt_1 = brt_1 + bft_1;
            rt = obj.r(:,id_br, id_bf);
            wt = obj.w(:,id_br, id_bf);
            eft = obj.e.f;
            
            % Future
            delta1 = obj.delta(:,:);
            rt1 = obj.r(:,:);
            wt1 = obj.w(:,:);
            qt1 = obj.q(:,:);
            zt1 = obj.z(:,:);
            brt1 = obj.br(:,:);
            bft1 = obj.bf(:,:);
            bgt1 = obj.bg(:,:);
            
            % In case of default, future interest rate and wages are the
            % default ones
            if obj.has_default
                rt1(~delta1) = obj.default.r(~delta1);
                wt1(~delta1) = obj.default.w(~delta1);
            end
            
            grid_r = obj.grid.r_aux;
            grid_f = obj.grid.f_aux;
            grid_g = obj.grid.b_g;
            l_grid_g = length(grid_g);
            
            
            num_g = obj.Ag + (obj.tc/(1+obj.tc))*...
                ((zt1.*bsxfun(@times,(1+rt1),grid_r) + ...
                wt1 - delta1.*qt1.*brt1) + ...
                (1+rt1).*(repmat(obj.e.f,1,l_grid_g) + ...
                bsxfun(@times,zt1,grid_f)) - delta1.*qt1.*bft1) + ...
                zt1.*(qt1.*bgt1 - repmat(grid_g,obj.n_states,1));
            
            % Disregard cases where the numerator is negative independently
            % of prices, and the (0,0) portfolio, which is the default
            % solution
            
            valid_g = all(num_g > 0);
            valid_g(1) = 0;
            valid_g = repmat(valid_g,[obj.n_states,1]);
            
            denom_g_0 = obj.Ag + (obj.tc/(1+obj.tc))*...
                ((1+rt).*brt_1 + wt + ...
                (1+rt).*(eft + bft_1)) - ...
                bgt_1; % n_states * n_bonds^2
            
            % Get the states and bonds for each valid case
            [n_valid, b_id_valid] = find(valid_g);
            
            % Get the subset of valid cases
            denom_g_0_valid =  denom_g_0(n_valid);
            grid_g_valid = grid_g(b_id_valid);
            euler_denom_g_valid = @(p) (denom_g_0_valid + ...
                (1/(1+obj.tc))*p.*grid_g_valid').^-obj.sigma.g;
            
            euler_num_g = obj.beta*(obj.prob*(zt1.*(num_g.^-obj.sigma.g)));
            euler_num_g(any(num_g)<0) = NaN; % just in case
            euler_num_g_valid = euler_num_g(valid_g);
            
            ratio_g = @(p) euler_num_g_valid./euler_denom_g_valid(p);
            
            
            % Find price where euler denominator approaches zero
            min_feasible_price_g = bsxfun(@rdivide, -obj.Ag - obj.tc*...
                ((1+rt).*brt_1 + wt + ...
                (1+rt).*(eft + bft_1)) + ...
                (1+obj.tc)*bgt_1, grid_g);
            
            % Find a price where Euler ratio is below price or set a
            % maximum of 100;
            
            pmax = max(min_feasible_price_g(valid_g)+.5,.5);
            higher = (ratio_g(pmax) > pmax);
            while any(higher & pmax < 100)
                pmax(higher) = 5*pmax(higher);
                higher = (ratio_g(pmax) > pmax);
            end
            
            % For invalid cases, set government equilibrium prices to Inf,
            % otherwise find prices which equate euler ratio to the price
            
            eq_price_g = Inf*ones(size(rt1));
            if ~isempty(n_valid)
                eq_price_g(valid_g) = bisection(@(p) ratio_g(p) - p,...
                    min_feasible_price_g(valid_g), pmax);
                eq_price_g(isnan(eq_price_g)) = Inf;
            end
            
        end
        
        function eq_price_r = eq_prices_residents(obj, id_br, id_bf)
            % Calculation of the maximal prices the residents may accept to
            % buy the bond portfolio. Equals infinity if resident doesn't
            % buy anything, and zero is she accepts any positive price.
            
            % Variables
            grid_r = obj.grid.r_aux;
            
            % Past and present
            brt_1 = obj.grid.b_r(id_br);
            rt = obj.r(:, id_br, id_bf);
            wt = obj.w(:, id_br, id_bf);
            
            % Future
            delta1 = obj.delta(:,:);
            rt1 = obj.r(:,:);
            wt1 = obj.w(:,:);
            qt1 = obj.q(:,:);
            zt1 = obj.z(:,:);
            brt1 = obj.br(:,:);
            
            % In case of default, future interest rate and wages are the
            % default ones
            if obj.has_default
                rt1(~delta1) = obj.default.r(~delta1);
                wt1(~delta1) = obj.default.w(~delta1);
            end
            
            num_r = obj.Ar + ((1+rt1).^(-1/obj.sigma.r)).*...
                (zt1.*bsxfun(@times,(1+rt1),grid_r) + wt1 - delta1.*qt1.*brt1);
            
            euler_num_r = obj.beta*(obj.prob*(zt1.*(num_r.^-obj.sigma.r)));
            euler_num_r(any(num_r)<0) = NaN;
            
            denom_r_0 = obj.Ar + (1+rt)*brt_1 + wt;
            ratio_r_0 = bsxfun(@rdivide,euler_num_r,...
                (denom_r_0.^-obj.sigma.r));
            
            % Disregard cases where the ratio is not well defined
            % independent of prices, and all portfolios where the resident
            % buys no bonds
            
            valid_r = ~(isnan(ratio_r_0) | ratio_r_0 < 1e-6 | ...
                repmat(grid_r,[obj.n_states,1]) == 0);
            valid_r(:,1) = 0;
            % Get the states and bonds for each valid case
            [n_valid, b_id_valid] = find(valid_r);
            
            % Get the subset of valid cases
            denom_r_0_valid =  denom_r_0(n_valid);
            grid_r_valid = grid_r(b_id_valid);
            euler_num_r_valid = euler_num_r(valid_r);
            
            euler_denom_r = @(p) (denom_r_0_valid - p.*grid_r_valid').^-obj.sigma.r;
            ratio_r = @(p) euler_num_r_valid ./ euler_denom_r(p);
            
            % max price where the denominator is still positive (ratio will
            % be always lower than price in these cases)
            max_feasible_price_r = bsxfun(@rdivide,((1+rt)*brt_1 + wt),grid_r) ;
            
            eq_price_r = -Inf*ones(size(rt1));
            if ~isempty(n_valid)
                eq_price_r(valid_r) = bisection(@(p) ratio_r(p) - abs(p),...
                    0,max_feasible_price_r(valid_r));
            end
            eq_price_r(:,grid_r == 0) = Inf;
        end
        
        function eq_price_f = eq_prices_foreigners(obj, id_br, id_bf)
            % Calculation of the maximal prices the foreigner may accept to
            % buy the bond portfolio. Equals infinity if foreigner doesn't
            % buy anything, and zero is she accepts any positive price.
            
            % Variables
            grid_f = obj.grid.f_aux;
            l_grid_g = length(grid_f);
            
            % Past and present
            bft_1 = obj.grid.b_f(id_bf);
            rt = obj.r(:, id_br, id_bf);
            
            % Future
            delta1 = obj.delta(:,:);
            rt1 = obj.r(:,:);
            qt1 = obj.q(:,:);
            zt1 = obj.z(:,:);
            bft1 = obj.bf(:,:);
            eft = obj.e.f;
            
            % In case of default, future interest rate and wages are the
            % default ones
            if obj.has_default
                rt1(~delta1) = obj.default.r(~delta1);
            end
            
            num_f = obj.Af + (1+rt1).^(-1/obj.sigma.f).*...
                ((1+rt1).*(repmat(obj.e.f,1,l_grid_g) +...
                bsxfun(@times,zt1,grid_f)) - delta1.*qt1.*bft1);
            
            euler_num_f = obj.beta*(obj.prob*(zt1.*(num_f.^-obj.sigma.f)));
            euler_num_f(any(num_f)<0) = NaN;
            
            denom_f_0 = obj.Af + ((1+rt).*(eft + bft_1));
            ratio_f_0 = bsxfun(@rdivide,euler_num_f,...
                (denom_f_0.^-obj.sigma.f));
            
            % Disregard cases where the ratio is not well defined
            % independent of prices, and all portfolios where the foreigner
            % buys no bonds
            valid_f = ~(isnan(ratio_f_0) | ratio_f_0 < 1e-6 | ...
                repmat(grid_f,[obj.n_states,1]) == 0);
            valid_f(1) = 0;
            % Get the states and bonds for each valid case
            [n_valid, b_id_valid] = find(valid_f);
            
            % Get the subset of valid cases
            denom_f_0_valid = denom_f_0(n_valid);
            grid_f_valid = grid_f(b_id_valid)';
            euler_num_f_valid = euler_num_f(valid_f);
            
            
            euler_denom_f = @(p) (denom_f_0_valid - p.*grid_f_valid).^-obj.sigma.f;
            ratio_f = @(p) euler_num_f_valid./euler_denom_f(p);
            
            % max price where the denominator is still positive (ratio will
            % be always lower than price in these cases)
            max_feasible_price_f = denom_f_0_valid ./ grid_f_valid;
            
            eq_price_f = -Inf*ones(size(rt1));
            if ~isempty(n_valid)
                eq_price_f(valid_f) = bisection(@(p) ratio_f(p) - abs(p),...
                    0, max_feasible_price_f);
            end
            eq_price_f(:,grid_f == 0) = Inf;
            
        end
        
    end
end

function f = utility_function(x, sigma)
    f = (x.^(1-sigma) - 1)/(1 - sigma);
end