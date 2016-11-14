classdef Economy
    properties
        % economy parameters
        beta
        sigma
        phi
        lambda
        tc
        alpha
        rho
        e
        prob
        n_states
        grid % structure that defines the grid of bonds, optional
        n_bonds
        default
        p_max
        
        % iteraction variables
        Vd
        Vnd
        Vo
        Wnd
        cr
        cf
        kr
        kf
        br
        bf
        bg
        z
        r
        w
        q
        g
    end
    
    properties
        extended_grid
    end
    
    methods
        
        function obj = Economy(param)
            % Setting parameters
            obj.p_max = param.p_max;
            obj.beta = param.beta;
            obj.sigma = param.sigma;
            obj.phi = param.phi;
            obj.lambda = param.lambda;
            obj.tc = param.tc;
            obj.alpha = param.alpha;
            obj.rho = param.rho;
            obj.e = param.e;
            obj.prob = param.prob;
            obj.n_states = param.n_states;
            obj.n_bonds = param.n_bonds;
            
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
            obj.br = reshape(kron(obj.grid.r_aux,[1,1,1]),...
                [obj.n_states,obj.n_bonds,obj.n_bonds]);
%             obj.bf = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);          %BONDS policy funtion for FOREIGNERS
            obj.bf = reshape(kron(obj.grid.f_aux,[1,1,1]),...
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
            
            obj.q = ones(obj.n_states,obj.n_bonds,obj.n_bonds);           %Price of Public Bond
            
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
                obj.default.g(n) = obj.tc*(obj.default.cr(n) + obj.default.cf(n));
                obj.default.W(n) =  Utility_Function(obj.default.cr(n),obj.sigma.r) +...
                    obj.lambda*Utility_Function(obj.default.cf(n),obj.sigma.f) ...
                    +  Utility_Function(obj.default.g(n),obj.sigma.g);
                
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
                next_Vo(i) = obj.Vo(n,loc_br_s,loc_bf_s);
                
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
            
            obj.g = obj.tc*(obj.cr + obj.cf) + obj.q.*obj.bg - ...
                (obj.extended_grid.b_r + obj.extended_grid.b_f);
            obj.g(obj.g<0) = 0;
            
            obj.Wnd = Utility_Function(obj.cr,obj.sigma.r) + ...
                obj.lambda*Utility_Function(obj.cf,obj.sigma.f) + ...
                Utility_Function(obj.g,obj.sigma.g);
            
            for i = 1:obj.n_bonds
                obj.Vnd(:,:,i) = obj.Wnd(:,:,i) + obj.beta*(obj.prob * next_Vo(:,:,i)); %%%
            end
            
            obj.Vd = obj.default.W + obj.beta*obj.prob*...
                (obj.phi * obj.Vo(:,1,1) + (1-obj.phi) * obj.Vd);
            
            % in case of no default
            obj.z = ones(obj.n_states,obj.n_bonds,obj.n_bonds);
            obj.Vo = obj.Vnd;
            
            % check where there is default
            def = bsxfun(@lt,obj.Vnd,obj.Vd);
            [def_states, ~] = find(def); % retrieving the states of all default occ.
            obj.z(def) = 0;
            obj.Vo(def) = obj.Vd(def_states);
            obj.r(def) = obj.default.r(def_states);
            obj.w(def) = obj.default.w(def_states);
            obj.cr(def) = obj.default.cr(def_states);
            obj.cf(def) = obj.default.cf(def_states);
            obj.g(def) = obj.default.g(def_states);
            obj.q(def) = 0;
            obj.br(def) = 0;
            obj.bf(def) = 0;
            obj.bg(def) = 0;
            obj.kr(def) = 0;
            obj.kf(def) = obj.e.f(def_states);
        end
        
        function [p, br_s, bf_s, bg_s] = Solution(obj, n, id_br, id_bf)
            
            % Default (corner) solution
            
            p = 0;
            br_s = 0;
            bf_s = 0;
            bg_s = 0;
            
            % VARIABLES NEEDED
            
            %Grid
            
            grid_r = obj.grid.r_aux;
            grid_f = obj.grid.f_aux;
            grid_g = obj.grid.b_g;
            l_grid_g = length(grid_g);
            
            % Past and present
            
            brt_1 = obj.br(n, id_br, id_bf);
            bft_1 = obj.bf(n, id_br, id_bf);
            bgt_1 = obj.bg(n, id_br, id_bf);
            rt = obj.r(n, id_br, id_bf);
            wt = obj.w(n, id_br, id_bf);
            eft = obj.e.f(n);
            
            % Future
            
            bgt1 = obj.bg(:,:);
            rt1 = obj.r(:,:);
            wt1 = obj.w(:,:);
            qt1 = obj.q(:,:);
            zt1 = obj.z(:,:);
            brt1 = obj.br(:,:);
            bft1 = obj.bf(:,:);
            
            % ALGORITHM
            
            % government
            
            num_g = (obj.tc/(1+obj.tc))*...
                ((zt1.*bsxfun(@times,(1+rt1),grid_r) + ...
                wt1 - zt1.*qt1.*brt1) + ...
                ((1+rt1).*(repmat(obj.e.f,1,l_grid_g) + ...
                bsxfun(@times,zt1,grid_f)) - zt1.*qt1.*bft1)) + ...
                zt1.*(qt1.*bgt1 - repmat(grid_g,obj.n_states,1));
            
            valid_g = all(num_g > 0);
            valid_g(1) = 0;
            
            % min price where the denominator is still positive
            min_feasible_price_g = ( - obj.tc*...
                ((1+rt)*brt_1 + wt + ...
                (1+rt)*(eft + bft_1)) + ...
                (1+obj.tc)*bgt_1 ) ./ grid_g';
            
            
            denom_g_0 = (obj.tc/(1+obj.tc))*...
                (((1+rt)*brt_1 + wt) + ...
                ((1+rt)*(eft + bft_1))) - ...
                bgt_1;
                        
            grid_g_valid = grid_g(valid_g);
            zt1_g_valid = zt1(:,valid_g);
            num_g_valid = num_g(:,valid_g);
            
            denom_g = @(p) denom_g_0 + ...
                (1/(1+obj.tc))*p.*grid_g_valid';
            
            ratio_g = @(p) Euler_ratio(obj, n, zt1_g_valid,...
                num_g_valid, denom_g(p),obj.sigma.g);
            
            pmax = obj.p_max;
            while any(ratio_g(pmax) < pmax & pmax < 1000)
                pmax = 10*pmax;
            end
                        
            eq_price_g = Inf*ones(1,l_grid_g);
            if any(valid_g)
                eq_price_g(valid_g) = max(bisection(@(p) ratio_g(p) - abs(p),...
                    min_feasible_price_g(valid_g), repmat(pmax,sum(valid_g),1)),0);
            end
            
            % residents
            
            num_r = ((1+rt1).^(-1/obj.sigma.r)).*...
                (zt1.*bsxfun(@times,(1+rt1),grid_r) + wt1 - zt1.*qt1.*brt1);
            
            % max price where the denominator is still positive
            max_feasible_price_r = ( (1+rt)*brt_1 + wt ) ./ grid_r' ;
            
            denom_r_0 = ((1+rt)*brt_1 + wt);
            ratio_r_0 = Euler_ratio(obj, n, zt1,...
                num_r, denom_r_0,obj.sigma.r);
            
            valid_r = ~(isnan(ratio_r_0) | ratio_r_0 < 1e-6 | grid_r' == 0);
            valid_r(1) = 0;
            
            grid_r_valid = grid_r(valid_r);
            num_r_valid = num_r(:,valid_r);
            zt1_r_valid = zt1(:,valid_r);
            
            denom_r = @(p) (denom_r_0 - p.*grid_r_valid');
            ratio_r = @(p) Euler_ratio(obj, n, zt1_r_valid,...
                num_r_valid, denom_r(p),obj.sigma.r);
            
            
            eq_price_r = -Inf*ones(1,l_grid_g);
            if any(valid_r)
                eq_price_r(valid_r) = bisection(@(p) ratio_r(p) - abs(p),...
                    0,max_feasible_price_r(valid_r));
            end
            eq_price_r(grid_r == 0) = 1e6;
            
            % foreigners
            
            num_f = (1+rt1).^(-1/obj.sigma.f).*...
                ((1+rt1).*(repmat(obj.e.f,1,l_grid_g) +...
                bsxfun(@times,zt1,grid_f)) - zt1.*qt1.*bft1);
            
            % max price where the denominator is still positive
            max_feasible_price_f = ( (1+rt)*(eft + bft_1) ) ./ grid_f';
            
            denom_f_0 = ((1+rt)*(eft + bft_1));
            ratio_f_0 = Euler_ratio(obj, n, zt1,...
                num_f, denom_f_0, obj.sigma.f);
            
            valid_f = ~(isnan(ratio_f_0) | ratio_f_0 < 1e-6 | grid_f' == 0);
            valid_f(1) = 0;
            
            grid_f_valid = grid_f(valid_f);
            num_f_valid = num_f(:,valid_f);
            zt1_f_valid = zt1(:,valid_f);
            
            denom_f = @(p) (denom_f_0 - p.*grid_f_valid');
            ratio_f = @(p) Euler_ratio(obj, n, zt1_f_valid,...
                num_f_valid, denom_f(p), obj.sigma.f);
            
            
            eq_price_f = -Inf*ones(1,l_grid_g);
            if any(valid_f)
                eq_price_f(valid_f) = bisection(@(p) ratio_f(p) - abs(p),...
                    0, max_feasible_price_f(valid_f));
            end
            eq_price_f(grid_f == 0) = 1e6;
            
            % finding market equilibrium
            
            prices = min(eq_price_f,eq_price_r);
            prices(eq_price_g > prices) = 0;
            
            [~, sorted_bonds] = sort(grid_g,2,'descend');
            for i = sorted_bonds
                if ( prices(i) > 0 )
                    p = prices(i);
                    br_s = grid_r(i);
                    bf_s = grid_f(i);
                    bg_s = grid_g(i);
                    break
                end
            end
        end
        
        function ratio = Euler_ratio(obj,n, zt1, x, y, sig)
            
            % Important remark: the 'x' refers to the NUMERATOR while the 'y' refers to
            % the DENOMINATOR. The former is inouted as a 'n_states x length_grid'
            % matriz and the later as a '1 x length_grid' vector.
            
            % VARIABLES NEEDED
            
            %Parameters
            probt = obj.prob(n,:);
            
            % PROGRAM
            
            %ERROR FUNTION
            %This function identifies where there is a problem with consumption (it
            %can never be negative), for both vectors: 'x' and 'y'. They will be part
            %of the numerator and denominator, respectively, of the euler equations.
            %The first part identifies the coordinates where the denominator has
            %negative values. The second part finds if there is any entry, for every
            %state of nature in the future, where consumption might be negative; the
            %whole column will be shut down.
            
            error_c = (y < 0) | any(x < 0)';
            
            numerator = obj.beta*(probt*(zt1.*(x.^-sig)))';
            
            denominator = y.^(-sig);
            
            %RATIO FUNCTION
            %This function calculates the MRS. But if there is any negative consumtion
            %on any of the parts used to calculate this ratio, will show up a 'NaN'
            %(not a number) on that coordinate, for the '0/0' that would be calculated.
            %The Matlab won't consider this coordinates with 'NaN' when trying to find
            %the policy functions.
            
            ratio = numerator./denominator;
            ratio(error_c) = NaN;
            
        end
        
    end
    
end