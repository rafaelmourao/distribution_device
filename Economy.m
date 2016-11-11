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
        grid % structure that defines the grid of bonds
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
            
            obj.extended_grid.b_r = repmat(obj.grid.b_r,obj.n_states,1,obj.n_bonds);
            obj.extended_grid.b_f = repmat(reshape(obj.grid.b_f,1,1,obj.n_bonds),...
                obj.n_states,obj.n_bonds);
            obj.extended_grid.e_f = repmat(obj.e.f,1,obj.n_bonds,obj.n_bonds);
            
            % Value functions
            
            obj.Vd = zeros(obj.n_states,1);                                 %Value function in case of DEFAULT
            obj.Vnd = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);          %Value function in case of NO DEFAULT
            obj.Vo = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);           %Value function for the OPTION of DEFAULT
            
            %% POLICY FUNCTIONS
            
            %Consumers
            
            obj.cr = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);          %CONSUMPTION policy funtion for RESIDENTS
            obj.cf = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);          %CONSUMPTION policy funtion for FOREIGNERS
            
            obj.kr = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);          %CAPITAL policy funtion for RESIDENTS
            obj.kf = repmat(obj.e.f,1,obj.n_bonds,obj.n_bonds);           %CAPITAL policy funtion for FOREIGNERS
            
            obj.br = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);          %BONDS policy funtion for RESIDENTS
            obj.bf = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);          %BONDS policy funtion for FOREIGNERS
            
            %Government
            obj.bg = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);     %     %BONDS policy funtion for the GOVERNMENT
            obj.g = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);            %PUBLIC EXPENDITURE policy funtion for the GOVERNMENT
            obj.z = (rand(obj.n_states,obj.n_bonds,obj.n_bonds) >.1);       %DEFAULT policy funtion for the GOVERNMENT
            
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
            
            
            parfor(i=1:n_b_states,nworkers)   %RESIDENTS bonds from previous period
                
                [n, id_br, id_bf] = ind2sub(n_b_states_size,i)
                
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
            
            [status_int, p_int, br_int, bf_int, bg_int] = Interior_Solution(obj, n, id_br, id_bf);
            
            if (~status_int)
                
                [status_cor, p_cor, br_cor, bf_cor, bg_cor] = Corner_Solution(obj, n, id_br, id_bf);
                
            end
            
            if ( status_int || ~status_cor )
                p = p_int;
                br_s = br_int;
                bf_s = bf_int;
                bg_s = bg_int;
            else
                p = p_cor;
                br_s = br_cor;
                bf_s = bf_cor;
                bg_s = bg_cor;
            end
            
        end
        
        function [status, p, br, bf, B] = Interior_Solution(obj, n, id_br, id_bf)
            % Parameters
            
            epsilon = .1;
            
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
            
            %YOU DON'T NEED TO CONSIDER 'zt' IN THE DENOMINADOR, SINCE
            %THE BOND MARKET MUST BE OPEN FOR THIS CALCULATION DO BE DONE
            %YOU DO NEED 'zt1'
            
            denom_g = @(p) (obj.tc/(1+obj.tc))*...
                (((1+rt)*brt_1 + wt - p*grid_r') + ...
                ((1+rt)*(eft + bft_1) - p*grid_f')) + ...
                (p*grid_g' - bgt_1);
            
            num_g = (obj.tc/(1+obj.tc))*...
                ((zt1.*bsxfun(@times,(1+rt1),grid_r) + ...
                wt1 - zt1.*qt1.*brt1) + ...
                ((1+rt1).*(repmat(obj.e.f,1,l_grid_g) + ...
                bsxfun(@times,zt1,grid_f)) - zt1.*qt1.*bft1)) + ...
                zt1.*(qt1.*bgt1 - repmat(grid_g,obj.n_states,1));
            
            ratio_g = @(p) Euler_ratio(obj, n, zt1, num_g, denom_g(p),obj.sigma.g);
            euler_g = @(p) abs(ratio_g(p) - p);
            
            denom_r = @(p) ((1+rt)*brt_1 + wt - p*grid_r');
            num_r = ((1+rt1).^(-1/obj.sigma.r)).*...
                (zt1.*bsxfun(@times,(1+rt1),grid_r) + wt1 - zt1.*qt1.*brt1);
            ratio_r = @(p) Euler_ratio(obj, n, zt1, num_r, denom_r(p),obj.sigma.r);
            euler_r = @(p) abs(p - ratio_r(p));
            
            
            denom_f = @(p) ((1+rt)*(eft + bft_1) - p*grid_f');
            num_f = (1+rt1).^(-1/obj.sigma.f).*...
                ((1+rt1).*(repmat(obj.e.f,1,l_grid_g) +...
                bsxfun(@times,zt1,grid_f)) - zt1.*qt1.*bft1);
            ratio_f = @(p) Euler_ratio(obj, n, zt1,num_f,denom_f(p),obj.sigma.f);
            euler_f = @(p) abs(p - ratio_f(p));
            
            f = @(p) min(euler_g(p) + euler_r(p) + euler_f(p));
            
            [p, sum_euler] = fminbnd(f,0,obj.p_max);
            
            [~, b_star] = f(p);
            
            B = grid_g(b_star);
            br = grid_r(b_star);
            bf = grid_f(b_star);
            
            status = (sum_euler < epsilon);
            
        end
        
        function [status, p, br, bf, B] = ...
                Corner_Solution(obj, n, id_br, id_bf)
            
            [status, p, br, bf, B, pr, pf, pg] = ...
                Corner_Solution_1(obj, n, id_br, id_bf);
            
            if (status)
                return
            end
            
            [status, p, br, bf, B] = Corner_Solution_2(obj, n, id_br, id_bf);
            
            if (status)
                return
            end
            
            [status, p, br, bf, B] = Corner_Solution_3(obj, n, id_br, id_bf);
            
            if (status)
                return
            end
            
            [status, p, br, bf, B] = Corner_Solution_others(obj, ...
                pr, pf, pg);
        end
        
        function [status, p, br, bf, B, pr, pf, pg] = ...
                Corner_Solution_1(obj, n, id_br, id_bf)
            
            B = 0;
            br = 0;
            bf = 0;
            
            % Future
            rt1 = obj.r(:,1,1);
            wt1 = obj.w(:,1,1);
            qt1 = obj.q(:,1,1);
            zt1 = obj.z(:,1,1);
            brt1 = obj.br(:,1,1);
            bft1 = obj.bf(:,1,1);
            bgt1 = obj.bg(:,1,1);
            brt_1 = obj.br(n, id_br, id_bf);
            bft_1 = obj.bf(n, id_br, id_bf);
            bgt_1 = obj.bg(n, id_br, id_bf);
            rt = obj.r(n, id_br, id_bf);
            wt = obj.w(n, id_br, id_bf);
            eft = obj.e.f(n);
            
            denom_r = ((1+rt)*brt_1 + wt);
            num_r = ((1+rt1).^(-1/obj.sigma.r)).*...
                (wt1 - zt1.*qt1.*brt1);
            
            pr = Euler_ratio(obj,n,zt1,num_r,denom_r,obj.sigma.r);
            
            denom_f = (1+rt)*(bft_1 + eft);
            num_f = ((1+rt1).^(-1/obj.sigma.f)).*...
                (bsxfun(@times,(1+rt1),obj.e.f) - zt1.*qt1.*bft1);
            
            pf = Euler_ratio(obj,n,zt1,num_f,denom_f,obj.sigma.f);
            
            denom_g = (obj.tc/(1+obj.tc))*(((1+rt)*brt_1 + wt) +...
                ((1+rt)*(bft_1 + eft))) - bgt_1;
            
            num_g = (obj.tc/(1+obj.tc))*((wt1 - zt1.*qt1.*brt1) +...
                (bsxfun(@times,(1+rt1),obj.e.f) - zt1.*qt1.*bft1)) +...
                zt1.*qt1.*bgt1;
            
            pg = Euler_ratio(obj, n, zt1, num_g, denom_g,obj.sigma.g);
            
            p = max(pr,pf);
            
            status = (p < pg);
            
        end
        
        function [status, p, br, bf, B] = Corner_Solution_2(obj, n, id_br, id_bf)
            
            epsilon = .1;
            br = 0;
            
            grid_f = obj.grid.b_f;
            l_grid_f = length(grid_f);
            
            % Past and current
            bgt_1 = obj.bg(n, id_br, id_bf);
            bft_1 = obj.bf(n, id_br, id_bf);
            brt_1 = obj.br(n, id_br, id_bf);
            rt = obj.r(n, id_br, id_bf);
            wt = obj.w(n, id_br, id_bf);
            eft = obj.e.f(n);
            
            % Future
            bft1 = squeeze(obj.bf(:,1,:));
            %crt = obj.cr(n, id_br, id_bf);
            rt1 = squeeze(obj.r(:,1,:));
            wt1 = squeeze(obj.w(:,1,:));
            qt1 = squeeze(obj.q(:,1,:));
            zt1 = squeeze(obj.z(:,1,:));
            brt1 = squeeze(obj.br(:,1,:));
            bgt1 = squeeze(obj.bg(:,1,:));
            %crt1 = squeeze(obj.cr(:,1,:));
            
            % Euler equations
            
            denom_g = @(p) (obj.tc/(1+obj.tc))*((1+rt)*brt_1 + wt) + ... %Resident consumption when br=0
                (obj.tc/(1+obj.tc))*((1+rt)*(eft + bft_1) - p*grid_f') -...
                bgt_1 + p*grid_f';
            num_g = (obj.tc/(1+obj.tc))*(wt1 - zt1.*qt1.*brt1) +... %Resident consumption when br=0
                (obj.tc/(1+obj.tc))*((1+rt1).*(repmat(obj.e.f,1,l_grid_f) +...
                repmat(grid_f,obj.n_states,1)) - zt1.*qt1.*bft1) - ...
                bsxfun(@times,zt1,grid_f) + qt1.*zt1.*bgt1;
            ratio_g = @(p) Euler_ratio(obj,n,zt1,num_g,denom_g(p),obj.sigma.g);
            euler_g = @(p) abs(ratio_g(p) - p);
            
            denom_f = @(p) (1+rt)*(eft + bft_1) - p*grid_f';
            num_f = ((1+rt1).^(-1/obj.sigma.f)).*...
                ((1+rt1).*(repmat(obj.e.f,1,l_grid_f) +...
                bsxfun(@times,zt1,grid_f) - zt1.*qt1.*bft1));
            ratio_f = @(p) Euler_ratio(obj,n,zt1,num_f,denom_f(p),obj.sigma.f);
            euler_f = @(p) abs(p - ratio_f(p));
            
            % objective function
            f = @(p) min(euler_g(p) + euler_f(p));
            
            [p, sum_euler] = fminbnd(f,0,obj.p_max);
            [~, b_star] = f(p);
            B = grid_f(b_star);
            bf = grid_f(b_star);
            
            rt1_c2 = obj.r(:,1,b_star);
            wt1_c2 = obj.w(:,1,b_star);
            qt1_c2 = obj.q(:,1,b_star);
            zt1_c2 = obj.z(:,1,b_star);
            brt1_c2 = obj.br(:,1,b_star);
            
            denom_r_c2 = (1+rt)*brt_1 + wt;
            num_r_c2 = ((1+rt1_c2).^(-1/obj.sigma.r)).*...
                (wt1_c2 - zt1_c2.*qt1_c2.*brt1_c2);
            
            pr_0_c2 = Euler_ratio(obj, n, zt1_c2, num_r_c2, denom_r_c2,obj.sigma.r);
            
            status = (sum_euler < epsilon & pr_0_c2 <= p);
            
        end
        
        function [status, p, br, bf, B] = Corner_Solution_3(obj, n, id_br, id_bf)
            %This is the case where the FOREIGN investor chooses zero bonds:'bf = 0'
            
            bf = 0;
            epsilon = .1;
            
            %Grid
            eft = obj.e.f(n);
            grid_r = obj.grid.b_r;
            
            %State
            bft_1 = obj.bf(n, id_br, id_bf);
            brt_1 = obj.br(n, id_br, id_bf);
            bgt_1 = obj.bg(n, id_br, id_bf);
            rt = obj.r(n, id_br, id_bf);
            wt = obj.w(n, id_br, id_bf);
            
            
            %Government
            rt1 = obj.r(:,:,1);
            wt1 = obj.w(:,:,1);
            qt1 = obj.q(:,:,1);
            zt1 = obj.z(:,:,1);
            brt1 = obj.br(:,:,1);
            bft1 = obj.bf(:,:,1);
            bgt1 = obj.bg(:,:,1);
            
            %% ALGORITHM
            
            denom_g = @(p) (obj.tc/(1+obj.tc))*...
                ((1+rt)*brt_1 + wt - p*grid_r') +...
                (obj.tc/(1+obj.tc))*(1+rt)*(bft_1 + eft) -... Foreigner consumption when bf=0
                bgt_1 + p*grid_r';
            num_g = (obj.tc/(1+obj.tc))*...
                (bsxfun(@times,(1+rt1),grid_r) + wt1 - zt1.*qt1.*brt1) + ...
                (obj.tc/(1+obj.tc))*(bsxfun(@times,(1+rt1),obj.e.f) - zt1.*qt1.*bft1) - ... Foreigner consumption when bf=0
                bsxfun(@times,zt1,grid_r) + zt1.*qt1.*bgt1;
            ratio_g = @(p) Euler_ratio(obj, n, zt1, num_g, denom_g(p),obj.sigma.g);
            euler_g = @(p) abs(ratio_g(p) - p);
            
            
            denom_r = @(p) (1+rt)*brt_1 + wt - p*grid_r';
            num_r = (1+rt1).^(-1/obj.sigma.r).*...
                (zt1.*bsxfun(@times,(1+rt1),grid_r) + wt1 - zt1.*qt1.*brt1);
            ratio_r = @(p) Euler_ratio(obj,n,zt1,num_r,denom_r(p),obj.sigma.r);
            euler_r = @(p) abs(p - ratio_r(p));
            
            f = @(p) min(euler_g(p) + euler_r(p));
            
            [p, sum_euler] = fminbnd(f,0,obj.p_max);
            
            [~, b_star] = f(p);
            
            B = grid_r(b_star);
            br = grid_r(b_star);
            
            rt1_c3 = obj.r(:,b_star,1);
            qt1_c3 = obj.q(:,b_star,1);
            zt1_c3 = obj.z(:,b_star,1);
            bft1_c3 = obj.bf(:,b_star,1);
            
            denom_f_c3 = (1+rt)*(obj.e.f(n) + bft_1);
            num_f_c3 = (1+rt1_c3).^(-1/obj.sigma.f).*...
                (bsxfun(@times,(1+rt1_c3),obj.e.f) - zt1_c3.*qt1_c3.*bft1_c3);
            
            pf_0_c3 = Euler_ratio(obj, n, zt1_c3, num_f_c3, denom_f_c3, obj.sigma.f);
            
            status = (sum_euler < epsilon && pf_0_c3 < p);
            
            
        end
        
        function [status, p, br, bf, B] = ...
                Corner_Solution_others(~, pr_c1, pf_c1, pg_c1 )
            
            %CASE 4: B = 0 and br = 0
            B = 0;
            br = 0;
            bf = 0;
            pr_0_c4 = pr_c1;
            pf_0_c4 = pf_c1;
            pg_0_c4 = pg_c1;
            p_c4 = pf_0_c4;
            status = (p_c4 < pg_0_c4) && (p_c4 > pr_0_c4);
            
            if status
                return
            end
            
            
            %CASE 5: B = 0 and bf = 0
            B = 0;
            br = 0;
            bf = 0;
            pr_0_c5 = pr_c1;
            pf_0_c5 = pf_c1;
            pg_0_c5 = pg_c1;
            p = pr_0_c5;
            status = (p < pg_0_c5) && (p > pf_0_c5);
            
            if status
                return
            end
            
            %CASE 6: br = 0 and bf = 0
            B = 0;
            br = 0;
            bf = 0;
            pr_0_c6 = pr_c1;
            pf_0_c6 = pf_c1;
            pg_0_c6 = pg_c1;
            p = pg_0_c6;
            status = (p > pr_0_c6) && (p > pf_0_c6);
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