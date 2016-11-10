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
            obj.grid = param.grid;
            
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
            obj.bg = 1e-10*ones(obj.n_states,obj.n_bonds,obj.n_bonds);     %BONDS policy funtion for the GOVERNMENT
            obj.g = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);           %PUBLIC EXPENDITURE policy funtion for the GOVERNMENT
            obj.z = ones(obj.n_states,obj.n_bonds,obj.n_bonds);            %DEFAULT policy funtion for the GOVERNMENT
            
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
            
            nbonds = obj.n_bonds;
            nstates = obj.n_states;
            
            p = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);
            br_s = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);
            bf_s = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);
            bg_s = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);
            next_Vo = zeros(obj.n_states,obj.n_bonds,obj.n_bonds);
            
            parfor(id_br=1:nbonds,nworkers)   %RESIDENTS bonds from previous period

                for id_bf = 1:nbonds                    %FOREIGNERS bonds from previous period
                    
                    for n = 1:nstates                          %State of Nature
                        
                        [p(n,id_br,id_bf), br_s(n,id_br,id_bf),...
                            bf_s(n,id_br,id_bf), bg_s(n,id_br,id_bf)] = ...
                            Solution(obj, n, id_br, id_bf);
                        
                        loc_br_s = ( obj.grid.b_r == br_s(n,id_br,id_bf) );
                        loc_bf_s = ( obj.grid.b_f == bf_s(n,id_br,id_bf) );
                        next_Vo(n,id_br,id_bf) = obj.Vo(n,loc_br_s,loc_bf_s);
                        
                    end
                    
                end
                
            end
            
            obj.q = p;                              %Equilibrium price for the Bonds' market
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
            
            obj.cf = (1/(1+obj.tc))*((1+obj.r).*obj.kf...
                - obj.q.*obj.bf);
            
            obj.g = obj.tc*(obj.cr + obj.cf) + obj.q.*obj.bg - ...
                (obj.extended_grid.b_r + obj.extended_grid.b_f);
            
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
            extended_Vd = repmat(obj.Vd,1,obj.n_bonds,obj.n_bonds);
            def = (obj.Vnd < extended_Vd);
            obj.z(def) = 0;
            obj.Vo(def) = extended_Vd(def);
        end
        
        function [p, br_s, bf_s, bg_s] = Solution(obj, n, id_br, id_bf)
            
            [status_int, p_int, br_int, bf_int, bg_int] = Interior_Solution(obj, n, id_br, id_bf);
            
            if ~status_int
                
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
            p_max = 10;
            
            % VARIABLES NEEDED
            
            %Grid
            grid_r = obj.grid.r_aux;
            grid_f = obj.grid.f_aux;
            grid_g = obj.grid.b_g;
            l_grid_g = length(grid_g);
            
            brt_1 = obj.br(n, id_br, id_bf);
            bft_1 = obj.bf(n, id_br, id_bf);
            bgt_1 = obj.bg(n, id_br, id_bf);
            rt = obj.r(n, id_br, id_bf);
            wt = obj.w(n, id_br, id_bf);
            eft = obj.e.f(n);
            %Government
            bgt1 = obj.bg(:,:);
            %Investors
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
            
            denom_g = @(p) obj.tc*(((1+rt)*brt_1 + wt - p*grid_r') + ...
                ((1+rt)*(eft + bft_1) - p*grid_f')) + ...
                (p*grid_g' - bgt_1);
            num_g = (obj.tc*((bsxfun(@times,(1+rt1).*zt1,grid_r) + wt1 - zt1.*qt1.*brt1) + ...
                ((1+rt1).*(repmat(obj.e.f,1,l_grid_g) + bsxfun(@times,zt1,grid_f)) - zt1.*qt1.*bft1)) + ...
                zt1.*(qt1.*bgt1 - repmat(grid_g,obj.n_states,1)));
            ratio_g = @(p) Euler_ratio(obj, n, zt1, num_g, denom_g(p),obj.sigma.g);
            euler_g = @(p) abs(ratio_g(p) - p);
            
            denom_r = @(p) ((1+rt)*brt_1 + wt - p*grid_r');
            num_r = ((1+rt1).^(-1/obj.sigma.r)).*...
                (bsxfun(@times,(1+rt1).*zt1,grid_r) + wt1 - zt1.*qt1.*brt1);
            ratio_r = @(p) Euler_ratio(obj, n, zt1, num_r, denom_r(p),obj.sigma.r);
            euler_r = @(p) abs(p - ratio_r(p));
            
            
            denom_f = @(p) ((1+rt)*(eft + bft_1) - p*grid_f');
            num_f = ((1+rt1).^(-1/obj.sigma.f)).*...
                ((1+rt1).*(repmat(obj.e.f,1,l_grid_g) + bsxfun(@times,zt1,grid_f)) - zt1.*qt1.*bft1);
            ratio_f = @(p) Euler_ratio(obj, n, zt1,num_f,denom_f(p),obj.sigma.f);
            euler_f = @(p) abs(p - ratio_f(p));
            
            f = @(p) min(euler_g(p) + euler_r(p) + euler_f(p));
            
            [p, sum_euler] = fminbnd(f,0,p_max);
            
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
            
            if status
                return
            end
            
            [status, p, br, bf, B] = Corner_Solution_2(obj, n, id_br, id_bf);
            
            if status
                return
            end
            
            [status, p, br, bf, B] = Corner_Solution_3(obj, n, id_br, id_bf);
            
            if status
                return
            end
            
            [status, p, br, bf, B] = Corner_Solution_others(obj, ...
                pr, pf, pg);
        end
        
        
        function [status, p, br, bf, B, pr, pf, pg] = ...
                Corner_Solution_1(obj, n, id_br, id_bf)
            
            B = 0;
            br = 0;
            bf = 1e-10;
            rt1 = obj.r(:,1,1);
            wt1 = obj.w(:,1,1);
            qt1 = obj.q(:,1,1);
            zt1 = obj.z(:,1,1);
            brt1 = obj.br(:,1,1);
            bft1 = obj.bf(:,1,1);
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
            
            bgt1 = obj.bg(:,1,1);
            
            denom_g = obj.tc*(((1+rt)*brt_1 + wt) + ((1+rt)*(bft_1 + eft))) - bgt_1;
            num_g = (obj.tc*((wt1 - zt1.*qt1.*brt1) + (bsxfun(@times,(1+rt1),obj.e.f) - zt1.*qt1.*bft1)) +...
                zt1.*qt1.*bgt1);
            
            pg = Euler_ratio(obj, n, zt1, num_g, denom_g,obj.sigma.g);
            
            p = max(pr,pf);
            
            status = (p < pg);
            
        end
        
        function [status, p, br, bf, B] = Corner_Solution_2(obj, n, id_br, id_bf)
            
            epsilon = .1;
            br = 0;
            p_max = 10;
            
            %variables for optimization
            
            grid_f = obj.grid.b_f;
            l_grid_f = length(grid_f);
            
            % State
            
            bgt_1 = obj.bg(n, id_br, id_bf);
            bft_1 = obj.bf(n, id_br, id_bf);
            brt_1 = obj.br(n, id_br, id_bf);
            rt = obj.r(n, id_br, id_bf);
            wt = obj.w(n, id_br, id_bf);
            eft = obj.e.f(n);
            
            % Investors
            
            bft1 = squeeze(obj.bf(:,1,:));
            % Government
            
            crt = obj.cr(n, id_br, id_bf);
            rt1 = squeeze(obj.r(:,1,:));
            qt1 = squeeze(obj.q(:,1,:));
            zt1 = squeeze(obj.z(:,1,:));
            bgt1 = squeeze(obj.bg(:,1,:));
            crt1 = squeeze(obj.cr(:,1,:));
            
            % Euler equations
            
            denom_g = @(p) obj.tc*(crt + ((1+rt)*(eft + bft_1) - p*grid_f')) - bgt_1 + p*grid_f';
            num_g = (obj.tc*(crt1 + ((1+rt1).*(repmat(obj.e.f,1,l_grid_f) +...
                repmat(grid_f,obj.n_states,1)) - zt1.*qt1.*bft1)) - ...
                repmat(grid_f,obj.n_states,1) + qt1.*zt1.*bgt1);
            ratio_g = @(p) Euler_ratio(obj,n,zt1,num_g,denom_g(p),obj.sigma.g);
            euler_g = @(p) abs(ratio_g(p) - p);
            
            denom_f = @(p) ((1+rt).^(-1/obj.sigma.f)).*((1+rt)*(eft + bft_1) - p*grid_f');
            num_f = ((1+rt1).^(-1/obj.sigma.f)).*...
                ((1+rt1).*(repmat(obj.e.f,1,l_grid_f) + repmat(grid_f,obj.n_states,1)) - zt1.*qt1.*bft1);
            ratio_f = @(p) Euler_ratio(obj,n,zt1,num_f,denom_f(p),obj.sigma.f);
            euler_f = @(p) abs(p - ratio_f(p));
            
            % objective function
            f = @(p) min(euler_g(p) + euler_f(p));
            
            [p, sum_euler] = fminbnd(f,0,p_max);
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
            
            if sum_euler > epsilon || pr_0_c2 > p
                
                status = 0;
                p = 0;
                B = 1e-10;
                bf = 1e-10;
                
            end
            
        end
        
        function [status, p, br, bf, B] = Corner_Solution_3(obj, n, id_br, id_bf)
            %This is the case where the FOREIGN investor chooses zero bonds:'bf = 0'
            
            bf = 0;
            status = 1;
            p_max = 10;
            epsilon = .1;
            
            %Grid
            grid_r = obj.grid.b_r;
            
            %State
            bft_1 = obj.br(n, id_br, id_bf);
            brt_1 = obj.br(n, id_br, id_bf);
            bgt_1 = obj.bg(n, id_br, id_bf);
            rt = obj.r(n, id_br, id_bf);
            wt = obj.w(n, id_br, id_bf);
            
            brt1 = squeeze(obj.br(:,1,:));
            
            %Government
            crt = obj.cr(n, id_br, id_bf);
            cft = obj.cf(n, id_br, id_bf);
            rt1 = obj.r(:,:,1);
            wt1 = obj.w(:,:,1);
            qt1 = obj.q(:,:,1);
            zt1 = obj.z(:,:,1);
            bgt1 = obj.bg(:,:,1);
            crt1 = obj.cr(:,:,1);
            cft1 = obj.cf(:,:,1);
            
            %% ALGORITHM
            
            denom_g = @(p) obj.tc*(((1+rt)*brt_1 + wt - p*grid_r') + cft) - bgt_1 + p*grid_r';
            num_g = ((obj.tc*((bsxfun(@times,(1+rt1),grid_r) + wt1 - zt1.*qt1.*brt1) + cft1) - ...
                repmat(grid_r,obj.n_states,1) + zt1.*qt1.*bgt1));
            ratio_g = @(p) Euler_ratio(obj, n, zt1, num_g, denom_g(p),obj.sigma.g);
            euler_g = @(p) abs(ratio_g(p) - p);
            
            
            denom_r = @(p) ((1+rt)*brt_1 + wt - p*grid_r');
            num_r = ((1+rt1).^(-1/obj.sigma.r)).*...
                (bsxfun(@times,(1+rt1),grid_r) + wt1 - zt1.*qt1.*brt1);
            ratio_r = @(p) Euler_ratio(obj,n,zt1,num_r,denom_r(p),obj.sigma.r);
            euler_r = @(p) abs(p - ratio_r(p));
            
            f = @(p) min(euler_g(p) + euler_r(p));
            
            [p, sum_euler] = fminbnd(f,0,p_max);
            
            [~, b_star] = f(p);
            
            B = grid_r(b_star);
            br = grid_r(b_star);
            
            rt1_c3 = obj.r(:,b_star,1);
            qt1_c3 = obj.q(:,b_star,1);
            zt1_c3 = obj.z(:,b_star,1);
            bft1_c3 = obj.bf(:,b_star,1);
            
            denom_f_c3 = (1+rt)*(obj.e.f(n) + bft_1);
            num_f_c3 = ((1+rt1_c3).^(-1/obj.sigma.f)).*...
                (bsxfun(@times,(1+rt1_c3),obj.e.f) - zt1_c3.*qt1_c3.*bft1_c3);
            
            pf_0_c3 = Euler_ratio(obj, n, zt1_c3, num_f_c3, denom_f_c3, obj.sigma.f);
            
            if sum_euler > epsilon || pf_0_c3 > p
                status = 0;
                p = 0;
                B = 1e-10;
                br = 0;
            end
            
            
        end
        
        function [status, p, br, bf, B] = ...
                Corner_Solution_others(obj, pr_c1, pf_c1, pg_c1 )
            
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
                                  
            %CALCULATING THE NUMERATOR
            %This functions computes the numerator itself, but for the coordinates
            %where there was negative consumption for either numerator or denominator,
            %there will be '0' instead.
            
            numerator = obj.beta*(probt*((zt1.*x).^(-sig)))';
            
            %FIX DENOMINATOR FUNTION
            %This function corrects the denominator, placing 'Inf' whenever there is
            %a negative consumption on either the numerator or denominator.
                        
            %CALCULATING THE DENOMINATOR
            %This functions computes the denomitor itself, but for the coordinates
            %where there was negative consumption for either numerator or denominator,
            %there will be '0' instead.
            
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