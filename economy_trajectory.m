function econ = economy_trajectory(Economy,id_br,id_bf,state,...
        num_periods, seed)
    
    if (exist('seed','var') && (seed > 0))
        rng(seed)
    end
    
    n_tot_bonds_size = size(Economy.Vo);
    states = markov_chain(Economy.prob,state,num_periods);
    idx_s = 0 * states;
    idx_s(1) = sub2ind(n_tot_bonds_size,state,id_br,id_bf);
    econ.delta = ones(1,num_periods);
    
    for i = 1:num_periods-1
        if (i == 1) || ( Economy.delta(idx_s(i)) && ...
                ( econ.delta(i-1) || rand(1) < Economy.phi ) )
            id_br_s = find(Economy.grid.b_r == Economy.br(idx_s(i)),1);
            id_bf_s = find(Economy.grid.b_f == Economy.bf(idx_s(i)),1);
        else
            id_br_s = 1;
            id_bf_s = 1;
            econ.delta(i) = 0;
        end
        idx_s(i+1) = sub2ind(n_tot_bonds_size,states(i+1),id_br_s,id_bf_s);
    end
    
    econ.idx = idx_s;
    econ.states = states;
    
    % finding the periods with default
    def_cases = find(~econ.delta);
    def_cases_states = states(def_cases);
    
    % case without default
    fields = {'q','br','bf','bg','r','w','cr','cf','kr','kf','g'};
    for i = 1:length(fields)
        field_name = fields{i};
        econ.(field_name) = Economy.(field_name)(idx_s);
    end
    
    %cases with default
    fields = {'r','w','cr','cf','g'};
    for i = 1:length(fields)
        econ.(field_name)(def_cases) = ...
            Economy.default.(field_name)(def_cases_states);
    end
    econ.br(def_cases) = 0;
    econ.bg(def_cases) = 0;
    econ.bf(def_cases) = 0;
end

function chain = markov_chain(transition, first_state, length_chain)
    chain = [first_state, zeros(1,length_chain-1)];
    r = rand(1,length_chain);
    cum_probs = cumsum(transition,2);
    if any(cum_probs(:,end) ~= 1)
        error('Invalid transition matrix, lines must sum to 1!')
    end
    for i = 2:length(chain)
        chain(i) = find(cum_probs(chain(i-1),:) > r(i), 1);
    end
end



