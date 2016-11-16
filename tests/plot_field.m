function plot_field(Economy,field)
    br = reshape(Economy.grid.r_aux,Economy.n_bonds,Economy.n_bonds);
    bf = reshape(Economy.grid.f_aux,Economy.n_bonds,Economy.n_bonds);
    for i=1:Economy.n_states
        figure
        surf(br,bf,squeeze(Economy.(field)(i,:,:)))
        title([field,' state ' num2str(i) ', e.f = ' num2str(Economy.e.f(i))])
        xlabel('r')
        ylabel('f')
    end
        