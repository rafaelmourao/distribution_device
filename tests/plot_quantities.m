function plot_quantities(Economy)
    for i=1:Economy.n_states
        br = reshape(Economy.grid.r_aux,Economy.n_bonds,Economy.n_bonds);
        bf = reshape(Economy.grid.f_aux,Economy.n_bonds,Economy.n_bonds);
        figure
        subplot(2,2,[1,2])
        surf(br,bf,squeeze(Economy.bg(i,:,:)))
        title(['Quantities, state ' num2str(i) ', e.f = ' num2str(Economy.e.f(i))])
        xlabel('r')
        ylabel('f')
        subplot(2,2,3)
        surf(br,bf,squeeze(Economy.br(i,:,:)))
        title('Resident')
        xlabel('r')
        ylabel('f')
        subplot(2,2,4)
        surf(br,bf,squeeze(Economy.bf(i,:,:)))
        title('Foreigner')
        xlabel('r')
        ylabel('f')
    end
        