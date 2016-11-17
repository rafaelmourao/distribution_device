% plot figures
addpath('../')
addpath('../plots')

param2

iter = Economy(param);

epsilon = 1e-3;                                     %Tolerance level
dist = 100;                                         %Distance between previous and current price and bond functions
t = 1;                                              %Number of interations

while dist > epsilon && t <= 10000
    tic
    t = t+1;
    
    old_iter = iter;
    iter = iter.update();

    time = toc;
    dist = max(abs(iter.Vo(:) - old_iter.Vo(:)));
    fprintf('Iter: %d, distance: %.6f, time: %.2f seconds\n',t,dist,time) 
    disp('Prices:')
    disp(iter.q(1:10))
    fprintf('Default proportion: %.4f\n',1-mean(iter.z(:))) 
    
end

plot_prices(iter)
plot_quantities(iter)
plot_default(iter)