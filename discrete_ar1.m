function [S, P] = discrete_ar1(n, mu, gamma, nu, kappa)

    % grid for shock
    nu_s = nu/sqrt(1-gamma^2);
    mu_s = mu/(1-gamma);
    smin = mu_s - kappa*nu_s;
    smax = mu_s + kappa*nu_s;
    sstep = (smax - smin)/(n-1);
    S = exp(linspace(smin,smax,n))';

    % transition matrix
    normarg1 = (smin - mu - gamma*log(S))/nu + 0.5*sstep/nu;
    P(:,1) = 0.5 + 0.5*erf(normarg1/sqrt(2));
    for j = 2:(n-1)
        normarg1 = (log(S(j)) - mu - gamma*log(S))/nu + 0.5*sstep/nu;
        normarg2 = (log(S(j)) - mu - gamma*log(S))/nu - 0.5*sstep/nu;
        P(:,j) = 0.5*erf(normarg1/sqrt(2)) - 0.5*erf(normarg2/sqrt(2));
    end
    P(:,n) = 1 - sum(P(:,1:(n-1)),2);
end