function prior = prior_select(data,x)
spec = data.prior_type;
switch spec
    case 'multi_no_corr'
        n = data.npar;% pm = spec.pm;
        M_shift = reshape(x(1:n),[n,1]);
        pm = data.pm;
        V_list = reshape(x(n+1:end),[n,1]); C_shift = diag(V_list(1:n));
        C = data.cov; Cov = 2 * pm * C .* C_shift + (1-pm) * C; 
      
        M = data.mean; mu = 2 * pm *  M .* M_shift + (1-pm) * M;
        %M = data.mu; mu =  M + 2 * M_shift - 1;
        D = det(Cov);
        prior.priorfun = @(th) (2*pi)^(-n/2) * sqrt(D) * exp(-1/2 * (th - mu)'* (Cov \ (th-mu)));
        prior.logprior = @(th)  n *log(2*pi) - 2 * log(sqrt(D)) + (th - mu)'* (Cov \ (th-mu));
end