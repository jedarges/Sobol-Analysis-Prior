function prior = is_prior_select(data,type)

switch type
    case 'gaussian'
         n = data.npar;
         mu = data.mean;
         Cov = data.cov_is;
         D = det(Cov);

         prior.priorfun = @(th) (2*pi)^(-n/2) * sqrt(D) * exp(-1/2 * (th - mu)'* (Cov \ (th-mu))); 
         prior.logprior = @(th)  n * log(2*pi) - 2 * log(sqrt(D)) + (th - mu)'* (Cov \ (th-mu));
end