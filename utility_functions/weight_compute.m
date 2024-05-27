function weights = weight_compute(Xi_hyp,chain,data)

iter = max(size(chain));
N = max(size(Xi_hyp));
weights = zeros(iter,N);
npar = data.npar;
pm = data.pm;

mu_is = data.mean;
cov_is = data.cov_is;

mu_nom = data.mean;
cov_nom = data.cov;

% IS weight portion
W_is = sum((chain - mu_is').^2 * diag(1 ./ diag(cov_is)),2);
for k = 1:N
    x = Xi_hyp(:,k);

    M_shift = reshape(x(1:npar),[npar,1]);
    M = mu_nom; M = 2 * pm *  M .* M_shift + (1-pm) * M;

    V_list = reshape(x(npar+1:end),[npar,1]); C_shift = diag(V_list(1:npar));
    C = cov_nom; C = 2 * pm * C .* C_shift + (1-pm) * C;  

    W_k =  sum((chain - M').^2 * diag(1./diag(C)),2);
    weights(:,k) = exp(-1/2 * (W_k - W_is));
end