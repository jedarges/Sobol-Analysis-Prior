function [Mean,Var,AtMean] = eval_linear(Xi_hyp,data)

G_noise = data.cov_noise; 

N = size(Xi_hyp,2); npar = data.npar;
Mean = zeros(N,1); Var = zeros(N,1); AtMean = zeros(N,1);

mu_nom = data.mean;
cov_nom = data.cov; 
pm = data.pm;
H = data.H; 
y = data.y;

HG_noise = H' * (G_noise \ H);
Hy = H' * (G_noise \ y);

for k = 1:N
    x = Xi_hyp(:,k);
    M_shift = reshape(x(1:npar),[npar,1]);
    M = mu_nom; M = 2 * pm *  M .* M_shift + (1-pm) * M;

    V_list = reshape(x(npar+1:end),[npar,1]); C_shift = diag(V_list(1:npar));
    C = cov_nom; C = 2 * pm * C .* C_shift + (1-pm) * C;  

    Gamma_post = (HG_noise + C \ eye(npar)) \ eye(npar);
    m_post = Gamma_post * (Hy + C\M);

    atmu = m_post' * m_post;
    mu = trace(Gamma_post) + atmu;
    atmu = m_post' * m_post;
    E2 = 2*trace(Gamma_post * Gamma_post) + trace(Gamma_post)^2 + (2*trace(Gamma_post) + m_post' * m_post) * (m_post' * m_post) + 4 * m_post' * Gamma_post * m_post;
    V = E2 - mu^2;

    Mean(k) = mu;
    Var(k) = V;
    AtMean(k) = atmu;
end
