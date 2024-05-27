function data = line_fitting_data

ndim = 1; 
npar = ndim + 1;
coeff = (-1).^(1:(ndim+1))' .* (1:(ndim+1))'; % true coefficients

t_data = [0; .5; 1.5; 2]; % Data meaurements at times t
d_sz = length(t_data);
t_data = [ones(d_sz,1) t_data]; 

%% Instance of sampled noise added to data
noise =   [-0.5587;
    0.1784;
    -0.1969;
    0.5864];
y_data = t_data * coeff + noise;

%% Construct covariance for proposal distribution
H = t_data;
J = H' * H;

mopt = J \ (t_data' * y_data);

n = length(y_data); p = length(mopt);
rss = norm(H*mopt - y_data)^2;
mse = rss / (n-p);
G = mse * inv(J);

%% Nominal prior hyperparameters
mu = [1; 1];

Cov = [1 0;
    0 1];

%% Sum of squares
ssfun = @(params,data) sum((t_data * params - y_data).^2);

%% Quantity of interest
qoi_non = @(params) norm(params)^2;
%qoi_lin = @(params) dot(params,[1 5]);


%% Data
data.mean = mu;% prior mean
data.cov = Cov; % prior covariance
data.ssfun = ssfun; % ssfun
data.prop = G;% proposal
data.npar = npar;% parameter dimension
data.nhyp = 2* npar;% hyperparameter dimension
data.qoi = qoi_non;% QoI
data.cov_is = 1.5^2 * eye(npar);
data.cov_noise = eye(4);
data.H = H;
data.y = y_data;
data.sd = 2.3^2 / npar;


