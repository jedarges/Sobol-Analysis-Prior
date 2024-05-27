function data = seir_data

%%
t_data = (3 * (1:15) + 30)'; % measurement times t
nt = length(t_data);


%% True parameters
mu = 5.48e-5; % birth rate
beta = 1 / 2.5; % Infection rate
sigma = 1/3; % incubation/latency rate
gamma = 1/7; % removal/Recovery rate
params = [mu beta sigma gamma];

npar = 4; % parameter dimension

%% R_0
qoi = @(th) ((exp(th(2)))/(exp(th(4))+exp(th(1)))) * (exp(th(3))/(exp(th(3))+exp(th(1))));
disp(beta / (gamma + mu) * sigma / (sigma + mu))

%% Initial conditions for SEIR compartments
N = 1000; E0 = 0; I0 = 1; R0 = 0; S0 = N - E0 - I0 - R0;
Y0 = [S0 E0 I0 R0];
meas.Y0 = Y0; meas.N = N;

% noise = [0.3714   -0.2256    1.1174   -1.0891    0.0326    0.5525    1.1006    1.5442    0.0859   -1.4916 -0.7423   -1.0616    2.3505   -0.6156    0.7481];
%% Create noise data measurements for infections at times t
[t2,Y2] = ode45(@(t,y) SEIR_rhs(t,y,params,N),[0 100],Y0);% Y01 = Y0; Y01(2) = []; [t1,Y1] = ode45(@(t,y) SIR_rhs3(t,y,params),[0 1],Y01);

s = spline(t2,Y2(:,3),t_data);
rng(1)
noise = 30; % noise variance
y_data = s + noise * randn(nt,1);


%% Sum of squares function
meas.xdata = t_data;
meas.ydata = y_data;
meas.range = [0,100];

ssfun = @(th) SEIR_ss([exp(th)],meas);


%% Construct covariance for proposal distribution
n_guess = -1 * ones(1,4);

options = optimoptions(@fminunc,'MaxFunctionEvaluations',10000,'Display','off');
[q_opt,~] = fminunc(ssfun,n_guess,options);

[t,Y] = ode45(@(t,y) SEIR_rhs(t,y,[exp(q_opt)],N),[0 100],Y0);
I = spline(t,Y(:,3),t_data);

Res = (I - y_data)';
sigma2_est = (1/(nt-npar))*Res*Res';
evals = zeros(6,nt);
evals(1,:) = y_data; evals(2,:) = I;


h = 1e-6; Sens = zeros(npar,nt);
for i=1:npar
    q = q_opt; q(i) = q(i) + h;
    [t,Ypert] = ode45(@(t,y) SEIR_rhs(t,y,[exp(q)],N),[0 100],Y0);
    Ipert = spline(t,Ypert(:,3),t_data);
    Sens(i,:) = (I - Ipert) / h;
    evals(i+2,:) = Ipert;
end

V = sigma2_est*inv(Sens*Sens');

%% Nominal prior hyperparameters
Cov = eye(npar);
mu = [-10; -1.5; -1.5; -1.5];


%% Data
data.cov_is = 4 * eye(npar);
data.mean = mu;% prior mean
data.cov = Cov; % prior covariance
data.ssfun =  @(th) ssfun(th') / (noise^2); % ssfun
data.prop = V / 1e6;% proposal
data.npar = npar;% parameter dimension
data.nhyp = 2 * npar;% hyperparameter dimension
data.qoi = qoi;% QoI
data.sd = 1 / npar;



