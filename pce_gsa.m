function [s1,sT,pcdata] = pce_gsa(Xd_train,Yd_train,choose_example)

if size(Yd_train,1) < size(Yd_train,2)
    Yd_train = Yd_train';
end
if size(Xd_train,1) < size(Xd_train,2)
    Xd_train = Xd_train';
end

t_set= Yd_train;
Xd = Xd_train;

nmeas = size(t_set,1);
nhyp = size(Xd,2);

switch choose_example
    case 'line_fitting'
        nord = 10; 
        tau = [1e-6 1e-4 1e-3 .01]; tsz = length(tau);
    case 'seir'
        nord = 5;
        tau = [.001 .005 .01 .05]; tsz = length(tau);
end

%S1 = zeros(nhyp,1); ST1 = zeros(nhyp,1);
%sigs = [.3 .42 .425 .425 .4211];
%%
[multiIndex, ~] = uq_initMultiIndex(nord, nhyp);
nPCTerms = size(multiIndex,1);
psiMultiDSq = zeros(nPCTerms,1);
for i = 1:nPCTerms
    mind = multiIndex(i,:);
    psiMultiDSq(i) = prod(1 ./ (2*mind + 1.0));
end
pcdata.ndim = nhyp;
pcdata.multiIndex = multiIndex;
pcdata.nPCTerms = nPCTerms;
pcdata.pc_type = 'LEGENDRE';
pcdata.psiMultiDSq = psiMultiDSq;
pcdata.nord = nord;

Psi_eval = zeros(nmeas,nPCTerms);

for i = 1:nmeas
    xi = Xd(i,:); x = 2 * xi - 1;
    for l = 1:nPCTerms
        mind = multiIndex(l,:);
        vals = zeros(nhyp,1);
        for n = 1:nhyp
            vals(n) = leg_eval(mind(n),x(n));
        end
        Psi_eval(i,l) = prod(vals);
    end
end

C_iter = 10;
c_step = round(nmeas/C_iter);
%% PCE
loo = zeros(tsz,C_iter); %err = zeros(1,sz);
for c = 1:C_iter
    H = Psi_eval; H((1+(c-1)*c_step):(c*c_step),:) = []; Hout = Psi_eval((1+(c-1)*c_step):(c*c_step),:);
    Y = t_set; Y((1+(c-1)*c_step):(c*c_step)) = []; yout = t_set((1+(c-1)*c_step):(c*c_step));
    G = H'*H; G_dim = length(G);
    opts = spgSetParms('verbosity',0);
    for i=1:tsz
        tau_i = tau(i);
        [x,] = spgl1(H,Y,0,tau_i,[],opts);
        loo(i,c) = sum((Hout * x - yout).^2);
        disp(loo(i,c));
    end
end
loo = sum(loo,2);
ind = loo == min(loo);
sigma = tau(ind); disp(sigma)
%sigma = sigs(q);
opts = spgSetParms('verbosity',0);
H = Psi_eval; Y = t_set;
[pce,~] = spgl1(H,Y,0,sigma,[],opts);

sobol1 = zeros(nhyp,1); sobolT1 = zeros(nhyp,1);
for i = 1:nhyp
    sobol1(i) = pce_sobol_fo(pce,pcdata,i);
    sobolT1(i) = pce_sobol_tot(pce,pcdata,i);
end

%S1(:,q) = sobol1; ST1(:,q) = sobolT1;
% S_mcmc(:,p,Q) = sobol1; ST_mcmc(:,p,Q) = sobolT1;
%pce_eval(q,:) = uq_sample(pcdata, pce, 1e5);

pcdata.pce = pce;
s1 = sobol1; sT = sobolT1;
