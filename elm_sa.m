function [S1,ST,elmmodel] = elm_sa(Xd_train,Yd_train,Xd_ver,Yd_ver)
rng('default')
if size(Yd_train,1) < size(Yd_train,2)
    Yd_train = Yd_train';
end
if size(Xd_train,1) < size(Xd_train,2)
    Xd_train = Xd_train';
end
if size(Yd_ver,1) < size(Yd_ver,2)
    Yd_ver = Yd_ver';
end
if size(Xd_ver,1) < size(Xd_ver,2)
    Xd_ver = Xd_ver';
end

out_train = Yd_train;
in_train = Xd_train;

out_ver = Yd_ver;
in_ver = Xd_ver;

p_list = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95 0.99];
%p_list = [0.95];
trials = 1;

%% ELM Settings
model_data.nneurons = round(2 * length(out_train)/3); % Number of hidden layer neurons, i.e., ELM basis functions
model_data.lambda = 1e-5; 
p_sz = length(p_list);

%% SW training
rel = zeros(p_sz,1);
for j = 1:p_sz
    model_data.p = p_list(j);
    for l = 1:trials
        [W,bias,beta] = elm_train_model(in_train,out_train,model_data);
        elm_valid = exp(in_ver * W + bias) * beta;
        rel(j) = rel(j) + norm(out_ver- elm_valid) / norm(out_ver);
    end
end

p_opt = p_list(rel == min(rel));
model_data.p = p_opt;
%disp(p_opt)
[W,bias,beta] = elm_train_model(in_train,out_train,model_data);

elmmodel = @(theta) exp(theta * W + bias) * beta;
%elm_valid = exp(in_valid * W + bias) * beta;
%error = norm(Yd_valid - elm_valid) / norm(Yd_valid);

[S1,ST,~] = elm_sobol_inds(W,beta,bias);