function [W,bias,beta] = elm_train_model(Xd,Yd,model_data)
%Yd = reshape(Yd,'',1); 
s_sz = size(Yd,1);
%Xd = reshape(Xd,s_sz,''); 
ndim = size(Xd,2);

%% User choice parameters
p = model_data.p;
lambda = model_data.lambda;
N = model_data.nneurons;

%% Model parameters
phi = @(x) exp(x);


%% Create ELM
% Weight matrix
W = randn(ndim,N);
Z = rand(ndim,N) > p;
W = W .* Z;
% bias vector
bias = randn(1,N);

H = phi(Xd * W + bias);

H = [H; lambda * eye(N)];
Yd = [Yd; zeros(N,1)];
[U,S,V] = svd(H); Sinv = diag(1./diag(S' * S));
z = Sinv * S' * U' * Yd;
beta = V' \ z;

end