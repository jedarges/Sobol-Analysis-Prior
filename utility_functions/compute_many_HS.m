function [w,q,q2] = compute_many_HS(Xi,data,chain,ql,batch_mcmc,batch)
    
%M = length(qoi_chain);
Ns = max(size(Xi));
Nc = max(size(chain));

iter = round(Ns / batch);
iter_c = round(Nc / batch_mcmc);

q = zeros(Ns,1); q2 = zeros(Ns,1);
w = zeros(Ns,1);

ql2 = ql.^2;


for k = 1:iter
    Xi_k = Xi(:,((k - 1) * batch + 1):(k * batch));

    for l = 1:iter_c
        chain_k = chain(((l-1)*batch_mcmc+1):l*batch_mcmc,:);
        q_k = ql(((l-1)*batch_mcmc+1):l*batch_mcmc);
        q_k2 = ql2(((l-1)*batch_mcmc+1):l*batch_mcmc);

        weights_A = weight_compute(Xi_k,chain_k,data);
        w((1 + batch * (k-1)) : (k * batch)) = w((1 + batch * (k-1)) : (k * batch)) + sum(weights_A,1)';
        q((1 + batch * (k-1)) : (k * batch)) = q((1 + batch * (k-1)) : (k * batch)) + weights_A' * q_k;
        q2((1 + batch * (k-1)) : (k * batch)) = q2((1 + batch * (k-1)) : (k * batch)) + weights_A' * q_k2;
    end
end

