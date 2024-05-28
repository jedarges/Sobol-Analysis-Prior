for l = 1:PF_iter
    chain_k = chain_reduced(((l-1)*PF_step+1):l*PF_step,:);
    q_k = qoi_chain_reduced(((l-1)*PF_step+1):l*PF_step);
    
    [w_int,q_int,q2_int] = compute_many_HS(A,data,chain_k,q_k,b_mcmc,b_pf);

    w_A = w_A + w_int;
    q_A = q_A + q_int;
    q2_A = q2_A + q2_int;

    fhat_mean_A = q_A ./ w_A;
    fhat_var_A = q2_A ./ w_A - fhat_mean_A.^2;

    [w_int,q_int,q2_int] = compute_many_HS(B,data,chain_k,q_k,b_mcmc,b_pf);

    w_B = w_B + w_int;
    q_B = q_B + q_int;
    q2_B = q2_B + q2_int;

    fhat_mean_B = q_B ./ w_B;
    fhat_var_B = q2_B ./ w_B - fhat_mean_B.^2;

    Fhatmean_array_Aj = zeros(Ns,nhyp);
    Fhatvar_array_Aj = zeros(Ns,nhyp);

    for j = 1:nhyp
        Aj = A; % Hybrid matrix created for each input
        Aj(j,:)=B(j,:);

        [w_int,q_int,q2_int] = compute_many_HS(Aj,data,chain_k,q_k,b_mcmc,b_pf);

        w_Aj(:,j) = w_Aj(:,j) + w_int;
        q_Aj(:,j) = q_Aj(:,j) + q_int;
        q2_Aj(:,j) = q2_Aj(:,j) + q2_int;
        
        Fhatmean_array_Aj(:,j) = q_Aj(:,j) ./  w_Aj(:,j);
        Fhatvar_array_Aj(:,j) = q2_Aj(:,j) ./  w_Aj(:,j) - Fhatmean_array_Aj(:,j).^2;
    end

    var_mean = (var(fhat_mean_A) + var(fhat_mean_B))/2;
    var_var = (var(fhat_var_A) + var(fhat_var_B))/2;

    for j = 1:nhyp
        fhat_mean_Aj = Fhatmean_array_Aj(:,j);
        fhat_var_Aj = Fhatvar_array_Aj(:,j);

        % Mean
        SFhat_mean_fo(j,l) = sum(fhat_mean_B .* (fhat_mean_Aj - fhat_mean_A)) / (Ns * var_mean);
        SFhat_mean_tot(j,l) = (1/(2 * Ns * var_mean)) * sum((fhat_mean_A - fhat_mean_Aj).^2);
        % Var

        SFhat_var_fo(j,l) = sum(fhat_var_B .* (fhat_var_Aj - fhat_var_A)) / (Ns * var_var);
        SFhat_var_tot(j,l) = (1/(2 * Ns * var_var)) * sum((fhat_var_A - fhat_var_Aj).^2);
    end
  save("pfconv_gsa.mat","A","B","SFhat_mean_fo","SFhat_mean_tot","SFhat_var_fo","SFhat_var_tot")
end