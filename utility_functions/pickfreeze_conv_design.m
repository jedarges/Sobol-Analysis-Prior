for l = 1:step_iter
    temp_mean_A = fhat_mean_A(1:(step_pf*l));
    temp_mean_B = fhat_mean_B(1:(step_pf*l));
    temp_var_A = fhat_var_A(1:(step_pf*l));
    temp_var_B = fhat_var_B(1:(step_pf*l));

    var_mean = (var(temp_mean_A) + var(temp_mean_B))/2;
    var_var = (var(temp_var_A) + var(temp_var_A))/2;

    for j = 1:nhyp
        temp_mean_Aj =  Fhatmean_array_Aj(1:(step_pf*l),j);
        temp_var_Aj =  Fhatvar_array_Aj(1:(step_pf*l),j);

        Mpropl = (2 * l * step_pf);

       Sfo_mean_pf(j,l) = (norm(temp_mean_B)^2 + norm(temp_mean_Aj)^2) / Mpropl / var_mean - (sum(temp_mean_B + temp_mean_Aj) / Mpropl)^2 / var_mean;
       % Sfo_mean_pf(j,l) = sum(temp_mean_B .* (temp_mean_Aj - temp_mean_A)) / (l * step_pf * var_mean);
        Stot_mean_pf(j,l) = (1/(2 * (l*step_pf) * var_mean)) * norm(temp_mean_A - temp_mean_Aj)^2;


         Sfo_var_pf(j,l) = (norm(temp_var_B)^2 + norm(temp_var_Aj)^2) / (Mpropl * var_var) - (sum(temp_var_B + temp_var_Aj) / (2*l * step_pf))^2 / var_var;
      %  Sfo_var_pf(j,l) = temp_var_B' * (temp_var_Aj - temp_var_A) / (l * step_pf * var_var);
        Stot_var_pf(j,l) = (1/(2 * (l*step_pf) * var_var)) * norm(temp_var_A - temp_var_Aj)^2;
    end
end

