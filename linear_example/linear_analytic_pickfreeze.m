switch choose_example
    case 'line_fitting'
        [Fmean_array_A,Fvar_array_A,~] = eval_linear(A,data);
        [Fmean_array_B,Fvar_array_B,~] = eval_linear(B,data);

        Fmean_array_Aj = zeros(Ns,nhyp);
        Fvar_array_Aj = zeros(Ns,nhyp);
        for j = 1:nhyp
            Aj = A; % Hybrid matrix created for each input
            Aj(j,:)=B(j,:);
            [Fmean_Aj,Fvar_Aj,~] = eval_linear(Aj,data);
            Fmean_array_Aj(:,j) = Fmean_Aj;
            Fvar_array_Aj(:,j) = Fvar_Aj;
        end

        SF_mean_fo = zeros(nhyp,1); SF_mean_tot = zeros(nhyp,1);
        SF_var_fo = zeros(nhyp,1); SF_var_tot = zeros(nhyp,1);

        f_mean_A = Fmean_array_A; f_mean_B = Fmean_array_B;
        f_var_A = Fvar_array_A; f_var_B = Fvar_array_B;

        var_mean = (var(f_mean_A) + var(f_mean_B))/2;
        var_var = (var(f_var_A) + var(f_var_B))/2;
        for j = 1:nhyp
            f_mean_Aj = Fmean_array_Aj(:,j);
            f_var_Aj = Fvar_array_Aj(:,j);

            % Mean
            SF_mean_fo(j) = sum(f_mean_B .* (f_mean_Aj - f_mean_A)) / (Ns * var_mean);
            SF_mean_tot(j) = (1/(2 * Ns * var_mean)) * sum((f_mean_A - f_mean_Aj).^2);
            % Var

            SF_var_fo(j) = sum(f_var_B .* (f_var_Aj - f_var_A)) / (Ns * var_var);
            SF_var_tot(j) = (1/(2 * Ns * var_var)) * sum((f_var_A - f_var_Aj).^2);
        end
       
    otherwise
end