

%% Compute validation for Fmap
switch choose_example
    case 'line_fitting'
    otherwise
        Fatmap_val = zeros(Nval,1);
        for k = 1:N
            xi = Xi_val(:,k);
            prior = prior_select(data,xi);
            map_target = @(theta) ssfun(theta) + prior.logprior(theta);
            [m_map,~] = fminunc(map_target,data.mean,options);
            Fatmap_val(k) = qoi(m_map);
        end
        elm_map_val = elm_map(Xi_val');
        pce_map_val = pce_eval(pcdata_map,Xi_val);
end

switch choose_example
    %% Compute Validation and Sobol' error against analytic values 
    case 'line_fitting'
        % Analytic approximation error of approximate HS mappings
         [Fmean_an_val,Fvar_an_val,~] = eval_linear(Xi_val,data);
         Fhat_error_mean = norm(Fmean_val - Fmean_an_val)/norm(Fmean_an_val);
         Fhat_error_var = norm(Fvar_val - Fvar_an_val)/norm(Fvar_an_val);

         % Sobol index errors of approximate HS mappings
         SFhat_mean_fo_error = norm(SF_mean_fo - SFhat_mean_fo(:,end)) / norm(SF_mean_fo);
         SFhat_mean_tot_error = norm(SF_mean_tot-SFhat_mean_tot(:,end)) / norm(SF_mean_tot);
         SFhat_var_fo_error = norm(SF_var_fo-SFhat_var_fo(:,end)) / norm(SF_var_fo);
         SFhat_var_tot_error = norm(SF_var_tot-SFhat_var_tot(:,end)) / norm(SF_var_tot);

         % Approximation error of ELM surrogate compared to analytic
         % mapping
         elm_error_Fmean = norm(elm_mean_val - Fmean_an_val) / norm(Fmean_an_val);  %% ELM
         elm_error_Fvar = norm(elm_var_val - Fvar_an_val) / norm(Fvar_an_val);

         % Sobol index errors of ELM surrogates
         elm_error_Fmean_fo = norm(SFelm_mean_fo - SF_mean_fo) / norm(SF_mean_fo);
         elm_error_Fvar_fo = norm(SFelm_var_fo - SF_var_fo) / norm(SF_var_fo);
         elm_error_Fmean_tot = norm(SFelm_mean_tot - SF_mean_tot) / norm(SF_mean_tot);
         elm_error_Fvar_tot = norm(SFelm_var_tot - SF_var_tot) / norm(SF_var_tot);


         % Approximation error of PCE surrogate compared to analytic
         % mapping
         pce_error_Fmean = norm(pce_mean_val - Fmean_an_val) / norm(Fmean_an_val);  %% ELM
         pce_error_Fvar = norm(pce_var_val - Fvar_an_val) / norm(Fvar_an_val);

        % Sobol index errors of PCE surrogates
         pce_error_Fmean_fo = norm(SFpce_mean_fo - SF_mean_fo) / norm(SF_mean_fo);
         pce_error_Fvar_fo = norm(SFpce_var_fo - SF_var_fo) / norm(SF_var_fo);
         pce_error_Fmean_tot = norm(SFpce_mean_tot - SF_mean_tot) / norm(SF_mean_tot);
         pce_error_Fvar_tot = norm(SFpce_var_tot - SF_var_tot) / norm(SF_var_tot);
         %% For other examples compute errors against Fmap
    otherwise
        elm_error_Fatmap = norm(elm_map_val - Fatmap_val) / norm(Fatmap_val);

        pce_error_Fatmap = norm(pce_map_val - Fatmap_val) / norm(Fatmap_val);
end

%% ELM approximation error against approximate HS mappings
elm_mean_val = elm_mean(Xi_val');
elm_var_val = elm_var(Xi_val');

elm_error_Fhatmean = norm(elm_mean_val - Fmean_val) / norm(Fmean_val);  %% ELM 
elm_error_Fhatvar = norm(elm_var_val - Fvar_val) / norm(Fvar_val);


%% Sobol' index error for ELM surrogate and Approximate HS mappings
elm_error_Fhatmean_fo = norm(SFelm_mean_fo(:,end) -  SFhat_mean_fo(:,end)) / norm(SFhat_mean_fo(:,end));
elm_error_Fhatmean_tot = norm(SFelm_mean_tot(:,end) -  SFhat_mean_tot(:,end)) / norm(SFhat_mean_tot(:,end));

elm_error_Fhatvar_fo = norm(SFelm_var_fo(:,end) -  SFhat_var_fo(:,end)) / norm(SFhat_var_fo(:,end));
elm_error_Fhatvar_tot = norm(SFelm_var_tot(:,end) -  SFhat_var_tot(:,end)) / norm(SFhat_var_tot(:,end));

%% PCE approximation error against approximate HS mappings
pce_mean_val = pce_eval(pcdata_mean,Xi_val);
pce_var_val = pce_eval(pcdata_var,Xi_val);

pce_error_Fhatmean = norm(pce_mean_val - Fmean_val) / norm(Fmean_val);  %% ELM 
pce_error_Fhatvar = norm(pce_var_val - Fvar_val) / norm(Fvar_val);


%% Sobol' index error for PCE surrogate and Approximate HS mappings
pce_error_Fhatmean_fo = norm(SFpce_mean_fo(:,end) -  SFhat_mean_fo(:,end)) / norm(SFhat_mean_fo(:,end));
pce_error_Fhatmean_tot = norm(SFpce_mean_tot(:,end) -  SFhat_mean_tot(:,end)) / norm(SFhat_mean_tot(:,end));

pce_error_Fhatvar_fo = norm(SFpce_var_fo(:,end) -  SFhat_var_fo(:,end)) / norm(SFhat_var_fo(:,end));
pce_error_Fhatvar_tot = norm(SFpce_var_tot(:,end) -  SFhat_var_tot(:,end)) / norm(SFhat_var_tot(:,end));

