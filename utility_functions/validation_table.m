switch choose_example
    case 'line_fitting'
        HSMean = [Fhat_error_mean; elm_error_Fmean; pce_error_Fmean; elm_error_Fhatmean; pce_error_Fhatmean;
            SFhat_mean_fo_error; elm_error_Fmean_fo; pce_error_Fmean_fo; elm_error_Fhatmean_fo; pce_error_Fhatmean_fo;
            SFhat_mean_tot_error; elm_error_Fmean_tot; pce_error_Fmean_tot; elm_error_Fhatmean_tot; pce_error_Fhatmean_tot];

        HSVariance = [Fhat_error_var; elm_error_Fvar; pce_error_Fvar; elm_error_Fhatvar; pce_error_Fhatvar;
            SFhat_var_fo_error; elm_error_Fvar_fo; pce_error_Fvar_fo; elm_error_Fhatvar_fo; pce_error_Fhatvar_fo;
            SFhat_var_tot_error; elm_error_Fvar_tot; pce_error_Fvar_tot; elm_error_Fhatvar_tot; pce_error_Fhatvar_tot];

        Error_type = {'Error(F,F^M)','Error(F,F^M_{ELM})', 'Error(F,F^M_{PCE})','Error(F^M,F^M_{ELM})', 'Error(F^M,F^M_{PCE})','S_{err}(F,F^M)','S_{err}(F,F^M_{ELM})', 'S_{err}(F,F^M_{PCE})','S_{err}(F^M,F^M_{ELM})','S_{err}(F^M,F^M_{PCE})','S^{tot}_{err}(F,F^M)','S^{tot}_{err}(F,F^M_{ELM})', 'S^{tot}_{err}(F,F^M_{PCE})','S^{tot}_{err}(F^M,F^M_{ELM})','S^{tot}_{err}(F^M,F^M_{PCE})'};

        T = table(HSMean,HSVariance,'RowNames',Error_type);
        disp(T)


    otherwise
        HSMean = [NaN; NaN; elm_error_Fhatmean; pce_error_Fhatmean;
         elm_error_Fhatmean_fo; pce_error_Fhatmean_fo;
         elm_error_Fhatmean_tot; pce_error_Fhatmean_tot];

        HSVariance = [NaN; NaN; elm_error_Fhatvar; pce_error_Fhatvar;
            elm_error_Fhatvar_fo; pce_error_Fhatvar_fo;
            elm_error_Fhatvar_tot; pce_error_Fhatvar_tot];

        HSAtMAP = [elm_error_Fatmap; pce_error_Fatmap; NaN; NaN;
            NaN; NaN; NaN; NaN];


        Error_type = {'Error(F,F_{ELM})','Error(F,F_{PCE})','Error(F^M,F^M_{ELM})','Error(F^M,F^M_{PCE})','S_{err}(F^M,F^M_{ELM})','S_{err}(F^M,F^M_{PCE})','S^{tot}_{err}(F^M,F^M_{ELM})','S^{tot}_{err}(F^M,F^M_{PCE})'};
        T = table(HSMean,HSVariance,HSAtMAP,'RowNames',Error_type);
        disp(T)
end