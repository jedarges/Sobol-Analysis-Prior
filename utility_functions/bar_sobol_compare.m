switch choose_example
    case 'line_fitting'
        % F_mean
        figure;
        subplot(2,1,1)
        bar([SF_mean_fo(:,end) SFhat_mean_fo(:,end) SF_mean_fo(:,end) SFelm_mean_fo])
        legend('F_{mean}','F^M_{mean}','SW-ELM','PCE','Location','NorthWest')
        ylabel('First-order Sobol'' Index')
        set(gca,'XTickLabel',{'\mu_{b}','\mu_{m}','\sigma_b^2','\sigma_m^2'})
        set(gca,'fontsize',18)
        ylim([0,0.5])

        subplot(2,1,2)
        bar([SF_mean_tot(:,end) SFhat_mean_tot(:,end) SF_mean_tot(:,end) SFelm_mean_tot])
        legend('F_{mean}','F^M_{mean}','SW-ELM','PCE','Location','NorthWest')
        ylabel('Total Sobol'' Index')
        set(gca,'XTickLabel',{'\mu_{b}','\mu_{m}','\sigma_b^2','\sigma_m^2'})
        set(gca,'fontsize',18)
        ylim([0,0.5])

        %F_var
        figure;
        subplot(2,1,1)
        bar([SF_var_fo(:,end) SFhat_var_fo(:,end) SFelm_var_fo SFpce_var_fo])
       legend('F_{var}','F^M_{var}','SW-ELM','PCE','Location','NorthWest')
        ylabel('First-order Sobol'' Index')
        set(gca,'XTickLabel',{'\mu_{b}','\mu_{m}','\sigma_b^2','\sigma_m^2'})
        set(gca,'fontsize',18)
        ylim([0,0.6])

        subplot(2,1,2)
        bar([SF_var_tot(:,end) SFhat_var_tot(:,end) SFelm_var_tot SFpce_var_tot])
        legend('F_{var}','F^M_{var}','SW-ELM','PCE','Location','NorthWest')
        ylabel('Total Sobol'' Index')
        set(gca,'XTickLabel',{'\mu_{b}','\mu_{m}','\sigma_b^2','\sigma_m^2'})
        set(gca,'fontsize',18)
        ylim([0,0.6])
    case 'seir'
        % F_mean
        figure;
        subplot(2,1,1)
        bar([SFhat_mean_fo(:,end) SFelm_mean_fo SFpce_mean_fo])
        legend('F^M_{mean}','SW-ELM','PCE','Location','NorthEast')
        ylabel('First-order Sobol'' Index')
       set(gca,'XTickLabel', {'m_{log\mu}','m_{log\beta}','m_{log\sigma}','m_{log\gamma}','s^2_{log\mu}','s^2_{log\beta}','s^2_{log\sigma}','s^2_{log\gamma}'})
       set(gca,'fontsize',18)
        ylim([0,0.8])

        subplot(2,1,2)
        bar([SFhat_mean_tot(:,end) SFelm_mean_tot SFpce_mean_fo])
        legend('F^M_{mean}','SW-ELM','PCE','Location','NorthEast')
        ylabel('Total Sobol'' Index')
       set(gca,'XTickLabel', {'m_{log\mu}','m_{log\beta}','m_{log\sigma}','m_{log\gamma}','s^2_{log\mu}','s^2_{log\beta}','s^2_{log\sigma}','s^2_{log\gamma}'})
       set(gca,'fontsize',18)
        ylim([0,0.8])

        % F_var
        figure;
        subplot(2,1,1)
        bar([SFhat_var_fo(:,end) SFelm_var_fo SFpce_var_fo])
       legend('F^M_{var}','SW-ELM','PCE','Location','NorthEast')
        ylabel('First-order Sobol'' Index')
        set(gca,'XTickLabel', {'m_{log\mu}','m_{log\beta}','m_{log\sigma}','m_{log\gamma}','s^2_{log\mu}','s^2_{log\beta}','s^2_{log\sigma}','s^2_{log\gamma}'})
        set(gca,'fontsize',18)
        %ylim([0,0.6])

        subplot(2,1,2)
        bar([SFhat_var_tot(:,end) SFelm_var_tot SFpce_var_tot])
        legend('F^M_{var}','SW-ELM','PCE','Location','NorthEast')
        ylabel('Total Sobol'' Index')
        set(gca,'XTickLabel', {'m_{log\mu}','m_{log\beta}','m_{log\sigma}','m_{log\gamma}','s^2_{log\mu}','s^2_{log\beta}','s^2_{log\sigma}','s^2_{log\gamma}'})
        set(gca,'fontsize',18)
        %ylim([0,0.6])

        %F_map
        figure;
        subplot(2,1,1)
        bar([SFelm_map_fo SFpce_map_fo])
       legend('SW-ELM','PCE','Location','NorthEast')
        ylabel('First-order Sobol'' Index')
        set(gca,'XTickLabel', {'m_{log\mu}','m_{log\beta}','m_{log\sigma}','m_{log\gamma}','s^2_{log\mu}','s^2_{log\beta}','s^2_{log\sigma}','s^2_{log\gamma}'})
        set(gca,'fontsize',18)
        ylim([0,0.9])

        subplot(2,1,2)
        bar([SFelm_map_tot SFpce_map_tot])
        legend('SW-ELM','PCE','Location','NorthEast')
        ylabel('Total Sobol'' Index')
        set(gca,'XTickLabel', {'m_{log\mu}','m_{log\beta}','m_{log\sigma}','m_{log\gamma}','s^2_{log\mu}','s^2_{log\beta}','s^2_{log\sigma}','s^2_{log\gamma}'})
        set(gca,'fontsize',18)
        ylim([0,0.9])
    otherwise
        figure;
        subplot(2,1,1)
        bar([SFhat_mean_fo(:,end) SFelm_mean_fo SFpce_mean_fo])
        legend('F^M_{mean}','SW-ELM','PCE','Location','NorthEast')
        ylabel('First-order Sobol'' Index')
        set(gca,'fontsize',18)

        subplot(2,1,2)
        bar([SFhat_mean_tot(:,end) SFelm_mean_tot SFpce_mean_tot])
        legend('F^M_{mean}','SW-ELM','PCE','Location','NorthEast')
        ylabel('Total Sobol'' Index')
        set(gca,'fontsize',18)

        figure;
        subplot(2,1,1)
        bar([SFhat_var_fo(:,end) SFelm_var_fo SFpce_var_fo])
        legend('F^M_{var}','SW-ELM','PCE','Location','NorthEast')
        ylabel('First-order Sobol'' Index')
        set(gca,'fontsize',18)

        subplot(2,1,2)
        bar([SFhat_var_tot(:,end) SFelm_var_tot SFpce_var_tot])
        legend('F^M_{var}','SW-ELM','PCE','Location','NorthEast')
        ylabel('Total Sobol'' Index')
        set(gca,'fontsize',18)

        figure;
        subplot(2,1,1)
        bar([SFelm_map_fo SFpce_map_fo])
        legend('SW-ELM','PCE','Location','NorthEast')
        ylabel('Total Sobol'' Index')
        set(gca,'fontsize',18)

        subplot(2,1,2)
        bar([SFelm_map_tot SFpce_map_tot])
        legend('SW-ELM','PCE','Location','NorthEast')
        ylabel('Total Sobol'' Index')
        set(gca,'fontsize',18)
end