linetype = {'--','-.','-',':','--','-.','-','--','-.','-','--','-','-.','--','-','--','-.','-','--','-.','-','--','-','-.','--','-','--','-.','-','--','-.','-','--','-','-.','--','-','--','-.','-','--','-.','-','--'};
 switch choose_example
     case 'line_fitting'
         figure;
         hold on;
         for j = 1:4
            h(j) = plot(s_list,SFhat_mean_fo(j,:)',linetype{j},'linewidth',2.5);
         end
         h(5:8) = yline(SF_mean_fo(:,end));
         legend(h([1,2,3,4]),'\mu_b','\mu_m','V_b','V_m')
         set(gca,'fontsize',18)
         xlabel('MCMC Sample Size')
         ylabel('First-order Sobol'' Index')
         ylim([0,.5])

        figure;
        hold on;
        for j = 1:4
            h(j) = plot(s_list,SFhat_mean_tot(j,:)',linetype{j},'linewidth',2.5);
        end
        h(5:8) = yline(SF_mean_tot(:,end));
        legend(h([1,2,3,4]),'\mu_b','\mu_m','V_b','V_m')
        set(gca,'fontsize',18)
        xlabel('MCMC Sample Size')
        ylabel('Total Sobol'' Index')
        ylim([0,.5])

        figure;
        hold on;
        for j = 1:4
            h(j) = plot(s_list,SFhat_var_fo(j,:)',linetype{j},'linewidth',2.5);
        end
        h(5:8) = yline(SF_var_fo(:,end));
        legend(h([1,2,3,4]),'\mu_b','\mu_m','V_b','V_m')
        set(gca,'fontsize',18)
        xlabel('MCMC Sample Size')
        ylabel('First-order Sobol'' Index')
        ylim([0,.6])

                figure;
        hold on;
        for j = 1:4
            h(j) = plot(s_list,SFhat_var_tot(j,:)',linetype{j},'linewidth',2.5);
        end
        h(5:8) = yline(SF_var_tot(:,end));
        legend(h([1,2,3,4]),'\mu_b','\mu_m','V_b','V_m')
        set(gca,'fontsize',18)
        xlabel('MCMC Sample Size')
        ylabel('Total Sobol'' Index')
        ylim([0,.6])
     case 'seir' 
         % F_mean
         figure;
        subplot(2,1,1)
        hold on;
         for j = 1:4
            h(j) = plot(s_list,SFhat_mean_tot(j,:)',linetype{j},'linewidth',2.5);
         end
        legend(h,'m_{log\mu}','m_{log\beta}','m_{log\sigma}','m_{log\gamma}','Location','East')
        set(gca,'fontsize',18)
        xlabel('MCMC Sample Size')
        ylabel('Total Sobol'' Index')
       % ylim([0,.5])

        subplot(2,1,1)
        hold on;
         for j = 1:4
            h(j) = plot(s_list,SFhat_mean_tot(j+4,:)',linetype{j},'linewidth',2.5);
         end
        legend('s^2_{log\mu}','s^2_{log\beta}','s^2_{log\sigma}','s^2_{log\gamma}','Location','SouthEast')
        set(gca,'fontsize',18)
        xlabel('MCMC Sample Size')
        ylabel('Total Sobol'' Index')
       % ylim([0,.5])

         % F_var
         figure;
        subplot(2,1,1)
        hold on;
         for j = 1:4
            h(j) = plot(s_list,SFhat_var_tot(j,:)',linetype{j},'linewidth',2.5);
         end
        legend(h,'m_{log\mu}','m_{log\beta}','m_{log\sigma}','m_{log\gamma}','Location','NorthEast')
        set(gca,'fontsize',18)
        xlabel('MCMC Sample Size')
        ylabel('Total Sobol'' Index')
        set(gca,'fontsize',18)

       % ylim([0,.5])
        subplot(2,1,1)
        hold on;
        for j = 1:4
            h(j) = plot(s_list,SFhat_var_tot(j+4,:)',linetype{j},'linewidth',2.5);
         end
        legend('s^2_{log\mu}','s^2_{log\beta}','s^2_{log\sigma}','s^2_{log\gamma}','Location','SouthEast')
        set(gca,'fontsize',18)
        xlabel('MCMC Sample Size')
        ylabel('Total Sobol'' Index')
       % ylim([0,.5])
     otherwise
         figure;
         subplot(2,1,1)
         for j = 1:nhyp
             plot(s_list,SFhat_mean_tot(j,:)',linetype{j},'linewidth',2.5)
         end
         set(gca,'fontsize',18)
         xlabel('MCMC Sample Size')
         ylabel('Total Sobol'' Index')

         subplot(2,1,2);
         for j = 1:nhyp
             plot(s_list,SFhat_var_tot(j,:)',linetype{j},'linewidth',2.5)
         end
         set(gca,'fontsize',18)
         xlabel('MCMC Sample Size')
         ylabel('Total Sobol'' Index')
 end