         figure;
         subplot(2,1,1)
         for j = 1:nhyp
             plot(spf_set, Stot_mean_pf(j,:)',linetype{j},'linewidth',2.5)
         end
         set(gca,'fontsize',18)
         xlabel('Design Sample Size')
         ylabel('Total Sobol'' Index')

         subplot(2,1,2);
         for j = 1:nhyp
             plot(spf_set,Stot_var_pf(j,:)',linetype{j},'linewidth',2.5)
         end
         set(gca,'fontsize',18)
         xlabel('Design Sample Size')
         ylabel('Total Sobol'' Index')