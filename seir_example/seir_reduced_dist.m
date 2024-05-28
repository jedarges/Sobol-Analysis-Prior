switch choose_example
    case 'seir'
        Xi_R = rand(nhyp-3,N);
        Xi_red = zeros(nhyp,N);
        Xi_red(1:3,:) = Xi_R(1:3,:);
        Xi_red(6:7,:) = Xi_R(4:end,:);
        Xi_red(4:5,:) = 0.5 * ones(2,N);
        Xi_red(8,:) = 0.5 * ones(1,N);


        weights_red = weight_compute(Xi_red,chain_reduced,data);
        C_red = sum(weights_red,1); % approximated normalizing constants
        Fmean_red = (qoi_chain_reduced' * weights_red) ./ C_red;
        Fvar_red = ((qoi_chain_reduced.^2)' * weights_red) ./ C_red - Fmean_red.^2;

        Fatmap_red = zeros(N,1);
        for k = 1:N
            xi = Xi_red(:,k);
            prior = prior_select(data,xi);
            map_target = @(theta) ssfun(theta) + prior.logprior(theta);
            [m_map,~] = fminunc(map_target,data.mean,options);
            Fatmap_red(k) = qoi(m_map);
        end

        % Plot distributions of Fmean and Fvar
        [dmean_red,qmean_red] = ksdensity(Fmean_red);
        [dvar_red,qvar_red] = ksdensity(Fvar_red);
        [datmap_red,qatmap_red] = ksdensity(Fatmap_red);

        figure;
        hold on;
        plot(qmean,dmean,'linewidth',2.5)
        plot(qmean_red,dmean_red,'--','linewidth',2.5)
        set(gca,'fontsize',18)
        xlabel('F_{mean}(\xi)')
        ylabel('Density')
        legend('Full','Reduced')
        hold off;

        figure;
        hold on;
        plot(qvar,dvar,'linewidth',2.5)
        plot(qvar_red,dvar_red,'--','linewidth',2.5)
        set(gca,'fontsize',18)
        xlabel('F_{var}(\xi)')
        ylabel('Density')
        legend('Full','Reduced')
        hold off;

        [datmap,qatmap] = ksdensity(Fatmap);

        figure;
        hold on;
        plot(qatmap,datmap,'linewidth',2.5)
        plot(qatmap_red,datmap_red,'--','linewidth',2.5)
        set(gca,'fontsize',18)
        xlabel('F_{atMAP}(\xi)')
        ylabel('Density')
        legend('Full','Reduced')
        hold off;
    otherwise
end