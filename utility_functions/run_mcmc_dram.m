for iter = 1:tot
    p = m_old + (beta*gamma_0 + (1-beta)*gamma_n) * randn(npar,1);
    Pi_new = logtarget(p);
    alpha = Pi_new - Pi_old;
    rando = rand;
    if alpha < -2*log(rando)
        Pi_old = Pi_new;
        m_old = p; qoi_old = qoi(m_old);
        if iter > burnin+adapt_iter
            accept = accept + 1;
        end
    else  % Delayed rejection step
        Y_list = zeros(npar,DR+1); Y_list(:,1) = p;
        Pi_list = zeros(DR+1,1); Pi_list(1) = Pi_new;
        for k = 1:DR
            %  Yi = beta * (Y_list(k)+sc(k)*gamma_0*randn(npar,1)) + (1-beta) * (Y_list(k)+sc(k)*gamma_n*randn(npar,1));
            Yi = Y_list(:,k) + sc(k)*(beta*gamma_0+ (1-beta)*gamma_n) * randn(npar,1);
            Y_list(:,k+1) = Yi;
            Pi_new = logtarget(Yi);
            Pi_list(k+1) = Pi_new;
            if Pi_list(k+1) < Pi_old
                Pi_old = Pi_list(k+1);
                m_old = Yi;
                if iter > burnin+adapt_iter
                    accept = accept + 1;
                end
                qoi_old = qoi(m_old);
                break
            elseif Pi_list(k+1) > Pi_old
                break
            else
                Pi_st = min(Pi_list);
                pb = (exp(-.5 * Pi_list(k+1)) - exp(-.5 * Pi_st)) / (exp(-.5 * Pi_old) - exp(-.5 * Pi_st));
                rando = rand;
                if pb > rando
                    Pi_old = Pi_list(k+1);
                    if iter > burnin+adapt_iter
                        accept = accept + 1;
                    end
                    m_old = Yi;
                    qoi_old = qoi(m_old); break
                end
            end
        end
    end
    chain(iter,:) = m_old;
    qoi_chain(iter) = qoi_old;

    %% Adaptive Step
    if mod(iter,Ad) == 0 && iter > burnin && iter < burnin+adapt_iter
        Cnt_ad = Cnt_ad + 1;
        batch = Cnt_ad * Ad;
        beta = 1 / batch;
        Sign = sd * (cov(chain(1:iter,:)) + sqrt(eps) * eye(npar));
        gamma_n = chol(Sign)';
    end
    if mod(iter,Ad) == 0 && iter > burnin+adapt_iter
       save("pfconv_gsa.mat","chain","qoi_chain")
    end
end