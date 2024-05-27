
 function dy = SEIR_rhs(t,y,q,N)

  mu = q(1); % birth rate
  beta = q(2); % Infection rate
  sigma = q(3); % latency rate
  gamma = q(4); % Recovery rate

  
  dy = [mu*N - beta * y(1)*y(3)/N - mu*y(1); %dS
        beta*y(1)*y(3)/N - (sigma+mu)*y(2);  %dE
        sigma*y(2) - (gamma+mu)*y(3);        %dI
        gamma*y(3) - mu*y(4)];               %dR             
                                         