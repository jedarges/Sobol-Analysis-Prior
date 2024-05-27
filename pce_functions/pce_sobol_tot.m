function ST = pce_sobol_tot(pce,pcdata,i)
%% Total Sobol' index S_i^tot from PC coefficients

multiIndex = pcdata.multiIndex;
psiMultiDSq = pcdata.psiMultiDSq;

multiIndex = multiIndex(2:end,:);
u = pce(2:end).^2;
v = psiMultiDSq(2:end);

ind = multiIndex(:,i) > 0;

ST = dot(u(ind),v(ind))/dot(u,v);