function S = pce_sobol_fo(pce,pcdata,i)
%% First order Sobol' index S_i from PC coefficients


multiIndex = pcdata.multiIndex;
psiMultiDSq = pcdata.psiMultiDSq;

multiIndex = multiIndex(2:end,:);
u = pce(2:end).^2;
v = psiMultiDSq(2:end);

multiIndex(:,i) = [];
ind = sum(multiIndex,2) == 0;

S = dot(u(ind),v(ind))/dot(u,v);