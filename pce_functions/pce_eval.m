function y = pce_eval(pcdata,Xd)

t_set = Xd';

t_sz = size(t_set,1);

nhyp = pcdata.ndim;
multiIndex = pcdata.multiIndex;
nPCTerms = pcdata.nPCTerms;
psiMultiDSq = pcdata.psiMultiDSq;
nord = pcdata.nord;
pce = pcdata.pce;

rlz = zeros(t_sz,1); 
Psi_eval = zeros(t_sz,nPCTerms);
for i = 1:t_sz
    xi = t_set(i,:); x = 2 * xi - 1;
    for j = 1:nPCTerms
        mind = multiIndex(j,:);
        vals = zeros(nhyp,1);
        for k = 1:nhyp
            vals(k) = leg_eval(mind(k),x(k));
        end
        Psi_eval(i,j) = prod(vals);
    end
end

y = Psi_eval * (pce);

