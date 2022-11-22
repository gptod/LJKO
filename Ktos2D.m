function [Rs] = Ktos2D(ncell,nei,ind,edges,gradphi)

% upwind reconstruction operator from cells to edges

Rs = sparse(nei,ncell);

for e=1:nei
    K = edges(ind.internal(e),3);
    L = edges(ind.internal(e),4);
    % upwind
    Rs(e,K) = Rs(e,K) + real(gradphi(e)>=0);
    Rs(e,L) = Rs(e,L) + real(gradphi(e)<0);
end