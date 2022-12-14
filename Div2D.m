function [div] = Div2D(ncell,nei,ind,edges)

% divergence matrix

div = sparse(ncell,nei);

for e=1:nei
    K = edges(ind.internal(e),3);
    L = edges(ind.internal(e),4);
    div(K,e) = edges(ind.internal(e),6);
    div(L,e) = -edges(ind.internal(e),6);
end