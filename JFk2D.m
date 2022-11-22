function [JF] = JFk2D(ind,edges,cc,Ms,Mx,gradx,divx,uk,tau,ddE)


ncell = size(cc,1);
nei = length(ind.internal);
phik = uk(1:ncell);
rhok = uk(ncell+1:2*ncell);
% if the energy is not local this has to be modified
ddEk = spdiags(ddE(rhok,cc(:,1),cc(:,2)),0,ncell,ncell); 

gradphik = gradx*phik;
Rs = Ktos2D(ncell,nei,ind,edges,-gradphik);
RK = Mx\(Rs'*Ms);

% derivatives of the optimality conditions in phi
JFpp = -divx*spdiags((Rs*rhok),0,nei,nei)*gradx;
JFpr = Mx/tau-divx*spdiags(gradphik,0,nei,nei)*Rs;

% derivatives of the optimality conditions in rho
JFrp = Mx/tau+Mx*RK*(spdiags(gradphik,0,nei,nei)*gradx);
JFrr = -Mx*ddEk/tau;


JF.pp = JFpp; JF.pr = JFpr;
JF.rp = JFrp; JF.rr = JFrr;

  

end