function [F] = Fk2D(ind,edges,cc,Ms,Mx,gradx,divx,rhot,uk,tau,dE)

ncell = size(cc,1);
nei = length(ind.internal);
phik = uk(1:ncell);
rhok = uk(ncell+1:2*ncell);
dEk = dE(rhok,cc(:,1),cc(:,2));

gradphik = gradx*phik;
Rs = Ktos2D(ncell,nei,ind,edges,-gradphik);
RK = Mx\(Rs'*Ms);

% optimality conditions in phi
Fp = Mx*(rhok-rhot)/tau-divx*((Rs*rhok).*gradphik);

% optimality conditions in rho
Fr = Mx*(phik-dEk)/tau+0.5*Mx*RK*(gradphik).^2;

F.p = Fp;
F.r = Fr;


end

