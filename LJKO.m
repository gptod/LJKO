clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Upwind finite volume LJKO scheme %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% LOAD MESH
load('sq_mesh2')

%% DEFINITION OF THE GRADIENT FLOW:

% NB: if the energy is not local the construction of the hessian
% in JFk2D must be adapted!
% NB: in general for other energies the functions Fk2D and JFk2D
% may need to be tailored in the call of the energy term

% Fokker-Planck equation
% V =@(x,y) x.^2;
V =@(x,y) zeros(size(x));
E =@(rho,x,y) rho.*log(rho)-rho.*V(x,y);
dE =@(rho,x,y) log(rho)+1-V(x,y);
ddE =@(rho,x,y) 1./rho;

% % Porous media equation:
% % definition of the energy
% V =@(x,y) 0.5*((x-0.5).^2+(y-0.5).^2);
% m = 2; % m > 1
% E =@(rho,x,y) rho.^m/(m-1)+rho.*V(x,y); 
% dE =@(rho,x,y) m*rho.^(m-1)/(m-1)+V(x,y);
% ddE =@(rho,x,y) m*rho.^(m-2);


% INITIAL CONDITIONS:
initial =@(x,y) exp((-(x-0.5).^2-(y-0.5).^2)/0.05);
%initial =@(x,y) x<=0.5;
rho0 = initial(cc(:,1),cc(:,2)); 
mass = sum(rho0.*area); 


%% PARAMETERS

T = 1; %integration time
tau = 0.01; %initial time step 

esin = 1; % set 1 if the energy is singular for vanishing measures, 0 else
          % avoid the Newton method to compute intermediate negative values
eloc = 0; % set 1 to solve via the Schur complement, 0 else
          % convenient only if the hessian of the energy is simple to
          % invert (ex. if the energy is local and therefore the hessian is diagonal)

verb = 1; % verbosity level: {0,1}

% Newton scheme's parameters:
kmax = 15; %maximum number of iteration
kmin = 5; %minimum number of iteration
eps = 1e-6; % tolerance


%% INITIALIZATION OF THE SCHEME

% assemble matrices
Mx = spdiags(area,0,ncell,ncell);
ds = edges(ind.internal,5).*edges(ind.internal,6);
Ms = spdiags(ds,0,nei,nei);
divx = Div2D(ncell,nei,ind,edges);
gradx = -Ms\divx';

% initial conditions
rhot = rho0;
phit =zeros(ncell,1);
uk = [phit; rhot];


time = 0;
itn = 0;

Energy = sum(area.*E(rhot,cc(:,1),cc(:,2)));
ts = time;

%% START TIME INTEGRATION

while abs(time-T) > 1e-10
    
    time = time+tau;
    if time>T && abs(time-T)>1e-10
        time = time-tau;
        tau = T-time;
        time = time+tau;
    end
    itn = itn+1;

    itk = 0;
    errs = [];
    while true
        
        itk = itk+1;
        
        if itk>kmax %|| any(uk(ncell+1:end)<0)
            time = time-tau;
            tau = tau/2;
            time = time+tau;
            uk(1:ncell) = phit;
            uk(ncell+1:end) = rhot;
            itk = 0;
            fprintf('%10s %4i %21s \n','iteration ',itn,', decrease time step')
            fprintf('%21s %1.5e \n','actual time step is: ',tau)
            continue
        end
        
        Fk = Fk2D(ind,edges,cc,Ms,Mx,gradx,divx,rhot,uk,tau,dE);
        errs = [errs; norm([Fk.p;Fk.r],2)];
        
        if verb>0
            fprintf('%11s %4i %18s %4i %7s %1.4e \n','Time step: ',itn,'Newton iteration: ',itk,'Error: ',errs(end))
        end
        
        if errs(end) < eps
            if itk < kmin
                tau = 1.2*tau;
                fprintf('%10s %4i %20s \n','iteration ',itn,', increase time step')
                fprintf('%21s %1.5e \n','actual time step is: ',tau)
            end
            break
        end

        JFk = JFk2D(ind,edges,cc,Ms,Mx,gradx,divx,uk,tau,ddE);
        if eloc==1
            S = JFk.pp-JFk.pr*(JFk.rr\JFk.rp);
            dp = -S\(Fk.p-JFk.pr*(JFk.rr\Fk.r));
            dr = JFk.rr\(-Fk.r-JFk.rp*dp);
            dk = [dp;dr];
        else
            JFk = [JFk.pp JFk.pr; JFk.rp JFk.rr];
            dk = -JFk\[Fk.p;Fk.r];
        end

        alfak = 1;
        if esin==1
            while any(uk(ncell+1:2*ncell)+alfak*dk(ncell+1:2*ncell)<0)
                alfak = 0.5*alfak;
            end
        end
        uk = uk + alfak*dk;
        
    end
    
    fprintf('%11s %4i %22s %4i %7s %1.4e \n','Time step: ',itn,'Total Newton iterations: ',itk,'Error: ',errs(end))
    
    rho_all = [uk(ncell+1:2*ncell);rhot];
    
    phit = uk(1:ncell);
    rhot = uk(ncell+1:2*ncell);
    
    Energy = [Energy; sum(area.*E(rhot,cc(:,1),cc(:,2)))];
    ts = [ts; time];
    
    % plot
    Frho = scatteredInterpolant(cc,rhot);
    Frho.ExtrapolationMethod = 'nearest';
    Zrho = Frho(nodes(:,1),nodes(:,2));
    figure(1)
    trisurf(cells(:,2:end),nodes(:,1),nodes(:,2),Zrho)
    %axis([0 1 0 1 0 1.2])
    shading interp
    view([0 90])
    colormap('jet')
    colorbar
    axis square
    axis off

end

% Energy decay
figure(2)
%plot(ts,Energy,'-*')
semilogy(ts,Energy,'-*')




