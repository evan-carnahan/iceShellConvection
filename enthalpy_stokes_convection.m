function [H,T,str] = enthalpy_stokes_convection(loadFile,res,...
    Ra,pi_tide,viscF,rhoPrimeF,nonKappaF,dnonT_dnonHF,nonHF,nonTF,cmplx,T_s,iterArr)
% Authors: Evan Carnahan (evan.carnahan@utexas.edu), Marc "Grim" Hesse, Jake
% Jordan and Natalie "Wolf" Wolfenbarger
%
% Description: function for dimensionless enthalpy convection from Carnahan
% et al., 2020
%
% Variables: @ loadFile (str or char) - either "spinup" to load linear
%                       conductive profile or path to conductive profile
%            @ res (int) - grid resolution will be multiplied by aspect
%                       ration and level of anisotropy in y-direction below
%            @ Ra (float) - Rayleigh number of simulation
%            @ pi_tide (float) - tidal heating number of simulation
%            @ viscF (anon. func.) - dimensionless viscosity function
%            @ rhoPrimeF (anon. func.) - dimensionless bouyancy function
%            @ nonKappaF (anon. func) - dimensionless thermal conducitivity
%                                       function
%            @ dnonT_dnonHF (anon. func.) - dimensionless inverse volumetric 
%                                           heat capacity 
%            @ nonHF (anon. func.) - convert dimensionless temperature to
%                                   dimensionless enthalpy
%            @ nonTF (anon. func.) - convert dimensionless enthalpy to
%                                   dimensionless temeprature
%            @ cmplx (char) - string used in wrapper to create anon. func.,
%                             just used for saving
%            @ T_s (int) - surface temperature used in wrapper, just used for 
%                          saving
%            @ iterArr (array) - [minimum iterations, maximum iterations]
%
% Outputs: @ H (array) - dimensionless enthalpy
%          @ T (array) - dimensionless temperature
%          @ str (array) - stored energy
%

%%
% theta method for time stepping, 0.5 is Crank-Nichelson
theta = 0.5;
%Matricies are ill conditioned due to large viscosity contrast
warning off
%whether to generate plots through simulation
pl = true;
% Ra must be positive value
if Ra < 0
    Ra = abs(Ra);
end
% Assign max and minimum iterations
minIter = iterArr(1);
maxIter = iterArr(2);
% % maximum iterations to run simulation
% it = maxIter+1;


%% Set up grid
% aspect ratio of simulation
lam = 2;
% anisotropy of grid in y-direction
anis = 2;

% Build cartesian grid for noperator building
Gridp = struct('xmin',{},'xmax',{},'Nx',{},'ymin',{},'ymax',{},'Ny',{});
Gridp(1).xmin = 0; Gridp(1).xmax = 1*lam; Gridp(1).Nx = res*lam;
Gridp(1).ymin = 0; Gridp(1).ymax = 1; Gridp(1).Ny = res*anis;

% Build staggerd stokes grid
Grid  = build_stokes_grid(Gridp);
% Build discrete stokes operators - see Carnahan 2020 supplement for
% explanation
[D,Edot,Dp,Gp,I,Gyy]=build_stokes_ops(Grid);

% Build identity and zero matricies for forming full linear operator
Zp = zeros(Grid.p.N);
Ip = speye(Grid.p.N);

% Build grids for plotting
[X,Y] = meshgrid(Grid.p.xc,Grid.p.yc);
[Xc,Yf] = meshgrid(Grid.p.xc,Grid.p.yf);

% Variable indicies
linInds = find(Gyy > 0);
[row,~] = ind2sub(size(Gyy),linInds);

%% load in starting profiles
% Temperature perturbation
tempPert = cos(2*pi*Grid.p.xc);
fullPert = repmat(tempPert',Grid.p.Ny,1);

% spinup conductive profile or load in conductive profile and perturb
if strcmp(loadFile,'spinup')
    T = 1 - Y;
    T = T(:);
    H = nonHF(T);
else
    load(loadFile,'T');
    tempPert =  cos(2*pi*Grid.p.xc); %cos(pi*Grid.p.xc); %
    fullPert = 0.001*repmat(tempPert',Grid.p.Ny,1);
    T = T + fullPert(:);
    % change enthalpy in accordance with temperature perturbation
    H = nonHF(T);
    
end
Hplot = reshape(H,Grid.p.Ny,Grid.p.Nx);

%% build enthalpy boundary condition - based off of set temepratures
T1 = 1;
T0 = 0;
Param = struct('H',{},'g',{},'dof_dir',{});
Param(1).H = struct('dof_dir',{},'dof_f_dir',{},'g',{},'dof_neu',{},'dof_f_neu',{},'qb',{},'dof_out',{});
% Convert temperature to enthalpy BC
H1 = nonHF(T1);
H0 = nonHF(T0);
Param.H(1).dof_dir = [Grid.p.dof_ymin;Grid.p.dof_ymax];
Param.H(1).dof_f_dir = [Grid.p.dof_f_ymin;Grid.p.dof_f_ymax];
Param.H(1).g = [H1*ones(length(Grid.p.dof_ymin),1); H0*ones(length(Grid.p.dof_ymax),1)];
Param.H(1).dof_neu = [Grid.p.dof_xmin;Grid.p.dof_xmax];
Param.H(1).dof_f_neu = [Grid.p.dof_f_xmin;Grid.p.dof_f_xmax];
Param.H(1).qb = 0*[Grid.p.dof_f_xmin;Grid.p.dof_f_xmax];
[BH,NH,fn_H] = build_bnd(Param.H,Grid.p,Ip);
Param.H(1).dof_out = [Grid.p.dof_ymin];

%% build boundary conditions
Param(1).dof_dir =  [...
    Grid.x.dof_xmax;...           %set x_max x-vel
    Grid.x.dof_xmin;...           %set x_min x-vel
    Grid.x.N+Grid.y.dof_ymin;...  %set y_min y-vel
    Grid.x.N+Grid.y.dof_ymax;...  %set y_max y-vel
    Grid.p.Nf+1];                 %set pressure

Param(1).g = 0*Param.dof_dir;
Param.g(end) = 0;
B = I([Param.dof_dir],:);
N = I;
N(:,[Param.dof_dir]) = [];

for i =1:(maxIter+1)
    Tplot= reshape(T,Grid.p.Ny,Grid.p.Nx);
    %% density forcing claculation
    rhoPrime = rhoPrimeF(T);
    rhoPlot = reshape(rhoPrime,Grid.p.Ny,Grid.p.Nx);
    rhoFace = diag(comp_mean(rhoPlot,1,1,Grid.p));
    rhoY = rhoFace(Grid.p.Nfx+1:Grid.p.Nf);
    fsVec = -Ra*rhoY;
    % forcing function applies on y-velocities
    fs = [zeros(Grid.p.Nfx,1); fsVec; zeros(Grid.p.N,1)];
    
    %% Temperature dependent viscosity
    
    %Gxx variable viscosity matrix
    nxxVec = zeros(Grid.x.Nfx,1);
    nxxVec(Grid.x.Ny+Grid.p.dof) = Tplot;
    
    %Gyy variable viscosity matrix
    nyyVec = zeros(Grid.y.Nfy,1);
    nyyVec(row) = Tplot;
    % Harmonic average of temperature to corners
    ncVec = comp_mean_corners(Tplot,-1,Grid.p);
    
    % full vector of relavent temperatures for viscosity calculation
    tempVec = [nxxVec;nyyVec;ncVec];
    viscVec = viscF(tempVec);
    viscMat = spdiags(viscVec,0,length(viscVec),length(viscVec));
    %% linear operator and stokes solution
    % form divergence matrix and linear operator
    tau = D*2*viscMat*Edot;
    L = [tau,Gp;
        Dp,Zp];
    
    % solve linear operator
    u = solve_lbvp(L,fs,B,Param.g,N);
    
    % make vectors of stokes solution
    vx = u(1:Grid.p.Nfx);
    vy = u(Grid.p.Nfx+1:(Grid.p.Nfx+Grid.p.Nfy));
    vm = [vx;vy];
    vmax= max((abs(vm)));
    p  = u((Grid.p.Nfx+Grid.p.Nfy+1):end);
    
    %% Tidal heating forcing
    tideVisc = viscF(T);
    R_prime = 2*tideVisc./(1+tideVisc.^2);
    R_tide = R_prime*pi_tide;
    
    %% Advection Diffusion Coefficients
    % calculate thermal conductivity
    kappaPrime = nonKappaF(T);
    kappaPrimePlot = reshape(kappaPrime,Grid.p.Ny,Grid.p.Nx);
    kappaFace = comp_mean(kappaPrimePlot,1,1,Grid.p);
    
    % calculate inverse volumetric specific heat
    dT_dH = dnonT_dnonHF(H);
    dT_dH_plot = reshape(dT_dH,Grid.p.Ny,Grid.p.Nx);
    dT_dHFace = comp_mean(dT_dH_plot,1,1,Grid.p);
    
    % Courant-Fredrichs-Lewy condition for time step
    dt = min([0.5*Grid.p.dx^2, Grid.p.dx/vmax,0.5*Grid.p.dy^2, Grid.p.dy/vmax])*0.8;
    tVec(i) = dt;
    
    %% Advection Diffusion Solution
    HOld = H;
    % upwinded advective flux
    AH = build_adv_op(vm,H,dt,Gp,Grid.p,Param.H,'mc');
    
    % Explicit and implicit operators for time stepping
    L_T_E = Ip - (1-theta)*dt*(Dp*AH-Dp*kappaFace*dT_dHFace*Gp);
    L_T_I = Ip + (theta)*dt*(Dp*AH-Dp*kappaFace*dT_dHFace*Gp);
    
    % Add in source terms
    %     RHS_T = L_T_E*H + (R_tide+fn_H)*dt;
    RHS_T = L_T_E*H + R_tide*dt;
    
    % Solve enthalpy transport problem
    H = solve_lbvp(L_T_I,RHS_T,BH,Param.H.g,NH);
    %% Switch back to temperature
    T =  nonTF(H);
    %Set temperature back to melting point
    T(T > 1) = 1;
    
    %% redimensionaltion of fluxes in domain
    % residual of transport problem
    res_trans = @(H,HOld) (L_T_I*H - L_T_E*HOld)/dt - R_tide;
    flux_trans = @(H,HOld) (AH - Gp)*(theta*HOld+(1-theta)*H);
    % corrected flux on boundary
    q_nonH = comp_flux_gen(flux_trans,res_trans,H,Grid.p,Param.H,HOld);
    
    % surface heat flux
    surfHeatF = q_nonH(Grid.p.dof_f_ymax);
    totSurfHeatFlux(i) = sum(surfHeatF .* Grid.p.A(Grid.p.dof_f_ymax));
    
    % bottom heat flux
    botHeatF = q_nonH(Grid.p.dof_f_ymin);
    totBotHeatFlux(i) = sum(botHeatF .* Grid.p.A(Grid.p.dof_f_ymin));
    
    % tidal heat flux
    totTide(i) = sum(R_tide.*Grid.p.V);
    
    % stored heat flux
    totStr(i) = sum((H - HOld).*Grid.p.V)/dt;
    
    % total flux- should sum to ~0
    totFlux(i) = totTide(i) - totStr(i) + totBotHeatFlux(i) - totSurfHeatFlux(i);
    
    % potential plots
    if pl == true
        if mod(i,20) == 0
            figure(2)
            [PSI,~,~] = comp_streamfun(vm,Grid.p);
            set(gcf, 'Position', [50 50 1500 600]);
            
            % plot dimensionless temperature and streamfunction
            subplot(3,3,1)
            cla;
            hold on
            axis equal
            contourf(X,Y,Tplot,40,'linestyle','none'),view(2),hold on
            c = colorbar('NorthOutside');
            cStrVal = linspace(min(PSI(:)),max(PSI(:)),10);
            contour(X,Y,Tplot,'r','LevelList',1)
            contour(Grid.p.xf,Grid.p.yf,PSI,'k','LevelList',cStrVal);
            caxis([min(Tplot(:)) max(Tplot(:))]);
            c.Label.String = 'Temperature, 1';
            xlabel('x-dir, 1')
            ylabel('z-dir, 1')
            
            % Plot dimensionless enthalpy
            subplot(3,3,2)
            cla;
            hold on
            axis equal
            contourf(X,Y,reshape(H,Grid.p.Ny,Grid.p.Nx),40,'linestyle','none'),view(2),hold on
            c = colorbar('NorthOutside');
            xlabel('x-dir, 1')
            ylabel('z-dir, 1')
            c.Label.String = 'Enthalpy, 1';
           
            % Plot dimensionless tidal heating
            subplot(3,3,3)
            cla;
            hold on
            axis equal
            contourf(X,Y,reshape(R_prime,Grid.p.Ny,Grid.p.Nx),40,'linestyle','none'),view(2),hold on
            c = colorbar('NorthOutside');
            xlabel('x-dir, 1')
            ylabel('z-dir, 1')
            c.Label.String = 'Tidal energy, 1';
            
            % plot dimensionless thermal conductivity
            subplot(3,3,4)
            cla;
            hold on
            axis equal
            contourf(X,Y,kappaPrimePlot,40,'linestyle','none'),view(2),hold on
            c = colorbar('NorthOutside');
            xlabel('x-dir, 1')
            ylabel('z-dir, 1')
            c.Label.String = 'Thermal cond, 1';
            
            % plot dimensionless inverse volumetric heat capacity
            subplot(3,3,5)
            cla;
            hold on
            axis equal
            contourf(X,Y,reshape(dT_dH,Grid.p.Ny,Grid.p.Nx),40,'linestyle','none'),view(2),hold on
            c = colorbar('NorthOutside');
            xlabel('x-dir, 1')
            ylabel('z-dir, 1')
            c.Label.String = 'dT/dH, 1';
            
            % plot dimensionless viscosity
            subplot(3,3,6)
            cla;
            axis equal
            cenVisc = reshape(viscVec(Grid.p.Ny+1:Grid.p.Ny+Grid.p.N),Grid.p.Ny,Grid.p.Nx);
            contourf(X,Y,log10(cenVisc),40,'linestyle','none'),view(2),hold on
            c = colorbar('NorthOutside');
            c.Label.String = 'Log. Viscosity, 1';
            xlabel('x-dir, 1')
            ylabel('z-dir, 1')
            
            % Plot each energy flux
            subplot(3,3,7)
            cla;
            hold on
            plot(totBotHeatFlux)
            plot(-totSurfHeatFlux)
            plot(-totStr)
            plot(totTide)
            plot(totFlux)
            ylabel('Total energy fluxs (positive into domain), 1')
            legend('Bottom flux','Surface flux','Stored','Tidal','Total','Location','NorthWest');
                        
            % plot stored energy
            subplot(3,3,8)
            cla;
            plot(cumsum(tVec),totStr)
            ylabel('Stored energy flux, 1');
            xlabel('t_{prime}');
            
             % Plot velocity over iterations
            subplot(3,3,9)
            cla;
            hold on
            plot(vRmsSW)
            ylabel('RMS velcotiy')
        end
    end
    
    %% make arrays for output
    vx(Grid.x.dof_xmax) = [];
    vy(Grid.y.dof_ymax) = [];
    vMagSW = sqrt(vx.^2+vy.^2);
    vRmsSW(i) = sum(vMagSW)/length(vMagSW);
    vMax(i) = max(vMagSW(:));
    
    % criteria for end of run -> used to check for onset of convection
    % would need to be changed to energy storage cut off for convective
    % steady state
    if (i > minIter && vRmsSW(end) < 1e-6 || i > maxIter)
        str = totStr;
        % path for saved file
        save(join(['convSS_res' num2str(res) '_Ra' num2str(Ra) ...
            '_tidal' num2str(pi_tide) '_Ts' num2str(T_s) '_' cmplx '.mat'],""),'H','T','str','totBotHeatFlux','totSurfHeatFlux','vx','vy','tVec','vRmsSW','surfHeatF','botHeatF','vMax');
        return
    end
end
end