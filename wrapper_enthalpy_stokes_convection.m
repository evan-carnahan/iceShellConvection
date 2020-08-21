% Authors: Evan Carnahan (evan.carnahan@utexas.edu), Marc "Grim" Hesse, Jake
% Jordan and Natalie "Wolf" Wolfenbarger
% Description: Wrapper to run enthalpy_stokes_convection.m

%%
clear all; close all
res = 20;
T_t = 100;
cmplx = 'WbCpNl';
pi_tide = 2.5;

%% temperature calculations
T_b = 273; %K
DT = T_b - T_t;

%% Specific heat - this is just used internally 
a = 185;
b = 7.037;
c_fun = @(nonT) a+b*(DT*nonT+T_t); %J/(kg K)
c_pi = c_fun(1); % J/(kg K)

%% Density variation - non-linear and linear
if strcmp(cmplx(5:6),'Nl')
    rho_fun = @(nonT) ice_density(nonT*DT+T_t,611.657);
    rho_b = rho_fun(1);
elseif strcmp(cmplx(5:6),'Ln')
    alpha = 1.6e-4; %1/K
    rho_b = 917; %kg/m^3
    rho_fun = @(nonT) rho_b*(1-alpha*DT*(nonT-1));
end
rho_t = rho_fun(0);
DRho = rho_b - rho_t;
% function for nondimensional denisty
rhoPrime_fun = @(nonT) (rho_fun(nonT) - rho_b)/DRho;

%% Thermal conductivity
if strcmp(cmplx(1:2),'He')
    % Hobbs, 1974 
    kappa_fun = @(nonT) 0.4685+488.12./(nonT*DT+T_t); %W/(m K)
elseif strcmp(cmplx(1:2),'Rb')
    % Rabin, 2001 
    kappa_fun = @(nonT) 2135./(nonT*DT+T_t).^1.235;%W/(m K)
elseif strcmp(cmplx(1:2),'Cn')
    % Constant
    kappa_fun = @(nonT) 2.26*ones(size(nonT));
elseif strcmp(cmplx(1:2),'Wb')
    % Model from Carnahan, 2020
     kappa_fun = @(nonT) 612./(nonT*DT+T_t);
end
kappa_b = kappa_fun(1); % W/(m K)
% function for nondimensional thermal conductivity
nonKappa_fun = @(nonT) kappa_fun(nonT)/kappa_b;

%% Viscosity
E_a = 50e3;% J mol^-1 
R = 8.314; % J K^-1 mol^-1 
Apar = E_a/R/T_b;
barrVisc = @(nonT) exp(Apar*(T_b./(DT.*nonT+T_t)-1));

%% Enthalpy conversions

% for constant specific heat
if strcmp(cmplx(3:4),'Cn')
    % convert dimensionless temperature to dimensionless enthalpy
    nonH_fun = @(nonT) nonT - 1;
    % convert dimensionless enthalpy to dimensionless temperature
    nonT_fun = @(nonH) nonH + 1;
    % calculate dimesnionless inverse volumetric heat capacity
    dnonT_dnonH_fun = @(nonH) ones(size(nonH));    
% for linear specific heat
elseif strcmp(cmplx(3:4),'Cp')
    % convert dimensionless temperature to dimensionless enthalpy
    nonH_fun = @(nonT) 1/(c_pi*DT)*(a*DT*(nonT-1)+b/2*((nonT*DT+T_t).^2-T_b^2));
    % convert dimensionless enthalpy to dimensionless temperature
    nonT_fun = @(nonH) (-a-b*T_t+sqrt(c_pi)*sqrt(2*b*DT*nonH + c_pi))/(b*DT);
    % calculate dimesnionless inverse volumetric heat capacity
    dnonT_dnonH_fun = @(nonH) max(1./sqrt(1+(2*b*DT*nonH)./c_pi),1);
end

%% first spin up to SS conductive profile
Ra = 0; % for spin up of conductive profile
loadFile = 'spinup';
iterArray = [1e3,1e3+1];
[H,T,str] = enthalpy_stokes_convection(loadFile,...
    res,Ra,pi_tide,barrVisc,rhoPrime_fun,nonKappa_fun,dnonT_dnonH_fun,nonH_fun,nonT_fun,cmplx,T_t,iterArray); 

%% perturb and convect
Ra = 8e6;
loadFile = join(['convSS_res' num2str(res) '_Ra0_tidal'...
    num2str(pi_tide) '_Ts' num2str(T_t) '_' cmplx '.mat'],"");
iterArray = [100,1e4];
[H,T,str] = enthalpy_stokes_convection(loadFile,...
    res,Ra,pi_tide,barrVisc,rhoPrime_fun,nonKappa_fun,dnonT_dnonH_fun,nonH_fun,nonT_fun,cmplx,T_t,iterArray);



