% Save system parameters - the same for all models (MFL, MFLGE, MFLGEGT)

function [parL] = savePars()

V_upperAirway = 0.105;  %L 0.105L for Upper Airway , 0.129L for MeshTrachea airway, 0.1571L for cylinder > assume zero since Vd = 0.15 - V_upperAirway
RuniversalGC = 8.31447;        % universal gas constant J/mol/K
T = 37 + 273.15;    % K
P = 101325;         % Pa
CONVERT_V2n = @(T,P) (RuniversalGC.*T)/P; % m^3/mol

ventilationType = 'breathing'; % or 'cardiogenic' or 'breathing' or 'apnea'

% cardiogenic ventialtion
% f_cos = 1/4.3;  % Hz (14 breaths per minute or 4.3s duration breaths)
% A_cos = (0.45); % L/s
f_cos = 1;    % Hz
A_cos = (1e-5)*1000; % L/s
%A_cos = 0; % L/s <- no cardiogenic oscillations
R_cos = 1.0760; % mmHg s/L <- when optimising for V_A amplitudes to match with CFD and when P_L is a function of P_m
%R_cos = 0.21;  % mmHg s/L <- when optimising for V_A amplitudes to match with CFD and when P_L was not a function of P_m
%R_cos = 1.0765;% mmHg s/L <- when optimising for Q_A amplitudes to match with CFD
parL.V_T = 0.4;          % L
parL.omega = (2.*pi)./5; % rad/s
Patm = 760;
radius = 0.01; % m
Length = 0.5; %m length pipe
Vcyl = (pi*(radius.^2))*Length;
mu_o = 2.03e-5; % kg/ms
mu_c = 1.47e-5; % kg/ms
mu_n = 1.73e-5; % kg/ms
%P_m = 760.003611; % Compare with Fluent single time-step
parL.P_m = 760;      % 760 mmHg
parL.p_w = 47;       % mmHg
parL.E = 2.5;        % mmHg/L
parL.R = 1;          % mmHg s/L
parL.Dc = 7.08e-3;   % L/s/mmHg
parL.Do = 3.5e-4;    % L/s/mmHg
parL.kp = 2.5; % Unit: mmHg/m - Conversion constant
parL.PL0 = 4.5; % Unit: mmHg - Constant related to pleural pressure
% Dbc = 3.16e-5;   % mol/s/L
% Dbo = 1.56e-5;   % mol/s/L
% ^- in Ben Tal 2010 (possibly calculated at T = 273.15 and P = 1 atm) does
% not significantly alter results despite Dbc being one order of magnitude
% lower in this listing than that calculated below
parL.Dbc = (parL.Dc*0.001)/(CONVERT_V2n(T,P)); % mol/s/L
parL.Dbo = (parL.Do*0.001)/(CONVERT_V2n(T,P)); % mol/s/L
parL.k1=2; % Unit: 1/s - Recoil rate constant of muscle
parL.k2=1; % Unit: m/s - Conversion constant
V_m = 0.105;    % L <- VdUA or upper airway dead-space volume < for variable Pm attempt 1

parL.V_0 = 0;
parL.M_o = ((-0.25/60)*0.001)/(CONVERT_V2n(T,P)); % L/min -> L/s -> m^3/s -> mol/s
parL.M_c = ((0.2/60)*0.001)/(CONVERT_V2n(T,P)); % L/min -> L/s -> m^3/s -> mol/s
parL.f_mo = 0.21;    % mouth fraction O2 
parL.f_mc = 0;       % mouth fraction CO2
%f_mo = 0.136800;     % mouth fraction O2  % Compare with Fluent single time-step
%f_mc = 0.052630;     % mouth fraction CO2 % Compare with Fluent single time-step
% f_mo = 0.1368;        % mouth fraction O2  % Compare with Fluent cylinder
% f_mc = 0.05263;       % mouth fraction CO2  % Compare with Fluent cylinder
% % >> Assume the fractions of CO2 and O2 are the same within the cylinder to the initial concentrations
f_ao = 0.1368;  % alveolar fraction O2
f_ac = 0.05263; % alveolar fraction CO2

parL.V_p = 0.07;     % L
%V_s = 0.2;      % L
parL.V_s = 4.93;      % L
parL.sigma_o = 1.4e-6; % mol/L/mmHg
parL.sigma_c = 3.3e-5; % mol/L/mmHg
parL.Th = 2e-3;      % mol/L
parL.delta = 10^(1.9); %
parL.l2 = 164e3;     % L/s/mol
parL.h = 10^(-7.4);  % mol/L
parL.r2 = 0.12;      % 1/s
parL.Q_CO = 0.0833;  % (5/60) L/s range:4-8L/min
parL.V_d = 0.15;% - V_upperAirway; % L
%pKa = 6.1;      % negative logarithm (base 10) of the acid dissociation constant of carbonic acid 
pKa = -log10(parL.r2./parL.l2);

% Diffusion rate of CO2 into body tissues from "The rate of rise of PaCO2 in the apneic anesthetized patient"
sigma_t = ((1.3*0.001*0.001)/(60*60))/(CONVERT_V2n(T,P)); % ml/mmHg/kg/hour -> m^3/mmHg/kg/s -> mol/mmHg/kg/s
BodyMass = 70; % kg

% for o2 dissociation only
% rR = 96.56e6;   % 1/mol/s
% lR = 26.8;      % 1/s
% rT = 6e6;       % 1/mol/s
% lT = 600;       % 1/s
% rl = 14;        % 1/s
% ll = 2.4e9;     % 1/s
% KR = rR./lR;
% KT = rT/lT;
% L = ll/rl;
parL.KR = 3.6*10^6; % L/mol
parL.KT = 10*10^3; % L/mol
parL.L = 171.2*10^6;
% ^- above listed in BenTal model - does not significantly alter results
save('SystemParameters');

return