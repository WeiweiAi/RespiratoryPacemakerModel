% Inital variables

function [x0] = initialVariables()
% x0 = column vector of initial variables
% Initial values for constant pH model during apnea using steady state
% breathing output as the initial apnea conditions
x1 = 759.9839;  % P_A (alveolar mmHg)
x2 = 0.1349;    % f_o (alveolar)
x3 = 0.0651;    % f_c (alveolar)
x4 = 0.1746;    % f_do (dead-space)
x5 = 0.0288;    % f_dc (dead-space)
x6 = 36.4882;   % p_Vo (venous mmHg)
x7 = 84.2121;   % p_Ao (arterial mmHg)
x8 = 53.9946;   % p_Vc (venous mmHg)
x9 = 46.8969;   % p_Ac (arterial mmHg)
x10 = 0.0327;   % z_V (venous mmHg)
x11 = 0.0314;   % z_A (arterial mmHg)
x12=0;			% x_m  Muscle displacement

x0 = [x1;x2;x3;x4;x5;x6;x7;x8;x9;x10;x11;x12];

return