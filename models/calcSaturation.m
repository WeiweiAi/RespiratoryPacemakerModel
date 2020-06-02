% Calculate O2-haemoglobin saturation function - Monod-Wyman-Changeux (MWC)

function [fpo] = calcSaturation(po,parL)

% A = L.*KT.*sigma_o.*po.*((1 + KT.*sigma_o.*po).^3);
% B = KR.*sigma_o.*po.*((1 + KR.*sigma_o.*po).^3);
% C = L.*(1 + KT.*sigma_o.*po).^4;
% D = (1 + KR.*sigma_o.*po).^4;
% fpo = (A + B)./(C + D);

A = parL.L.*parL.KT.*parL.sigma_o;
B = parL.KT.*parL.sigma_o;
C = parL.KR.*parL.sigma_o;

G = A.*po.*((1 + B.*po).^3) + C.*po.*((1 + C.*po).^3);
H = parL.L.*((1 + B.*po).^4) + (1 + C.*po).^4;
fpo = G./H;

return