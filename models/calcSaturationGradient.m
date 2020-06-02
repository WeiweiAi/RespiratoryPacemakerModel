% Calculate O2-haemoglobin binding curve gradient

% Check wolfram: differentiate y=((a*x*(1+b*x)^3) + (c*x*(1+ d*x)^3))/(e*(1+f*x)^4 + (1+ g*x)^4)

function [dfdpo] = calcSaturationGradient(po,parL)
% O2-haemoglobin binding
A = parL.L.*parL.KT.*parL.sigma_o;
B = parL.KT.*parL.sigma_o;
C = parL.KR.*parL.sigma_o;

G = A.*po.*((1 + B.*po).^3) + C.*po.*((1 + C.*po).^3);
H = parL.L.*((1 + B.*po).^4) + (1 + C.*po).^4;
Hsqrd = H.^2;
dGdpo = A.*((1 + B.*po).^2).*(1 + 4.*B.*po) + C.*((1 + C.*po).^2).*(1 + 4.*C.*po);
dHdpo = 4.*(parL.L.*B.*((1 + B.*po).^3) + C.*((1 + C.*po).^3));
dfdpo = ((dGdpo.*H) - (G.*dHdpo))./(Hsqrd);
return