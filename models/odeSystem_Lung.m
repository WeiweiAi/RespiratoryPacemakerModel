% Calculates derivatives given x & t

function [dxdt] = odeSystem_Lung(t,x,parL,parRp,ventilationType,t0)

P_A = x(1);
f_o = x(2);
f_c = x(3);
f_do = x(4);
f_dc = x(5);
p_Vo = x(6);
p_Ao = x(7);
p_Vc = x(8);
p_Ac = x(9);
z_V = x(10);
z_A = x(11);
x_m=x(12);

[xAux] = calcAuxVars(t,x,parL,parRp,ventilationType,t0);

P_L = xAux(1);
dP_Ldt = xAux(2);
Q_A = xAux(3);
q = xAux(4);
p_ao = xAux(5);
p_ac = xAux(6);
V_A = xAux(7);
f_o_cos = xAux(8);
f_c_cos = xAux(9);
f_do_cos = xAux(10);
f_dc_cos = xAux(11);
dfpVodpVo = xAux(12);
dfpAodpAo = xAux(13);
S_Ao=xAux(14);
S_Vo=xAux(15);
R_p=xAux(16);


p_bo = p_Ao; % p_b = partial pressure in pulmonary capillaries which does the flux
p_bc = p_Ac;


dx_mdt=-parL.k1.*x_m+parL.k2.*R_p;

dP_Adt = (parL.P_m.*parL.E.*Q_A)./P_A + dP_Ldt;

df_odt = (1./V_A).*((f_o_cos-f_o).*q + parL.Do.*(p_bo - p_ao) - f_o.*(parL.Dc.*(p_bc - p_ac) + parL.Do.*(p_bo - p_ao)));
df_cdt = (1./V_A).*((f_c_cos-f_c).*q + parL.Dc.*(p_bc - p_ac) - f_c.*(parL.Do.*(p_bo - p_ao) + parL.Dc.*(p_bc - p_ac)));

df_dodt = (1./parL.V_d).*(f_do_cos - f_do).*(abs(q));
df_dcdt = (1./parL.V_d).*(f_dc_cos - f_dc).*(abs(q));

% with blood flux from haemoglobin O2 binding
dp_Aodt = (1./parL.V_p).*((parL.sigma_o + 4.*parL.Th.*dfpAodpAo).^(-1)).*(parL.Dbo.*(p_ao - p_Ao) + parL.Q_CO.*parL.sigma_o.*(p_Vo - p_Ao) + parL.Q_CO.*4.*parL.Th.*(S_Vo - S_Ao));
dp_Vodt = (1./parL.V_s).*((parL.sigma_o + 4.*parL.Th.*dfpVodpVo).^(-1)).*(parL.M_o + parL.Q_CO.*parL.sigma_o.*(p_Ao - p_Vo) + parL.Q_CO.*4.*parL.Th.*(S_Ao - S_Vo));


dp_Acdt = (parL.Q_CO./parL.V_p).*(p_Vc - p_Ac) + (parL.Dbc./(parL.V_p.*parL.sigma_c)).*(p_ac - p_Ac) + (parL.delta.*parL.l2.*z_A.*parL.h)./(parL.sigma_c) - (parL.delta.*parL.r2.*p_Ac);
dp_Vcdt = (parL.Q_CO./parL.V_s).*(p_Ac - p_Vc) + (parL.M_c./(parL.V_s.*parL.sigma_c)) + (parL.delta.*parL.l2.*z_V.*parL.h)./(parL.sigma_c) - (parL.delta.*parL.r2.*p_Vc);


% M_c is a function of P_ac
%M_c_adjusted = M_c - (sigma_t*p_Vc*BodyMass); % p_body_tissue ~ p_Vc 
%dp_Vcdt = (Q_CO./V_s).*(p_Ac - p_Vc) + (M_c_adjusted./(V_s.*sigma_c)) + (delta.*l2.*z_V.*h)./(sigma_c) - (delta.*r2.*p_Vc);


% with blood flux from bicarbonate
dz_Adt = (parL.delta.*parL.r2.*parL.sigma_c.*p_Ac) - (parL.delta.*parL.l2.*z_A.*parL.h) + (parL.Q_CO./parL.V_p).*(z_V - z_A);
dz_Vdt = (parL.delta.*parL.r2.*parL.sigma_c.*p_Vc) - (parL.delta.*parL.l2.*z_V.*parL.h) + (parL.Q_CO./parL.V_s).*(z_A - z_V);

dxdt = [dP_Adt; df_odt; df_cdt; df_dodt; df_dcdt; dp_Vodt; dp_Aodt; dp_Vcdt; dp_Acdt; dz_Vdt; dz_Adt;dx_mdt];

return