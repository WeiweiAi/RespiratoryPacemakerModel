% Calculate auxiliary variables given x & t

function [xAux] = calcAuxVars(t,x,parL,parRp,ventilationType,t0)

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

switch(ventilationType)
    case('intrinsic')
        taf  = parRp(1);
        tbf  = parRp(2)+parRp(1);
        a0 = parRp(3);
        a1 = parRp(4);
        c0 = parRp(5);
        tcf = parRp(6)+parRp(2)+parRp(1);
        if t<=taf+t0
            R_p=a0*(t-t0)./((t-t0)+a1);
        else
            if t>taf+t0 && t<=tbf+t0
                R_p= ( a0.*a1./(a1+taf).^2).*(t-t0)+(a0.*taf./(taf+a1)-( a0.*a1./(a1+taf).^2).*taf);
            else
                if t>tbf+t0 && t<=tcf+t0
                    R_p=(( a0.*a1./(a1+taf).^2)*tbf+(a0.*taf./(taf+a1)- ( a0.*a1./(a1+taf).^2).*taf))*exp(-c0*((t-t0)-tbf));
                else
                    R_p=0;
                end
            end
        end
        P_L=parL.P_m-parL.PL0-parL.kp.*x_m;
        dP_Ldt=parL.k1.*parL.kp.*x_m-parL.k2.*parL.kp.*R_p;
    case('pacing')
        InspirT_p=parRp(1);
        a_p=parRp(2);
        b_p=parRp(3);
        if t<=InspirT_p+t0
            R_p=a_p+(b_p-a_p)/InspirT_p*(t-t0);
        else
            R_p=0;
        end
        P_L=parL.P_m-parL.PL0-parL.kp.*x_m;
        dP_Ldt=parL.k1.*parL.kp.*x_m-parL.k2.*parL.kp.*R_p;
    case ('cos')
        P_L = parL.P_m - (parL.R.*parL.omega.*parL.V_T./2).*sin(parL.omega.*t) - parL.E.*(2.5 - (parL.V_T./2).*cos(parL.omega.*t));
        dP_Ldt = (-parL.V_T.*parL.omega./2).*(parL.R.*parL.omega.*cos(parL.omega.*t) + parL.E.*sin(parL.omega.*t));
        R_p=dP_Ldt./(-parL.kp)+parL.k1.*x_m;
end


p_bo = p_Ao; % p_b = partial pressure in pulmonary capillaries which does the flux
p_bc = p_Ac;

p_ao = f_o.*(P_A - parL.p_w);
p_ac = f_c.*(P_A - parL.p_w);

q = (parL.P_m - P_A)./parL.R;
Q_A = q + parL.Dc.*(p_bc - p_ac) + parL.Do.*(p_bo - p_ao);

V_A = (P_A - P_L)./parL.E + parL.V_0;

% f_o_cos & f_c_cos = transfer between alveoli and dead-space
% f_do_cos & f_dc_cos = transfer between dead-space and atmospheric gas

if(q > 0)     % inspiration (q)
    f_o_cos = f_do;
    f_c_cos = f_dc;
    f_do_cos = parL.f_mo;
    f_dc_cos = parL.f_mc;
elseif(q < 0) % expiration (q)
    f_o_cos = f_o;
    f_c_cos = f_c;
    f_do_cos = f_o;
    f_dc_cos = f_c;
else % q = 0 initially since P_A = P_m
    f_o_cos = f_o;
    f_c_cos = f_c;
    f_do_cos = f_o;
    f_dc_cos = f_c;
end

% O2-haemoglobin binding Monod-Wyman-Changeux
[dfpVodpVo] = calcSaturationGradient(p_Vo,parL);
[dfpAodpAo] = calcSaturationGradient(p_Ao,parL);

S_Ao = calcSaturation(p_Ao,parL);
S_Vo = calcSaturation(p_Vo,parL);

xAux = [P_L;dP_Ldt;Q_A;q;p_ao;p_ac;V_A;f_o_cos;f_c_cos;f_do_cos;f_dc_cos;dfpVodpVo;dfpAodpAo;S_Ao;S_Vo;R_p];
return