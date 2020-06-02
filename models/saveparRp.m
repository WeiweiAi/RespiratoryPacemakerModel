function [parRp] = saveparRp(ventilationType)
switch(ventilationType)
    case('intrinsic')
        parRp(1) = 0.8; %taf
        parRp(2) =0.8; %tbf-taf
        parRp(3) =1; % a0
        parRp(4) =0.5;% a1
        parRp(5) =3; %c0
        parRp(6) =1; %tcf-tbf-taf
    case('pacing')
        parRp(1)=2.5; %InspirT_p, the inspiration duration for pacing, second
        parRp(2)=0.5; % Phrenic nerve jump,a_p
        parRp(3)=1.1; % The max phrenic activity, b_p
        case('cos')
        parRp=[];       
end
return