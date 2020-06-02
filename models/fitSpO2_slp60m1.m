% O2 Saturation Fitting Code
% Weiwei Ai
% Date: 26/2/2020
clear;
% Get the filename
path=['..' filesep 'patients' filesep];
fname=[path 'slp60m1_2.mat'];
load(fname);
% the segment to fit
istart=1; % Index of the first segment;
iendt=191;% Index of the last segment;
istarttdata=tdata(DSense(istart));
iendtdata=tdata(DSense(iendt));
itdata=[tdata(DSense(istart):DSense(iendt))]; % The new time range; may be different from tdata
iitdata=[tdata(DSense(istart):DSense(iendt))]-11.5; % The time span shift forward 9.5s to deal with the SpO2 delay
tshift=find(iitdata>=tdata(1),1,'first');
iSpO2=SpO2(DSense(istart):DSense(iendt)); % The original SpO2 segment, i.e., the target to fit
iiSpO2=SpO2(DSense(istart)+tshift:DSense(iendt)+tshift);% similar shift 
iflow=flow(DSense(istart):DSense(iendt)); % The original flow segment,
%default Lung parameters
parL= savePars; % Get the parameters for the lung model;
x0=initialVariables; % Initialization
ventilationType='intrinsic'; %'intrinsic','pacing','cos'
% Initial parameter guesses
params0(1,1) = 0.28;% the M_o 
params0(1,2) = 3.9724;% the Q_CO
params0(1,3) = 25;%x0(6)p_Vo
params0(1,4) = 60;%x0(7)p_Ao
lb(1,1) = 0.128; % MO 0.128-0.349 L/min
lb(1,2) = 2.96; % CO 2.96-5.92 L/min
lb(1,3) = 25; %x0(7)p_Ao 50-100
lb(1,4) = 50; %x0(7)p_Ao 50-100
ub(1,1) = 0.349;
ub(1,2) = 5.92;
ub(1,3) = 40;
ub(1,4) = 100;

%x0(7)=60;
optionsfit = optimoptions('lsqcurvefit','OptimalityTolerance',1e-9,'StepTolerance',1e-9,'MaxFunctionEvaluations',2000); % for lsqcurvefit
options=odeset('RelTol',1e-6,'AbsTol',1e-9);

f = @(x,targetx)getLungModelSpO2(x,targetx,breaths,istart,iendt,x0,ibeta,parL,ventilationType,options);
ialpha = lsqcurvefit(f,params0,itdata,iiSpO2,lb,ub,optionsfit);
%ialpha=[0.28 3.9724 25 60];
%% fit the patient
parL.M_o=parL.M_o*ialpha(1,1)/0.25;
parL.M_c=parL.M_o*(-0.8);
parL.Q_CO=ialpha(1,2)/60;
x0(6)=ialpha(1,3);
x0(7)=ialpha(1,4);

[t,x,xAux,mSpO2]= combFittingSpO2(itdata,breaths,istart,iendt,x0,ibeta,parL,ventilationType,options); % get the fitted simulation output

q = xAux(:,4);
vq = interp1(t,q,itdata,'linear','extrap'); % interpolate q at tdata points
RMSEq =sqrt(mean((vq - iflow).^2));
RMSESpO2=sqrt(mean((mSpO2 - iiSpO2).^2));

save([path 'slp60m1_3.mat'],'t','x','xAux','parL','x0',...
    'ibeta','T','DSense','breaths','SpO2','tdata','flow',...
    'ialpha','itdata','iSpO2','iiSpO2','mSpO2','iflow','RMSEq','RMSESpO2');

function mSpO2=getLungModelSpO2(alpha0,itdata,breaths,istart,iendt,x0,ibeta,parL,ventilationType,options)
% simulate the lung model cycle by cycle
t=[];
x=[];
xAux=[];
% Update the parameter
parL.M_o=parL.M_o*alpha0(1,1)/0.25;
parL.M_c=parL.M_o*(-0.8);
parL.Q_CO=alpha0(1,2)/60;
x0(6)=alpha0(1,3);
x0(7)=alpha0(1,4);
for i=istart:iendt
    targetx  = breaths{i}(2,:); %time
    t0 = targetx(1);
    tf = targetx(end);
    tspan=[t0 tf];
    parRp=ibeta(i,:)';
    [ti, xi] = ode23t(@odeSystem_Lung, tspan, x0, options,parL,parRp,ventilationType,t0 );
    x0=xi(end,:); % The new initial state should be the last value of previous cycle simulation.
    xAuxi=[];% To remove the old values
    for tn=1:length(ti)
        xAuxi(tn,:) = calcAuxVars(ti(tn),xi(tn,:),parL,parRp,ventilationType,t0);
    end
    if i==istart
        t=ti;
        x=xi;
        xAux=[xAuxi];
    else % Remove the first repeated simulation instance
        t=[t;ti(2:end)];
        x=[x;xi(2:end,:)];
        xAux=[xAux;xAuxi(2:end,:)];
    end
end
S_Ao=xAux(:,14).*100;
vSpO2 = interp1(t,S_Ao,itdata,'linear','extrap'); % interpolate S_Ao at tdata points
[yupper,ylower] = envelope(vSpO2,150,'peak');
mSpO2 =(yupper+ylower)./2; % get the mean value
end

function [t,x,xAux,mSpO2]=combFittingSpO2(itdata,breaths,istart,iendt,x0,ibeta,parL,ventilationType,options)
% simulate the lung model cycle by cycle
t=[];
x=[];
xAux=[];
for i=istart:iendt
    targetx  = breaths{i}(2,:); %time
    t0 = targetx(1);
    tf = targetx(end);
    tspan=[t0 tf];
    parRp=ibeta(i,:)';
    [ti, xi] = ode23t(@odeSystem_Lung, tspan, x0, options,parL,parRp,ventilationType,t0 );
    x0=xi(end,:); % The new initial state should be the last value of previous cycle simulation.
    xAuxi=[];% To remove the old values
    for tn=1:length(ti)
        xAuxi(tn,:) = calcAuxVars(ti(tn),xi(tn,:),parL,parRp,ventilationType,t0);
    end
    if i==istart
        t=ti;
        x=xi;
        xAux=[xAuxi];
    else % Remove the first repeated simulation instance
        t=[t;ti(2:end)];
        x=[x;xi(2:end,:)];
        xAux=[xAux;xAuxi(2:end,:)];
    end
end
S_Ao=xAux(:,14).*100;
vSpO2 = interp1(t,S_Ao,itdata,'linear','extrap'); % interpolate S_Ao at tdata points
[yupper,ylower] = envelope(vSpO2,150,'peak');
mSpO2 =(yupper+ylower)./2; % get the mean value
end

