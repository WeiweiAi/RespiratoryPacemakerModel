% Phrenic Nerve Signal Fitting Code
% Author: Chad Eichler
% Date: 16/2/18
% Modified by Weiwei Ai
clear;
path=['..' filesep 'patients' filesep];
fname=[path 'slp60m1'];
[val,tdata,units,signal] = plotATM_mod(fname);
tstart=2800;%2800
tend=3600; %3600
nflow=5; % the index of air flow
flowshift=0.2686; % baseline shift 0.2686;
flowscale=1; % scale the q
nSpO2=7; % the index of SpO2 %
[breaths,T,tdata,flow,SpO2,DSense]=getPatientData(val,tdata,tstart,tend,nflow,flowshift,flowscale,nSpO2);
% Get the parameters for the lung model;
parL= savePars;
ventilationType='intrinsic'; %'intrinsic','pacing','cos'
parRp=saveparRp(ventilationType); % Get the parameters for the phrenic nerve model;
optionsfit = optimoptions('lsqcurvefit','OptimalityTolerance',1e-12,'FunctionTolerance',1e-12,'StepTolerance',1e-12,...
    'MaxFunctionEvaluations',1000); % for lsqcurvefit
x0=initialVariables; % Initialization
options=odeset('RelTol',1e-6,'AbsTol',1e-9); % for ode
n=length(breaths); %
istart=1;
iend=n;

ibeta=fitFlow(breaths,istart,iend,x0,parL,ventilationType,options,optionsfit); % get the fitted Rp parameters
[t,x,xAux]= combFitting(breaths,istart,iend,x0,ibeta,parL,ventilationType,options); % get the fitted simulation output

save([path 'slp60m1_2.mat'],'t','x','xAux','ibeta','T','DSense','breaths','n','SpO2','tdata','flow');

% fit the flow cycle by cycle
function ibeta=fitFlow(breaths,istart,iend,x0,parL,ventilationType,options,optionsfit)
% Initial parameter guesses
n=iend;
params0 = ones(n+1,6);
params0(:,1) = params0(:,1).*0.8;
params0(:,2) = params0(:,2).*0.8;
params0(:,3) = params0(:,3).*1;
params0(:,4) = params0(:,4).*0.5;
params0(:,5) = params0(:,5).*3;
params0(:,6) = params0(:,6).*1;%1
ibeta = params0;
lb = ones(n+1,6);
lb(:,1) = lb(:,1).*0.001;
lb(:,2) = lb(:,2).*0.001;
lb(:,3) = lb(:,3).*0.001;
lb(:,4) = lb(:,4).*0.001;
lb(:,5) = lb(:,5).*0.001;
lb(:,6) = lb(:,6).*0.001;
ub = ones(n+1,6);
ub(:,1) = ub(:,1).*3;
ub(:,2) = ub(:,2).*6;
ub(:,3) = ub(:,3).*3; % ub(:,3).*+inf;
ub(:,4) = ub(:,4).*10; %inf
ub(:,5) = ub(:,5).*10; %inf
ub(:,6) = ub(:,6).*2; %inf
for i=istart:iend
    targetx  = breaths{i}(2,:); %time
    targety = breaths{i}(1,:); % flow
    f = @(x,targetx)getLungModelFlow(x,targetx,x0,parL,ventilationType,options);
    [ibeta(i,:),resnorm(i),residual{i},exitflag(i),output{i}]= lsqcurvefit(f,params0(i,:),targetx,targety,lb(i,:),ub(i,:),optionsfit);
    RMSE(i) =sqrt(resnorm(i)/length(targety));
    if RMSE(i)>=0.05 % if the 1st fit error >=0.05, try another round (tend to decrease)
        resnorm0=resnorm(i);
        params1(3)=0.5;
        params1(1:2)=params0(i,1:2);
        params1(4:6)=params0(i,4:6);
        [ibetatemp,resnormtemp,residualtemp,exitflagtemp,outputtemp]= lsqcurvefit(f,params1,targetx,targety,lb(i,:),ub(i,:),optionsfit);
        RMSEtemp =sqrt(resnormtemp/length(targety));
        if RMSEtemp< RMSE(i) % if the 2nd fit is better than the 1st one; the best is the 2nd
            RMSE(i)= RMSEtemp;
            ibeta(i,:)=ibetatemp;
        end
        if RMSEtemp>=0.05 % if the 2nd fit is still >=0.05, try again (tend to decrease)
            params1(1:3)=params0(i,1:3)*0.5;
            params1(4:6)=params0(i,4:6)*1;
            [ibetatemp2,resnormtemp2,residualtemp2,exitflagtemp2,outputtemp2]= lsqcurvefit(f,params1,targetx,targety,lb(i,:),ub(i,:),optionsfit);
            RMSEtemp2 =sqrt(resnormtemp2/length(targety));
            if RMSEtemp2<RMSE(i) % if the 3rd fit is better than the best (2nd or 1st), the best is the 3rd
                RMSE(i)= RMSEtemp2;
                ibeta(i,:)=ibetatemp2;
            end
          if RMSEtemp2>=0.05 % if the 3rd fit is still >=0.05, try again (tend to increase)
            params1(1:3)=params0(i,1:3)*1.5;
            params1(4:6)=params0(i,4:6)*1.5;
            [ibetatemp3,resnormtemp3,residualtemp3,exitflagtemp3,outputtemp3]= lsqcurvefit(f,params1,targetx,targety,lb(i,:),ub(i,:),optionsfit);
            RMSEtemp3 =sqrt(resnormtemp3/length(targety));
            if RMSEtemp3<RMSE(i) % if the 4th fit is better than the best (3rd,2nd or 1st), the best is the 4th
                RMSE(i)= RMSEtemp3;
                ibeta(i,:)=ibetatemp3;
            end 
          end
        end
    end
end
end
%%
function vq=getLungModelFlow(parRp,tdata,x0,parL,ventilationType,options)
t0 = tdata(1);
tf = tdata(end);
tspan=[t0 tf];
[ti, xi] = ode23t(@odeSystem_Lung, tspan, x0, options,parL,parRp,ventilationType,t0 );
for tn=1:length(ti)
    xAuxi(tn,:) = calcAuxVars(ti(tn),xi(tn,:),parL,parRp,ventilationType,t0);
end
q = xAuxi(:,4);
vq = interp1(ti,q,tdata,'linear','extrap'); % interpolate q at tdata points
end

%% get patients ventilation patterns
function [breaths,T,tdata,flow,SpO2,DSense]=getPatientData(val,tdata,tstart,tend,nflow,flowshift,flowscale,nSpO2)
% Extract data
tori = tdata; % original time
% Specify the period of the signal
tstartindex=find(tdata>=tstart,1,'first');%
tendindex=find(tdata>=tend,1,'first');%
% Downsample every nth
n=10;
flowo = downsample(val(nflow,tstartindex:tendindex)./flowscale,n)+flowshift; % ucddb ./60 flowscale=60;'slp59m' flowshift-0.181; slp67xm flowshift-0.066
% smooth the data
windowSize = 7;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
flow=filter(b,a,flowo);
SpO2=downsample(val(nSpO2,tstartindex:tendindex),n);
tdata = downsample(tdata(tstartindex:tendindex),n);% may not start from 0;
% Find the start of inspiration
tempq=flow; % q
qthreshold=0.05;
tempq(tempq>=qthreshold)=1; % when q>0, q'=1;
tempq(tempq<qthreshold)=0; % when q<0, q'=0;
indexqin = find( diff(tempq)==1)-5; %start of inspiration, back 5 samples*0.04
if indexqin(1)<0
    indexqin(1)=1;
end
% Find volume
%     vol = cumtrapz(tdata,flow);
%     % Find local minima
%     [pks2,indexqin] = findpeaks(-vol,'MinPeakDistance',20,'MinPeakProminence',0.03);
%% remove very shallow breathing
nnq=0;
newindexqin=[];
for nq=1:length(indexqin)-1
    flowi=flow(indexqin(nq):indexqin(nq+1));
    if max(flowi)>=qthreshold*1.1 %0.02
        nnq=nnq+1;
        newindexqin(nnq)=indexqin(nq);
    end
end
newindexqin(nnq+1)=indexqin(nq+1);
DSense = newindexqin; % index of the beginning of inspiration => diaphragm event sensed => DSense
%     % Find volume
%     vol = cumtrapz(tdata,flow);
%     % Find local maxima
%     [pks,locs] = findpeaks(vol,'MinPeakDistance',0.5,'MinPeakProminence',0.03);
%     % Find local minima
%     [pks2,locs2] = findpeaks(-vol,'MinPeakDistance',0.5,'MinPeakProminence',0.03);
%    DSense = locs2; % local minima => beginning of inspiration => diaphragm event sensed => DSense
% Period Calculation
T = [tdata(DSense(1)),tdata(DSense(2:length(DSense)))-tdata(DSense(1:(length(DSense)-1)))];
% Define number of breaths to fit
n = round(length(T));
% Create breath airflow cell array
breaths = {};
for i = 1:n-1
    breaths{i} = [flow(DSense(i):(DSense(i+1)));tdata(DSense(i):(DSense(i+1)))];
end
end

function [t,x,xAux]= combFitting(breaths,istart,iend,x0,ibeta,parL,ventilationType,options)
t=[];
x=[];
xAux=[];
for i=istart:iend
    targetx  = breaths{i}(2,:); %time
    t0 = targetx(1);
    tf = targetx(end);
    tspan=[t0 tf];
    parRp=ibeta(i,:);
    [ti, xi] = ode23t(@odeSystem_Lung, tspan, x0, options,parL,parRp,ventilationType,t0 );
    x0=xi(end,:); % The new initial state should be the last value of previous simulation.
    xAuxi=[];% To remove the old values
    for tn=1:length(ti)
        xAuxi(tn,:) = calcAuxVars(ti(tn),xi(tn,:),parL,parRp,ventilationType,t0);
    end
    if i==1
        t=ti;
        x=xi;
        xAux=[xAuxi];
    else % Remove the first repeated simulation instance
        t=[t;ti(2:end)];
        x=[x;xi(2:end,:)];
        xAux=[xAux;xAuxi(2:end,:)];
    end
end
end
