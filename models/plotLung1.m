%Plot and compare
q = xAux(:,4);
S_Ao=xAux(:,14);
R_p=xAux(:,16);
margins=[0.02 0.03]; % space among subplots
nplot=3; % the number of subplots
figure;
m=1;
ax(1,m)=subtightplot(nplot,1,m,margins);
plot(ax(1,m),t,q,tdata,flow);
xrange=xlim(ax(1,m));
m=2;
ax(1,m)=subtightplot(nplot,1,m,margins);
plot(ax(1,m),t,R_p);
xlim(ax(1,m),xrange);
m=3;
ax(1,m)=subtightplot(nplot,1,m,margins);
plot(ax(1,m),t,S_Ao.*100,tdata,SpO2);
xlim(ax(1,m),xrange);

%Plot and compare after SpO2 fitting
q = xAux(:,4);
S_Ao=xAux(:,14);
p_Ao=x(:,7);
R_p=xAux(:,16);
margins=[0.02 0.03]; % space among subplots
nplot=4; % the number of subplots
figure;
m=1;
ax(1,m)=subtightplot(nplot,1,m,margins);
plot(ax(1,m),t,q,itdata,iflow);
ylabel(ax(1,m),'q (L/s)')
xrange=xlim(ax(1,m));
grid on
m=2;
ax(1,m)=subtightplot(nplot,1,m,margins);
plot(ax(1,m),t,R_p);
ylabel(ax(1,m),'R_p')
xlim(ax(1,m),xrange);
grid on
m=3;
ax(1,m)=subtightplot(nplot,1,m,margins);
plot(ax(1,m),t,S_Ao.*100,itdata,iiSpO2,itdata,mSpO2);
ylabel(ax(1,m),'O_2 saturation')
xlim(ax(1,m),xrange);

m=4;
ax(1,m)=subtightplot(nplot,1,m,margins);
plot(ax(1,m),t,p_Ao);
ylabel(ax(1,m),'p_Ao')

grid on

path='H:\GitHub\Cardiorespiratory\data\ISRUC-Sleep\sub2_2\2\2\';
name = '2.rec';
[hdr, record] =edfread([path name]);
t=0:1/hdr.samples(14):hdr.records;
nt=length(t);
figure;plot(t,record(14,1:nt),t,record(18,1:nt)/100)
