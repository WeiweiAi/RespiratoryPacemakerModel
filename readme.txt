*********************************
1. Fit patient slp60m
slp60m1.mat, slp60m1.info, the original data.
(1) 
Input: slp60m1.mat 
run fitLung_slp60m1.m, fit q with multiple starting points(4), q threshold 0.05;
output: slp60m1_2.mat
(2)
Input: slp60m1_2.mat
run fitSpO2_slp60m1.m: 

%ialpha=[0.28 3.9724 25 60];
%% fit the patient
parL.M_o=parL.M_o*ialpha(1,1)/0.25;
parL.M_c=parL.M_o*(-0.8);
parL.Q_CO=ialpha(1,2)/60;
x0(6)=ialpha(1,3);
x0(7)=ialpha(1,4);
 to fit SpO2 tshift=-11.5;
 
output: slp60m1_3.mat
slp60m1_3_fixed.mat, after fixing the lung model.

load('H:\GitHub\RespiratoryPacemakerModel\patients\slp60m1_3_fixed.mat')
q = xAux(:,4);
S_Ao=xAux(:,14);
p_Ao=x(:,7);
R_p=xAux(:,16);
rec=[itdata',iflow',iiSpO2',mSpO2'];
sim=[t,q,R_p,S_Ao.*100];
sim=sim(1:5:end,:);
csvwrite('sim0_fixed.csv',sim);
csvwrite('rec0_fixed.csv',rec);

**************************************
2. Fit patient slp67xm
slp67xm.mat, slp67xm.info, the original data.
slp67xm_4.mat: fit q with multiple starting points (3),using fitLung_slp67xm.m.
slp67xm_6: using  slp67xm_4.mat, fitSpO2_slp67xm.
ialpha=[0.128 2.96];
parL.M_o=parL.M_o*ialpha(1,1)/0.25;
parL.M_c=parL.M_o*(-0.8);
parL.Q_CO=ialpha(1,2)/60;
x0(7)=60;
 to fit SpO2
********************************************
3.Parameters for closed-loop validation:
parL.mat: Kp=0.01; no xAux0; fitted, the min, MO 0.128L/min CO 2.96 L/min
			MC=6.70592917856088E-5, Mo=-8.382411473201099E-5, Co=0.04933333333333333
parL0.mat: Kp=0.04; MC=6.70592917856088E-5, Mo=-8.382411473201099E-5, Co=0.04933333333333333
parL1.mat: Kp=0.04;typical MO,MC, CO (MO 0.25 L/min CO 5 L/min), MC=1.3097517926876717E-4,Mo=-1.6371897408595895E-4,Co=0.0833
parL2.mat: Ki=0.002;Kp=0.4;LRI=6.5;SAo_set=96;fitted, the min, MO 0.128L/min CO 2.96 L/min
parL3.mat: Ki=0.002;Kp=0.4;LRI=6.5;SAo_set=96;typical MO,MC, CO (MO 0.25 L/min CO 5 L/min), MC=1.3097517926876717E-4,Mo=-1.6371897408595895E-4,Co=0.0833
*********************************************************************************************************************************************
4.Closed-loop validation
closedLoop.slx 
simX0, fitted validation, parL.mat
simX1, use typical MO,MC, CO, parL1.mat
simX2, fitted validation, parL.mat, fixed pacing
simX3, use typical MO,MC, CO, parL1.mat,fixed pacing
simX4, use typical MO,MC, CO, parL1.mat,fixed pacing, LRI=5.858609467017760, Kp=0.014878141945917

closedLoop_SpO2.slx
simX5, parL2.mat, fixed pacing
simX6, parL2.mat, adaptive pacing
simX7, parL3.mat, fixed pacing
simX8, parL3.mat, adaptive pacing
*********************************************************
5.Devices parametrization:
closedLoop_SpO2L.slx, closedLoop_SpO2H.slx;
runClosedLoopSo2_30.m, runClosedLoopSo2_30_2.m, runClosedLoopSo2_30_3.m, runClosedLoopSo2_60.m, runClosedLoopSo2_60_2.m
****************************************************************************************************************************
6. Data for the figures in the paper
Fig. 2(a) rec0.csv, sim0.csv;
Fig. 2(b) sim0.csv;
Fig. 2(c) rec.csv, sim.csv;
Fig. 2(b) sim.csv;
Fig. 3 dataHDp1.csv;
Fig. 4(a) dataLD0.csv, dataLD1.csv;
Fig. 4(b) dataHD0.csv, dataHD1.csv;
Fig. 5(b) LD30.csv;
Fig. 5(c) SpO2L30.csv, SpO2H30.csv;