%%
clc;
close all;
clear;
currentPath = fileparts(mfilename('fullpath'));
%% 初始参数
%
[time,massFlow,Fre,massFlowE ] = getMassFlowData('N',4096,'isfindpeaks',1);
[massFlow.time,massFlow.massFlow,massFlow.Fre,massFlow.massFlowE] = getMassFlowData('N',4096,'isfindpeaks',1);
Fs = 1/(time(2)-time(1));
[los]=find(Fre>19 & Fre <21);
Fre(los) = Fre(los)./1.5;

temp = Fre<20 | (Fre>22&Fre<60);%Fre<80;
Fre = Fre(temp);

massFlowE = massFlowE(temp);

simValue = [9.17543970000000	10.3250122100000	10.1907085000000	10.6239057600000	10.4358613300000	9.80200586000000	9.64491675000000	9.80446094000000	9.65491846000000	9.01961939000000	7.78409082000000	5.75730152000000	2.91316907000000	1.46925610400000	2.39139527000000	3.15910156000000	3.73210498000000	4.09421813000000	4.19560510000000	4.19263013000000	4.09637684000000	3.97438465000000	3.74883410000000	3.37208447000000	2.76874720000000	1.93113226400000	0.929064820000000	4.43451025300000e-06];
XsimValue = [0:13,((0:13) + 13+2+0.115*2)];
%%

acousticVelocity = 335;%声速（m/s）
isDamping = 1;
coeffFriction = 0.03;
meanFlowVelocity = 16;
L1 = 13;%(m)
L2 = 13;
Lv = 2;
l = 0.115;%(m)缓冲罐的连接管长
Dv = 0.5;
sectionL1 = linspace(0,L1,14);
sectionL2 = linspace(0,L2,14);
Dpipe = 0.157;%管道直径（m）
X = [sectionL1,sectionL2 + 13 + Lv + 2 * l];

baseFrequency = 10;
multFreTimes = 8;
allowDeviation = 0.5;
%% 计算脉动压力

[pressure1,pressure2] = oneVesselPulsationCalc(massFlowE,Fre,time,L1,L2,Lv,l,Dpipe,Dv ...
	,sectionL1,sectionL2 ...
	,'a',acousticVelocity,'isDamping',isDamping,'friction',coeffFriction,'meanFlowVelocity',meanFlowVelocity);
pressure = [pressure1,pressure2];

dcpss = getDefaultCalcPulsSetStruct();
dcpss.calcSection = [0.4,0.5];
dcpss.fs = Fs;
dcpss.isHp = 0;
dcpss.f_pass = 7;%通过频率5Hz
dcpss.f_stop = 5;%截止频率3Hz
dcpss.rp = 0.1;%边带区衰减DB数设置
dcpss.rs = 30;%截止区衰减DB数设置

[plus,filterData] = calcPuls(pressure,dcpss);

figure
plot(X,plus./1000,'-r');
hold on;
plot(XsimValue,simValue,'-b');
set(gcf,'color','w');

% [plus,pressure] = fun_straightPipe(massFlowE1,L,Dpipe,opt,sectionL);
% plus = plus'./1000;
[multFreMag] = calcBaseFres(pressure,Fs,baseFrequency,multFreTimes,allowDeviation);

figure
plot(multFreMag(1,:)./1000);
hold on;
plot(multFreMag(2,:)./1000);
plot(multFreMag(3,:)./1000);
set(gcf,'color','w');


