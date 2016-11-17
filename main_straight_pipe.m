%%
clc;
close all;
clear;
currentPath = fileparts(mfilename('fullpath'));
%% 初始参数
%
[time,massFlow,Fre,massFlowE ] = getMassFlowData('N',4096,'isfindpeaks',1);
Fs = 1/(time(2)-time(1));
[los]=find(Fre>19 & Fre <21);
% Fre(los) = Fre(los)./1.5;

temp = Fre<60;%Fre<20 | (Fre>22&Fre<80);
Fre = Fre(temp);

massFlowE = massFlowE(temp);
% figure
% subplot(2,2,1)
% plot(detrend(massFlowRaw),'-b');
% title('原始数据');

% subplot(2,2,2)
% plot(FreRaw,AmpRaw,'-b');
% title('原始数据频谱');

% subplot(2,2,3)
% plot(detrend(massFlowRaw),'-r');
% hold on;
% plot(detrend(waveC),'-b');
% title('提取主频率数据波形');

% subplot(2,2,4)
% [FreTemp,AmpTemp] = fun_fft(detrend(waveC),Fs);
% plot(FreTemp,AmpTemp,'-r');

acousticVelocity = 335;%声速（m/s）
isDamping = 1;
coeffFriction = 0.04;
meanFlowVelocity = 14.5;
L=28.75;%L3(m)
Dpipe = 0.106;%管道直径（m）


%% 计算脉动压力
sectionL = linspace(0,28.75,30);
pressure = straightPipePulsationCalc(massFlowE,Fre,time,L,sectionL...
	,'d',Dpipe,'a',acousticVelocity,'isDamping',isDamping,'friction',coeffFriction,'meanFlowVelocity',meanFlowVelocity);
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
plot(plus./1000,'-r');
hold on;
%plot(realValue,'or');
set(gcf,'color','w');
% [plus,pressure] = fun_straightPipe(massFlowE1,L,Dpipe,opt,sectionL);
% plus = plus'./1000;
% [multFreMag] = calcBaseFres(pressure,Fs,baseFrequency,multFreTimes,allowDeviation);
% 
% figure
% plot(multFreMag(1,:)./1000);
% hold on;
% plot(multFreMag(2,:)./1000);
% plot(multFreMag(3,:)./1000);
% set(gcf,'color','w');


