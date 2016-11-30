%%
clc;
close all;
clear;
currentPath = fileparts(mfilename('fullpath'));
%% 初始参数
%
dynViscosity = 0.02 * 10^-3;%动力粘度 60度
rpm = 420;outDensity = 1.5608;multFre=[14,28,42];%环境25度绝热压缩到0.15MPaG的温度对应密度
Fs = 4096;
[massFlowRaw,time,~,opt.meanFlowVelocity] = massFlowMaker(0.25,0.098,rpm...
	,0.14,1.075,outDensity,'rcv',0.15,'k',1.4,'pr',0.15,'fs',Fs,'oneSecond',6);
[FreRaw,AmpRaw,PhRaw,massFlowERaw] = frequencySpectrum(detrend(massFlowRaw,'constant'),Fs);
FreRaw = [7,14,21,28,14*3];
massFlowERaw = [0.02,0.2,0.03,0.05,0.007];
% 提取主要频率
massFlowE = massFlowERaw;
Fre = FreRaw;

acousticVelocity = 345;%声速（m/s）
isDamping = 1;
coeffFriction = 0.02;
meanFlowVelocity = 25.5;
L=10.62;%L3(m)
Dpipe = 0.098;%管道直径（m）
isOpening = 0;
dcpss = getDefaultCalcPulsSetStruct();
dcpss.calcSection = [0.25,0.75];
dcpss.sigma = 2.8;
dcpss.fs = Fs;
dcpss.isHp = 0;
dcpss.f_pass = 7;%通过频率5Hz
dcpss.f_stop = 5;%截止频率3Hz
dcpss.rp = 0.1;%边带区衰减DB数设置
dcpss.rs = 30;%截止区衰减DB数设置

mach = meanFlowVelocity / acousticVelocity;
calcWay2 = 0;

dataCount = 2;
calcDatas{1,2} = 'x值';
calcDatas{1,3} = '压力脉动';
calcDatas{1,4} = '1倍频';
calcDatas{1,5} = '2倍频';
calcDatas{1,6} = '3倍频';
%% 计算脉动压力
notmach = 1;
isDamping = 0;
sectionL = [[2.5,3.5],[4.29,5.08],[5.87,6.37,6.87,7.37,7.87,8.37,8.87,9.37,9.87,10.37]];

pressure = straightPipePulsationCalc(massFlowE,Fre,time,L,sectionL...
	,'d',Dpipe,'a'...
    ,acousticVelocity,'isDamping',isDamping,'friction',coeffFriction,'meanFlowVelocity',meanFlowVelocity ...
    ,'notmach',notmach,'m',mach,'isOpening',isOpening...
    );
[plus,filterData] = calcPuls(pressure,dcpss);
multFreAmpValue_straightPipe = calcWaveFreAmplitude(pressure,Fs,multFre,'freErr',1);
calcDatas{dataCount,1} = sprintf('直管-无阻尼-无马赫数');
calcDatas{dataCount,2} = sectionL;
calcDatas{dataCount,3} = plus;
calcDatas{dataCount,4} = multFreAmpValue_straightPipe(1,:);
calcDatas{dataCount,5} = multFreAmpValue_straightPipe(2,:);
calcDatas{dataCount,6} = multFreAmpValue_straightPipe(3,:);
dataCount = dataCount + 1;

%%
notmach = 1;
isDamping = 1;
pressure = straightPipePulsationCalc(massFlowE,Fre,time,L,sectionL...
	,'d',Dpipe,'a'...
    ,acousticVelocity,'isDamping',isDamping,'friction',coeffFriction,'meanFlowVelocity',meanFlowVelocity ...
    ,'notmach',notmach,'m',mach,'isOpening',isOpening...
    ,'calcWay2',calcWay2,'density',outDensity,'dynViscosity',dynViscosity...
    );
[plus,filterData] = calcPuls(pressure,dcpss);
multFreAmpValue_straightPipe = calcWaveFreAmplitude(pressure,Fs,multFre,'freErr',1);
calcDatas{dataCount,1} = sprintf('直管-有阻尼-无马赫数');
calcDatas{dataCount,2} = sectionL;
calcDatas{dataCount,3} = plus;
calcDatas{dataCount,4} = multFreAmpValue_straightPipe(1,:);
calcDatas{dataCount,5} = multFreAmpValue_straightPipe(2,:);
calcDatas{dataCount,6} = multFreAmpValue_straightPipe(3,:);
dataCount = dataCount + 1;
%%
notmach = 0;
isDamping = 1;
pressure = straightPipePulsationCalc(massFlowE,Fre,time,L,sectionL...
	,'d',Dpipe,'a'...
    ,acousticVelocity,'isDamping',isDamping,'friction',coeffFriction,'meanFlowVelocity',meanFlowVelocity ...
    ,'notmach',notmach,'m',mach,'isOpening',isOpening...
    ,'calcWay2',calcWay2,'density',outDensity,'dynViscosity',dynViscosity...
    );
[plus,filterData] = calcPuls(pressure,dcpss);
multFreAmpValue_straightPipe = calcWaveFreAmplitude(pressure,Fs,multFre,'freErr',1);
calcDatas{dataCount,1} = sprintf('直管-有阻尼-有马赫数');
calcDatas{dataCount,2} = sectionL;
calcDatas{dataCount,3} = plus;
calcDatas{dataCount,4} = multFreAmpValue_straightPipe(1,:);
calcDatas{dataCount,5} = multFreAmpValue_straightPipe(2,:);
calcDatas{dataCount,6} = multFreAmpValue_straightPipe(3,:);
dataCount = dataCount + 1;
%%
massFlowE(4) = 0.1;
notmach = 0;
isDamping = 1;
pressure = straightPipePulsationCalc(massFlowE,Fre,time,L,sectionL...
	,'d',Dpipe,'a'...
    ,acousticVelocity,'isDamping',isDamping,'friction',coeffFriction,'meanFlowVelocity',meanFlowVelocity ...
    ,'notmach',notmach,'m',mach,'isOpening',isOpening...
    ,'calcWay2',calcWay2,'density',outDensity,'dynViscosity',dynViscosity...
    );
[plus,filterData] = calcPuls(pressure,dcpss);
multFreAmpValue_straightPipe = calcWaveFreAmplitude(pressure,Fs,multFre,'freErr',1);
calcDatas{dataCount,1} = sprintf('直管-有阻尼-有马赫数-二倍频提高一倍');
calcDatas{dataCount,2} = sectionL;
calcDatas{dataCount,3} = plus;
calcDatas{dataCount,4} = multFreAmpValue_straightPipe(1,:);
calcDatas{dataCount,5} = multFreAmpValue_straightPipe(2,:);
calcDatas{dataCount,6} = multFreAmpValue_straightPipe(3,:);
dataCount = dataCount + 1;
%%

ignoreHeader = 1;
%绘制压力脉动
figure 
plotDataCells(calcDatas,'xcol',2,'ycol',3,'legendcol',1,'ignoreHeader',ignoreHeader);
title('脉动压力峰峰值');
%绘制1倍频
figure
plotDataCells(calcDatas,'xcol',2,'ycol',4,'legendcol',1,'ignoreHeader',ignoreHeader);
title('压力1倍频');
%绘制2倍频
figure
plotDataCells(calcDatas,'xcol',2,'ycol',5,'legendcol',1,'ignoreHeader',ignoreHeader);
title('压力2倍频');
%绘制3倍频
% figure
% plotDataCells(calcDatas,'xcol',2,'ycol',6,'legendcol',1,'ignoreHeader',ignoreHeader);
% title('压力3倍频');

