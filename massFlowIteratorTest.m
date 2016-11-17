clc;
close all;
clear;
currentPath = fileparts(mfilename('fullpath'));
%% 压缩机基本参数
DCylinder=250/1000;%缸径m
dPipe=98/1000;%管内径m
crank=140/1000;%曲柄长
connectingRod=1.075;%连杆长度
k=1.4;%绝热指数
%% 计算参数
fs = 1/0.0035;
totalSecond = 5;
totalPoints = fs*totalSecond+1;
time = linspace(0,totalSecond,totalPoints);%0:1/fs:totalSecond;
periodIndex = time < (1/5);
periodIndex = find(periodIndex ==1);
%% 若认为质量流量为正弦情况
massFlowSin = 0.5.*abs(sin(10*pi*time));
[freSin,magSin] = frequencySpectrum(massFlowSin,fs);
figure
subplot(2,1,1)
plot(time(periodIndex),massFlowSin(periodIndex));
subplot(2,1,2)
plotSpectrum(freSin,magSin);
xlim([0,50]);
set(gcf,'color','w');
%% 可变参数


%% 改变转速对其影响
rcv = 0.4;% 相对余隙容积
pressureRadio=1.5;%压力比（排气压力/吸气压力）
outDensity =1.293*pressureRadio;% 1.9167;%%环境25度绝热压缩到0.2的温度对应密度
rpm = 200:10:500;
for i=1:length(rpm)
    [massFlowRPM{i},~,meanMassFlowRPM{i}] = massFlowMaker(DCylinder,dPipe,rpm(i)...
        ,crank,connectingRod,outDensity,'rcv',rcv,'k',k,'pr',pressureRadio,'fs',fs);
    while length(massFlowRPM{i})<totalPoints %生成1s的数据
        massFlowRPM{i} = [massFlowRPM{i},massFlowRPM{i}];
    end
    massFlowRPM{i} = massFlowRPM{i}(1:totalPoints);
end

figure
hold on;
for i=1:length(rpm)
    periodIndex = time < (1 / (rpm(i) / 60));
    periodIndex = find(periodIndex ==1);
    x = ones(1,length(periodIndex)).*rpm(i);
    y = time(periodIndex);
    z = massFlowRPM{i}(periodIndex);
   	plot3(x,y,z);
    xlabel('rpm');
    ylabel('time(s)');
    zlabel('mass flow(kg/s)');
end
grid on;
view(-68,56);
set(gcf,'color','w');

figure
hold on;
for i=1:length(rpm)
    bseFre = rpm(i)/60*2;
    [fre,amp] = frequencySpectrum(massFlowRPM{i},fs);
    ans = calcWaveFreAmplitude(massFlowRPM{i},fs,[bseFre,bseFre*2]);
    ampValue(:,i) = ans';
    index = fre < 50;
    index = find(index == 1);
    x = ones(1,length(index)).*rpm(i);
    y = fre(index);
    z = amp(index);
    
   	plot3(x,y,z);
    xlabel('rpm');
    ylabel('frequency(Hz)');
    zlabel('mass flow(kg/s)');
end
view(-117,44);
set(gcf,'color','w');

figure
hold on;
plot(rpm,ampValue(1,:),'r');
plot(rpm,ampValue(2,:),'b');
title(sprintf('转速和倍频间关系-余隙:%g,压比：%g',rcv,pressureRadio));
set(gcf,'color','w');


