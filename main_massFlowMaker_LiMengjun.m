clc;
close all;
clear;
currentPath = fileparts(mfilename('fullpath'));

rcv = 0.11;% 相对余隙容积
pressureRadio=1.5;%压力比（排气压力/吸气压力）
k=1.4;%绝热指数
DCylinder=326/1000;%缸径m
dPipe=98/1000;%管内径m
crank=140/1000;%曲柄长
connectingRod=1.075;%连杆长度
%50deg,0.15MPa(A),Density=1.617
outDensity = 1.617;%1.9167;%环境25度绝热压缩到0.2的温度对应密度
fs = 1/0.0035;
totalSecond = 5;
%%
totalPoints = fs*totalSecond+1;
rpm=300;%转速r/min
meanMassFlow_300 = 0;
[massFlow_300,time_300,meanMassFlow_300,meanFlowVelocity] = massFlowMaker(DCylinder,dPipe,rpm...
	,crank,connectingRod,outDensity,'rcv',rcv,'k',k,'pr',pressureRadio,'fs',fs);
disp(sprintf('300rpm:%g',meanMassFlow_300));
while length(massFlow_300)<totalPoints %生成1s的数据
    massFlow_300 = [massFlow_300,massFlow_300];
end
massFlow_300 = massFlow_300(1:totalPoints);
time_300 = linspace(0,totalSecond,totalPoints);%0:1/fs:totalSecond;

meanMassFlow_420 = 0;
rpm=420;%转速r/min
pressureRadio=0.15/0.1;%压力比（排气压力/吸气压力）
outDensity = 1.5608;
[massFlow_420,time_420,meanMassFlow_420,meanFlowVelocity] = massFlowMaker(DCylinder,dPipe,rpm...
	,crank,connectingRod,outDensity,'rcv',rcv,'k',k,'pr',pressureRadio,'fs',fs);
disp(sprintf('420rpm:%g',meanMassFlow_420));
while length(massFlow_420)<totalPoints
    massFlow_420 = [massFlow_420,massFlow_420];
end
massFlow_420 = massFlow_420(1:totalPoints);
time_420 = linspace(0,totalSecond,totalPoints);

mm = 0.5*abs(sin(pi*10*time_300));
figure
[fre,amp] = frequencySpectrum(mm,fs);
plotSpectrum(fre,amp,'isFill',1);
xlim([0,50]);
%%
[ time,massFlow,Fre,massFlowE ] = getMassFlowData();
figure
hold on;
plot(time_300,massFlow_300,'-r');
plot([0,time_300(end)],[meanMassFlow_300,meanMassFlow_300],'--r');
plot(time_420,massFlow_420,'-b');
plot([0,time_420(end)],[meanMassFlow_420,meanMassFlow_420],'--b');
plot(massFlow(1:3600,1),massFlow(1:3600,2),'--b');
set(gcf,'color','w');

figure
subplot(2,1,1)
plot(time_300,massFlow_300);
xlabel('time(s)');
ylabel('mass flow(kg/s)');
[fre,amp] = frequencySpectrum(massFlow_300,fs);
subplot(2,1,2)
plotSpectrum(fre,amp,'isFill',1);
xlim([0,100]);
set(gcf,'color','w');

figure
subplot(2,1,1)
plot(time_420,massFlow_420);
xlabel('time(s)');
ylabel('mass flow(kg/s)');
[fre,amp] = frequencySpectrum(massFlow_420,fs);
subplot(2,1,2)
plotSpectrum(fre,amp,'isFill',1);
xlim([0,100]);
set(gcf,'color','w');

% 导出txt
fid = fopen(fullfile(currentPath,sprintf('northZone_300rpm_0.1MPa-%.4g.txt',meanMassFlow_300)),'w');
if fid>0
    for i=1:length(massFlow_300)
        fprintf(fid,'%g,%g\r\n',time_300(i),massFlow_300(i));
    end
    fclose(fid);
end
fid = fopen(fullfile(currentPath,sprintf('northZone_420rpm_0.05MPa-%.4g.txt',meanMassFlow_420)),'w');
if fid>0
    for i=1:length(massFlow_420)
        fprintf(fid,'%g,%g\r\n',time_420(i),massFlow_420(i));
    end
    fclose(fid);
end


