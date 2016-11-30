clc;
close all;
clear;
currentPath = fileparts(mfilename('fullpath'));
%% 余隙容积对质量流量的影响

rcv = [0.05:0.01:0.35];% 相对余隙容积
pressureRadio=2;%压力比（排气压力/吸气压力）1.5-3.5
k=1.4;%绝热指数
DCylinder=250/1000;%缸径m
dPipe=98/1000;%管内径m
crank=140/1000;%曲柄长
connectingRod=1.075;%连杆长度
outDensity =1.293*pressureRadio;% 1.9167;%%环境25度绝热压缩到0.2的温度对应密度
fs = 1/0.0035;
totalSecond = 5;
rpm=300;%转速r/min
%%
totalPoints = fs*totalSecond+1;

for i=1:length(rcv)
    [massFlow{i},time{i},meanMassFlow{i}] = massFlowMaker(DCylinder,dPipe,rpm...
        ,crank,connectingRod,outDensity,'rcv',rcv(i),'k',k,'pr',pressureRadio,'fs',fs);
    while length(massFlow{i})<totalPoints %生成1s的数据
        massFlow{i} = [massFlow{i},massFlow{i}];
    end
    massFlow{i} = massFlow{i}(1:totalPoints);
    time{i} = linspace(0,totalSecond,totalPoints);%0:1/fs:totalSecond;
end

%%
figure
hold on;
for i=1:length(rcv)
    x = ones(1,length(massFlow{i})).*rcv(i);
    y = time{i};
    z = massFlow{i};
   	plot3(x,y,z);
    xlabel('rcv');
    ylabel('time(s)');
    zlabel('mass flow(kg/s)');
end
ylim([0,1]);
set(gcf,'color','w');

figure
hold on;
for i=1:length(rcv)
    baseFre = rpm/60*2;
    [fre,amp] = frequencySpectrum(massFlow{i},fs);
    ampValue(:,i) = calcWaveFreAmplitude(massFlow{i},fs,[baseFre,baseFre*2]);
    x = ones(1,length(fre)).*rcv(i);
    y = fre;
    z = amp;
   	plot3(x,y,z);
    xlabel('rcv');
    ylabel('frequency(Hz)');
    zlabel('mass flow(kg/s)');
end
ylim([0,100]);
set(gcf,'color','w');

figure
plot(rcv,ampValue(2,:),'-r');

%% 导出txt
% fid = fopen(fullfile(currentPath,sprintf('northZone_300rpm_0.1MPa-%.4g.txt',meanMassFlow_300)),'w');
% if fid>0
%     for i=1:length(massFlow_300)
%         fprintf(fid,'%g,%g\r\n',time_300(i),massFlow_300(i));
%     end
%     fclose(fid);
% end
% fid = fopen(fullfile(currentPath,sprintf('northZone_420rpm_0.05MPa-%.4g.txt',meanMassFlow_420)),'w');
% if fid>0
%     for i=1:length(massFlow_420)
%         fprintf(fid,'%g,%g\r\n',time_420(i),massFlow_420(i));
%     end
%     fclose(fid);
% end


