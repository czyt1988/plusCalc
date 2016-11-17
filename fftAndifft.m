massFlow = load(fullfile(currentPath,'mass_flow_0.1478_NorthZone.txt'));
time = massFlow(1:4096,1);
Fs = 1/(time(2)-time(1));
massFlowRaw = massFlow(1:4096,2);
[Fre,Mag,Ph,massFlowE] = fun_fft(detrend(massFlowRaw),Fs);
trData = changToWave(massFlowE,Fre,time);
figure
subplot(2,2,1)
plot(time,detrend(massFlowRaw));
title('质量流量');
subplot(2,2,2)
plot(Fre,Mag);
title('质量流量频谱');
subplot(2,2,[3,4])
hold on;
plot(time,detrend(massFlowRaw),'-b');
plot(time,detrend(abs(trData)),'-r');
title('对比');

set(gcf,'color','w');
