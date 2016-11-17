p = dataStruct.saMainFreFilterStruct.pressure(:,1);
p = p(1000:2500);
p(p>7.392 | p<1.361) = [];
f = dataStruct.rawData.Fre(:,1);
amp = dataStruct.rawData.Mag(:,1);
figure
subplot(2,1,1)
plot(0:(1/100):(1/100*(length(p)-1)),p);
xlabel('time (s)');
ylabel('pressure (kPa)');
subplot(2,1,2)
plotSpectrum(f,amp);