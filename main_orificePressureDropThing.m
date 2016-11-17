%孔板压力降分析
Dpipe = 0.106;
Dv1 = 0.45;
d = 0.106/2;
vPipe = 14.5;
density = 0.211;
vVessel = (Dpipe^2./Dv1^2)*vPipe;
kV = orificeKs(Dv1,d);
pressureDrop(1) = pressureDrop_kPa(vVessel,density,kV);
kP = orificeKs(Dpipe,d);
pressureDrop(2) = pressureDrop_kPa(vPipe,density,kP);


fprintf('管道流速:%g m/s 缓冲罐流速:%g m/s\n',vPipe,vVessel);
fprintf('管道孔板k:%g 缓冲罐孔板k:%g \n',kP,kV);
fprintf('管道孔板压降:%g kPa 缓冲罐孔板压降:%g  kPa\n',pressureDrop(2),pressureDrop(1));

%%在同等质量流量下，
d = 0.106/2;
dPipeOriginal = 0.106;
D = [0.6:0.01:4].*dPipeOriginal;
vOriginal=14.5;
v = (dPipeOriginal.^2./D.^2).*vOriginal;
k = orificeKs(D,d);
pressureDrop = pressureDrop_kPa(v,density,k);
figure
plot(D,pressureDrop);
set(gcf,'color','w');