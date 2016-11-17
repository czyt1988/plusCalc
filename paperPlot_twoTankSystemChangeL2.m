%% 绘制理论分析的双缓冲罐中间接管的间距变化的影响总管长保持不变
% 双罐的体积比例不一样对气流脉动的影响
% 
clc;
close all;
clear;
currentPath = fileparts(mfilename('fullpath'));
%  长度 L1     l    Lv1   l   L2  l    Lv2   l     L3
%              __________         __________
%             |          |       |          |
%  -----------|          |-------|          |-------------
%             |__________|       |__________|  
%  直径 Dpipe       Dv1    Dpipe       Dv2          Dpipe
[time,massFlow,Fre,massFlowE ] = getMassFlowData('N',4096,'isfindpeaks',1);
Fs = 1/(time(2)-time(1));
% [los]=find(Fre>19 & Fre <21);
% Fre(los) = Fre(los)./1.5;

temp =Fre>5 & Fre<20.1;%Fre<20 | (Fre>22&Fre<60);% Fre>9.5 & Fre<10.1;%Fre<20 | (Fre>22&Fre<80);
Fre = Fre(temp);

massFlowE = massFlowE(temp);
opt.acousticVelocity = 345;%声速
opt.isDamping = 1;%是否计算阻尼
opt.coeffFriction = 0.05;%管道摩察系数
opt.meanFlowVelocity = 14.6;%管道平均流速
opt.isUseStaightPipe = 1;%计算容器传递矩阵的方法

dcpss = getDefaultCalcPulsSetStruct();
dcpss.calcSection = [0.4,0.5];
dcpss.isHp = 0;
dcpss.f_pass = 7;%通过频率5Hz
dcpss.f_stop = 5;%截止频率3Hz
dcpss.rp = 0.1;%边带区衰减DB数设置
dcpss.rs = 30;%截止区衰减DB数设置
%%参数设定
L2 = 0:0.2:8;
Ltotal = 28;%总管长
for i=1:length(L2)
    paraDV(i).opt = opt;
    paraDV(i).L1 = 13;%L1(m)
    paraDV(i).L2 = L2(i);%缓冲罐中间连接管道的长度（m）
    paraDV(i).Dpipe = 0.157;%管道直径（m）
    paraDV(i).Dv1 = 0.5;
    paraDV(i).Dv2 = 0.5;    
    paraDV(i).l = 0.115;%0.115;%缓冲罐前管道的长度(m)   
    paraDV(i).Lv1 = 1;%[[0.157,0.25,0.5,0.75],[1:0.25:5]];%第一个缓冲罐的直径（m）
    paraDV(i).Lv2 = 1;%0.5;%第二个缓冲罐的直径（m）
    paraDV(i).L3 = Ltotal-L2(i)-paraDV(i).L1-4*paraDV(i).l-paraDV(i).Lv1-paraDV(i).Lv2;
    paraDV(i).sectionL1 = 0:1:paraDV(i).L1;
    paraDV(i).sectionL2 = 0:1:paraDV(i).L2;
    paraDV(i).sectionL3 = 0:1:paraDV(i).L3;
    VDV1(i) = (pi*paraDV(i).Dv1^2/4)*paraDV(i).Lv1+2*((3/5)*(pi*paraDV(i).Dv1^2/4)*paraDV(i).l);
    VDV2(i) = (pi*paraDV(i).Dv2^2/4)*paraDV(i).Lv2+2*((3/5)*(pi*paraDV(i).Dv2^2/4)*paraDV(i).l);
    VDV(i) = VDV1(i)+VDV2(i);
%     paraOV(i).L1 = 13;
%     paraOV(i).L2 = 13;
%     paraOV(i).Lv = 2;
%     paraOV(i).Dv = 0.5;
%     paraOV(i).l = 0.115;
%     paraOV(i).Dpipe = 0.157;%管道直径（m）
%     paraOV(i).sectionL1 = 0:1:paraOV(i).L1;
%     paraOV(i).sectionL2 = 0:1:paraOV(i).L2;
%     VOV(i) = (pi*paraOV(i).Dv^2/4)*paraOV(i).Lv+2*((3/5)*(pi*paraOV(i).Dv^2/4)*paraOV(i).l);
    %计算气流脉动
    [pressure1DV,pressure2DV,pressure3DV] = ...
                doubleVesselPulsationCalc(massFlowE,Fre,time,...
                paraDV(i).L1,paraDV(i).L2,paraDV(i).L3,...
                paraDV(i).Lv1,paraDV(i).Lv2,paraDV(i).l,paraDV(i).Dpipe,paraDV(i).Dv1,paraDV(i).Dv2,...
                paraDV(i).sectionL1,paraDV(i).sectionL2,paraDV(i).sectionL3,...
                'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,...
                'meanFlowVelocity',opt.meanFlowVelocity);
    plus1DV{i} = calcPuls(pressure1DV,dcpss);
    plus2DV{i} = calcPuls(pressure2DV,dcpss);
    plus3DV{i} = calcPuls(pressure3DV,dcpss);
    aheadMaxPlusDV(i) = max(plus1DV{i});%罐前最大气流脉动
    afterMaxPlusDV(i) = max(plus3DV{i});%罐后最大气流脉动

end
paraOV(1).L1 = 13;
paraOV(1).Lv = 2;
paraOV(1).l = 0.115;
paraOV(1).Dpipe = 0.157;
paraOV(1).Dv = 0.5;
paraOV(1).L2 = Ltotal - paraOV(1).L1 - 2*paraOV(1).l - paraOV(1).Lv;
paraOV(1).sectionL1 = 0:1:paraOV(1).L1;
paraOV(1).sectionL2 = 0:1:paraOV(1).L2;

[pressure1OV,pressure2OV] = oneVesselPulsationCalc(massFlowE,Fre,time,...
            paraOV(1).L1,paraOV(1).L2,paraOV(1).Lv,paraOV(1).l,paraOV(1).Dpipe,paraOV(1).Dv ...
            ,paraOV(1).sectionL1,paraOV(1).sectionL2 ...
            ,'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,'meanFlowVelocity',opt.meanFlowVelocity);
        plus1OV{1} = calcPuls(pressure1OV,dcpss);
        plus2OV{1} = calcPuls(pressure2OV,dcpss);
        aheadMaxPlusOV(1) = max(plus1OV{1});%罐前最大气流脉动
        afterMaxPlusOV(1) = max(plus2OV{1});%罐后最大气流脉动

%% 计算完成，进行绘图

%绘制罐前压力脉动
scale = 1;
figure
hold on;
h1 = plot(L2,aheadMaxPlusDV./1000,'--r');
set(h1,'Marker','x');
h2 = plot(L2,afterMaxPlusDV./1000,'--r');
set(h2,'Marker','<');
xy = axis;
hL1 = plot([xy(1,1:2)],[aheadMaxPlusOV(1),aheadMaxPlusOV(1)]./1000,'-b');
hL2 = plot([xy(1,1:2)],[afterMaxPlusOV(1),afterMaxPlusOV(1)]./1000,'-.b');
xlabel('L2 distance(m)');
ylabel('gas pulsation(kPa)');
legend([h1,h2],{'two tank system ahead','two tank system after'});
text(xy(2)-4,aheadMaxPlusOV(1)./1000+0.5,'single tank system ahead');
text(xy(2)-4,afterMaxPlusOV(1)./1000+0.5,'single tank system after');
box on;
set(gcf,'color','w');
set(gcf,'position',[200,200,400*scale,300*scale]);

