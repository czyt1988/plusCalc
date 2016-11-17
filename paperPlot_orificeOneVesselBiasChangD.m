%% 绘制缓冲罐内含孔板结构随着孔板孔径大小变化，缓冲罐后的气流脉动变化情况
clc;
close all;
clear;
currentPath = fileparts(mfilename('fullpath'));
%含孔板缓冲罐的气流脉动计算
%   Detailed explanation goes here
%           |  L2
%        l  | 
%  bias2 ___|_______________
%       |         |        |
%   Dv2 |     V2   d    V1 |  Dv1
%       |_________|________|
%                    l  |   bias1  
%                       |
%              inlet:   | L1 Dpipe 

massFlow = load(fullfile(currentPath,'mass_flow_0.1478_NorthZone.txt'));
N = 4096;
time = massFlow(1:N,1);
massFlowRaw = massFlow(1:N,2);
Fs = 1/(time(2)-time(1));
[FreRaw,AmpRaw,PhRaw,massFlowE] = fun_fft(detrend(massFlowRaw),Fs);
% 提取主要频率
[pks,locs] = findpeaks(AmpRaw);
Fre = FreRaw(locs);
massFlowE = massFlowE(locs);
temp = Fre < 100;%(Fre<29) | (Fre>30 & Fre < 100);
Fre = Fre(temp);
massFlowE = massFlowE(temp);
isDamping = 1;
%绘图参数
isXLength = 1;
isShowStraightPipe=1;%是否显示直管
isShowOnlyVessel=1;%是否显示无内件缓冲罐

opt.acousticVelocity = 345;%声速
opt.isDamping = 1;%是否计算阻尼
opt.coeffFriction = 0.05;%管道摩察系数
opt.meanFlowVelocity = 14.6;%管道平均流速
opt.isUseStaightPipe = 1;%计算容器传递矩阵的方法
density = 2.09;%气体密度
dcpss = getDefaultCalcPulsSetStruct();
dcpss.calcSection = [0.4,0.5];
dcpss.isHp = 0;
dcpss.f_pass = 7;%通过频率5Hz
dcpss.f_stop = 5;%截止频率3Hz
dcpss.rp = 0.1;%边带区衰减DB数设置
dcpss.rs = 30;%截止区衰减DB数设置
isUseStaightPipe = 1;
d = 0.106/3:0.01:0.4;
for i=1:length(d)
    para(i).opt = opt;
    para(i).L1 = 2.94;
    para(i).L2 = 10;
    para(i).Dpipe = 0.106;
    para(i).Dv1 = 0.45;
    para(i).Dv2 = 0.45;
    para(i).l = 0.115;
    para(i).Lv1 = 0.4755;
    para(i).Lv2 = 0.4755;
    para(i).d = d(i);
    para(i).bias1 = 0.3;
    para(i).bias2 = 0.3;
    para(i).sectionL1 = 0:1:para(i).L1;
    para(i).sectionL2 = 0:1:para(i).L2;
end
pressureDrop = [];
for i = 1:length(para)
    pressure1 = [];
    pressure2 = [];
    
    [pressure1,pressure2] = ...
        vesselBiasHaveOrificePulsationCalc(massFlowE,Fre,time,...
            para(i).L1,para(i).L2,...
            para(i).Lv1,para(i).Lv2,para(i).l,para(i).Dpipe,para(i).Dv1,para(i).Dv2,...
            para(i).d,para(i).bias1,para(i).bias2,...
            para(i).sectionL1,para(i).sectionL2,...
            'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,...
            'meanFlowVelocity',opt.meanFlowVelocity,'isUseStaightPipe',isUseStaightPipe);
    plus1{i} = calcPuls(pressure1,dcpss);
    plus2{i} = calcPuls(pressure2,dcpss);
    plus{i} = [plus1{i},plus2{i}];
    %计算缓冲罐流速
    vVessel = (para(i).Dpipe^2./para(i).Dv1^2)*opt.meanFlowVelocity;
    k = orificeKs(para(i).Dv1,para(i).d);
    pressureDrop(i) = pressureDrop_kPa(vVessel,density,k);
    if i == 1
        [pressure1OV,pressure2OV] = oneVesselBiasPulsationCalc(massFlowE,Fre,time,...
            para(i).L1,para(i).L2,...
            para(i).Lv1+para(i).Lv2,para(i).l,para(i).Dpipe,para(i).Dv1,para(i).bias1,para(i).bias2,...
            para(i).sectionL1,para(i).sectionL2,...
            'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,...
            'meanFlowVelocity',opt.meanFlowVelocity,'isUseStaightPipe',isUseStaightPipe);
        plus1OV = calcPuls(pressure1OV,dcpss);
        plus2OV = calcPuls(pressure2OV,dcpss);
        plusOV = [plus1OV,plus2OV];
        maxPlus1OV = max(plus1OV);
        maxPlus2OV = max(plus2OV);
    end
    if isempty(plus1{i})
        maxPlus1(i) = nan;
    else
        maxPlus1(i) = max(plus1{i});
    end

    if isempty(plus2{i})
        maxPlus2(i) = nan;
    else
        maxPlus2(i) = max(plus2{i});
    end

    
    newSectionL2 = para(i).L1 + 2*para(i).l+para(i).Lv1 + para(i).Lv2 + para(i).sectionL2;

    realLengthSection{i} = [para(i).sectionL1,newSectionL2];

    L_straight(i) = para(i).L1 + para(i).L2 + para(i).Lv1 + para(i).Lv2 + 2*para(i).l;
    temp = find(realLengthSection{i}>para(i).L1);
    sepratorIndex(i) = temp(1);
    %temp = fun_straightPipe(massFlowE1,L_straight(i),para(i).Dpipe,para(i).opt,realLengthSection{i});
    temp = straightPipePulsationCalc(massFlowE,Fre,time,L_straight(i),realLengthSection{i}...
	,'d',para(i).Dpipe,'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,'meanFlowVelocity',opt.meanFlowVelocity);
    temp = calcPuls(temp,dcpss);
    maxPlus1Straight(i) = max(temp(1:sepratorIndex(i)));
    maxPlus3Straight(i) = max(temp(sepratorIndex(i):end));
    plusStraight{i} = temp;
    temp2 = plus{i};
    temp(temp<1e-4) = 1;
    temp2(temp2<1e-4) = 1;
    reduceRate{i} = (temp - temp2)./temp;
    
    maxPlus1Straight(i) = max(temp(1:sepratorIndex(i)));
    maxPlus3Straight(i) = max(temp(sepratorIndex(i):end));

end

%% 绘图-孔径变化对罐前最大压力脉动的影响
figure
x = d.*1000;
y = maxPlus1./1000;
hold on;
%h(1) = plot(x,y,'-r');
[AX H1 H2] = plotyy(x,y,x,pressureDrop);
set(H2,'LineStyle','-.','LineWidth',1.5);
set(H1,'LineWidth',1.5);
set(get(AX(2),'Ylabel'),'string','pressure drop(kPa)');%设置y2对应坐标轴内容
ax = axis();
h(1) = plot([ax(1) ax(2)],[maxPlus1OV./1000 maxPlus1OV./1000],'--b');
h(2) = plot([0.106 0.106].*1000,[25 35],'--');
xlabel('orifice diameter(mm)');
ylabel(sprintf('max peak-to-peak\n pressure pulsation(kPa)'));
box on;
grid on;
set(gcf,'Units','pixels','position',[200,200,350,250]);
set(AX,'Units','pixels','position',[65,45,230,190]);
set(AX,'FontName','times new roman');
%export_fig(gcf,'d:\cloudDisk\share\【论文】孔板\图\孔板孔径变化对罐前气流脉动的影响3en.png','-m3.5','-painters');
set(gcf,'color','w');

%% 绘图-孔径变化对罐后最大压力脉动的影响
figure
x = d.*1000;
y = maxPlus2./1000;
hold on;
[AX H1 H2] = plotyy(x,y,x,pressureDrop);
set(H2,'LineStyle','-.','LineWidth',1.5);
set(H1,'LineWidth',1.5);
set(get(AX(2),'Ylabel'),'string','pressure drop(kPa)');%设置y2对应坐标轴内容
ax = axis();
h(1) = plot([ax(1) ax(2)],[maxPlus2OV./1000 maxPlus2OV./1000],'--b');
h(2) = plot([0.106 0.106].*1000,[0 7],'--');
xlabel('orifice diameter(mm)');
ylabel(sprintf('max peak-to-peak\n pressure pulsation(kPa)'));
box on;
grid on;
set(gcf,'Units','pixels','position',[200,200,350,250]);
set(AX,'Units','pixels','position',[65,45,230,190]);
set(AX,'FontName','times new roman');
set(gcf,'color','w');
%export_fig(gcf,'d:\cloudDisk\share\【论文】孔板\图\孔板孔径变化对罐后气流脉动的影响3en.png','-m3.5','-painters')