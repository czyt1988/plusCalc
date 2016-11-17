%% 绘制缓冲罐内含孔板结构随着孔板孔径大小变化，缓冲罐后的气流脉动变化情况
clc;
close all;
clear;
currentPath = fileparts(mfilename('fullpath'));
%含孔板缓冲罐的气流脉动计算
%   Detailed explanation goes here
%  L1  l    Lv1     Lv2   l  L2
%        __________________
%       |         |        |
% ------|     V1   d    V2 |-------
%       |_________|________|
%    Dpipe  Dv1    d   Dv2    Dpipe 

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
temp = Fre<100;%(Fre<29) | (Fre>30 & Fre < 100);
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
dcpss = getDefaultCalcPulsSetStruct();
dcpss.calcSection = [0.4,0.5];
dcpss.isHp = 0;
dcpss.f_pass = 7;%通过频率5Hz
dcpss.f_stop = 5;%截止频率3Hz
dcpss.rp = 0.1;%边带区衰减DB数设置
dcpss.rs = 30;%截止区衰减DB数设置
LTotal = 0.951;
Lv1 = 0.1:0.05:(LTotal-0.1);
d = 0.106/3:0.005:0.12;
d(15) = 0.106;%全流通
d(5) = 0.053;%0.5D流通
d(10) = 0.0795;%0.75D流通
for i=1:length(Lv1)%y轴
    for j=1:length(d)%x轴
        para(i,j).opt = opt;
        para(i,j).L1 = 2.94;
        para(i,j).L2 = 10;
        para(i,j).Dpipe = 0.106;
        para(i,j).Dv1 = 0.45;
        para(i,j).Dv2 = 0.45;
        para(i,j).l = 0.115;
        para(i,j).Lv1 = Lv1(i);
        para(i,j).Lv2 = LTotal - Lv1(i);
        para(i,j).d = d(j);
        para(i,j).sectionL1 = 0:1:para(i,j).L1;
        para(i,j).sectionL2 = 0:1:para(i,j).L2;
    end
end

for i = 1:size(para,1)
    for j = 1:size(para,2)
        pressure1 = [];
        pressure2 = [];

        [pressure1,pressure2] = ...
            vesselHaveOrificePulsationCalc(massFlowE,Fre,time,...
            para(i,j).L1,para(i,j).L2,...
            para(i,j).Lv1,para(i,j).Lv2,para(i,j).l,para(i,j).Dpipe,para(i,j).Dv1,para(i,j).Dv2,...
            para(i,j).d,...
            para(i,j).sectionL1,para(i,j).sectionL2,...
            'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,...
            'meanFlowVelocity',opt.meanFlowVelocity);
        plus1{i,j} = calcPuls(pressure1,dcpss);
        plus2{i,j} = calcPuls(pressure2,dcpss);
        plus{i,j} = [plus1{i,j},plus2{i,j}];
        if i == 1
            [pressure1OV,pressure2OV] = oneVesselPulsationCalc(massFlowE,Fre,time,...
                para(i,j).L1,para(i,j).L2,...
                para(i,j).Lv1+para(i,j).Lv2,para(i,j).l,para(i,j).Dpipe,para(i,j).Dv1,...
                para(i,j).sectionL1,para(i,j).sectionL2,...
                'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,...
                'meanFlowVelocity',opt.meanFlowVelocity);
            plus1OV = calcPuls(pressure1OV,dcpss);
            plus2OV = calcPuls(pressure2OV,dcpss);
            plusOV = [plus1OV,plus2OV];
            maxPlus2OV = max(plus2OV);
            maxPlus1OV = max(plus1OV);
        end
        if isempty(plus1{i,j})
            maxPlus1(i,j) = nan;
        else
            maxPlus1(i,j) = max(plus1{i,j});
        end

        if isempty(plus2{i,j})
            maxPlus2(i,j) = nan;
        else
            maxPlus2(i,j) = max(plus2{i,j});
        end


        newSectionL2 = para(i,j).L1 + 2*para(i,j).l+para(i,j).Lv1 + para(i,j).Lv2 + para(i,j).sectionL2;

        realLengthSection{i,j} = [para(i,j).sectionL1,newSectionL2];

        L_straight(i,j) = para(i,j).L1 + para(i,j).L2 + para(i,j).Lv1 + para(i,j).Lv2 + 2*para(i,j).l;
        temp = find(realLengthSection{i,j}>para(i,j).L1);
        sepratorIndex(i,j) = temp(1);
        %temp = fun_straightPipe(massFlowE1,L_straight(i,j),para(i,j).Dpipe,para(i,j).opt,realLengthSection{i,j});
        temp = straightPipePulsationCalc(massFlowE,Fre,time,L_straight(i,j),realLengthSection{i,j}...
        ,'d',para(i,j).Dpipe,'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,'meanFlowVelocity',opt.meanFlowVelocity);
        temp = calcPuls(temp,dcpss);
        maxPlus1Straight(i,j) = max(temp(1:sepratorIndex(i,j)));
        maxPlus3Straight(i,j) = max(temp(sepratorIndex(i,j):end));
        plusStraight{i,j} = temp;
        temp2 = plus{i,j};
        temp(temp<1e-4) = 1;
        temp2(temp2<1e-4) = 1;
        reduceRate{i,j} = (temp - temp2)./temp;

        maxPlus1Straight(i,j) = max(temp(1:sepratorIndex(i,j)));
        maxPlus3Straight(i,j) = max(temp(sepratorIndex(i,j):end));
    end
end

%% 绘图
%% 后端云图
if 1
    figure
    [X,Y] = meshgrid(d*1000,Lv1*1000);
    Z = (maxPlus2OV-maxPlus2)./maxPlus2OV*100;
    if 1
        [C1,h1] = contourf(X,Y,Z);
        clabel(C1,h1);
    else
        surfc(X,Y,Z);
    end
    xlabel('orifice diameter(mm)');
    ylabel('orifice location(mm)');
    set(gcf,'color','w');
    colorbar
    set(gcf,'position',[200,200,380,280]);
end


%% 前端云图
if 1
    figure
    [X,Y] = meshgrid(d*1000,Lv1*1000);
    Z = (maxPlus1OV-maxPlus1)./maxPlus1OV*100;
    if 1
        [C1,h1] = contourf(X,Y,Z);
        clabel(C1,h1);
    else
        surfc(X,Y,Z);
    end
    xlabel('orifice diameter(mm)');
    ylabel('orifice location(mm)');
    set(gcf,'color','w');
    colorbar
    set(gcf,'position',[200,200,380,280]);
end


%% 后端
if 0
    figure
    plot(Lv1*1000,maxPlus2(:,5)./1000,'LineWidth',1.5);text(800,maxPlus2(end,5)./1000-0.2,sprintf('53mm'),'FontName','times new roman');
    hold on;
    plot(Lv1*1000,maxPlus2(:,10)./1000,'--','LineWidth',1.5);text(800,maxPlus2(end,10)./1000-0.2,sprintf('79.5mm'),'FontName','times new roman');
    plot(Lv1*1000,maxPlus2(:,15)./1000,':','LineWidth',1.5);text(800,maxPlus2(end,15)./1000-0.2,sprintf('106mm'),'FontName','times new roman');
    text(250,5.8,sprintf('surge volume'),'FontName','times new roman');
    xlabel('orifice location(mm)');
    ylabel(sprintf('max peak-to-peak\n pressure pulsation(kPa)'));
    set(gcf,'position',[200,200,350,250]);
    set(gca,'FontName','times new roman');
    ax = axis();
    plot([ax(1),ax(2)],[maxPlus2OV,maxPlus2OV]./1000,'--');
    set(gcf,'color','w');
end

%% 前端
if 0
    figure
    plot(Lv1*1000,maxPlus1(:,5)./1000,'LineWidth',1.5);text(300,22.92,sprintf('53mm'),'FontName','times new roman');
    hold on;
    plot(Lv1*1000,maxPlus1(:,10)./1000,'--','LineWidth',1.5);text(200,24.46,sprintf('79.5mm'),'FontName','times new roman');
    plot(Lv1*1000,maxPlus1(:,15)./1000,':','LineWidth',1.5);text(250,26.67,sprintf('106mm'),'FontName','times new roman');
    text(250,29.5,sprintf('surge volume'),'FontName','times new roman');
    xlabel('orifice location(mm)');
    ylabel(sprintf('max peak-to-peak\n pressure pulsation(kPa)'));
    set(gcf,'position',[200,200,350,250]);
    set(gca,'FontName','times new roman');
    ax = axis();
    plot([ax(1),ax(2)],[maxPlus1OV,maxPlus1OV]./1000,'--');
    set(gcf,'color','w');
end
