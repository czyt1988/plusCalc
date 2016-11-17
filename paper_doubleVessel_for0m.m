%% 双容间距变化对比
% 与实验相同工况的模拟进行对比
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

%% 模拟结果
%1，2,3 -》 0m,0.5m,0.9m
simulationData{1} = [8.58886963000000,9.70715454000000,9.66618115000000,10.3086987300000,10.1641142600000,9.60220923000000,9.50955761000000,9.62627491000000,9.48254126000000,8.80504687000000,7.56495191000000,5.53023132000000,2.81892358700000,1.34184783900000,2.60702405000000,3.43060047000000,4.07793481000000,4.42835828000000,4.46505384000000,4.39963598000000,4.24979663000000,4.18714112000000,4.06717078000000,3.73728027000000,3.09321313000000,2.16947571000000,1.04654071000000,4.77363361700000e-06];
%% 质量流量获取
[time,massFlow,Fre,massFlowE ] = getMassFlowData('N',4096,'isfindpeaks',1);

[los]=find(Fre>19 & Fre <21);
Fre(los) = Fre(los)./1.5;

temp = Fre<80;%Fre<20 | (Fre>22&Fre<80);
Fre = Fre(temp);

massFlowE = massFlowE(temp);

%% 基本参数设置

isSaveXls = 0;
isXLength = 1;%如果X轴是距离将会按照实际尺寸来绘制。
isShowStraightPipe = 0;%是否加入直管对比
titleText = '进口管长1m总长26米，两缓冲罐间距不同对脉动的影响';
isShowL2Plus = 0;%不显示L2的脉动
isDamping = 1;
%%
inputData = [...
    %1  2  3  4     5   6   7  8  9  
    %L1,L2,L3,Dpipe,Dv1,Dv2,l,Lv1,Lv2
     13,0,13,0.157,0.5,0.5,0.115,1,1 ...
%     ;13,0.5,13,0.157,0.5,0.5,0.115,1,1 ...
%     ;13,0.9,13,0.157,0.5,0.5,0.115,1,1 ...
];
for i=1:size(inputData,1)
    xTick(i) = inputData(i,2);
    xTickLabel{i} = sprintf('间距%gm',inputData(i,2));
end
%opt.frequency = 10;%脉动频率
opt.acousticVelocity = 345;%声速
opt.isDamping = isDamping;%是否计算阻尼
opt.coeffFriction = 0.05;%管道摩察系数
opt.meanFlowVelocity = 14.6;%管道平均流速
opt.isUseStaightPipe = 1;%计算容器传递矩阵的方法


%% 计算不同缓冲罐间隔下的脉动压力

xlsDataToBeWrite = {};

dcpss = getDefaultCalcPulsSetStruct();
dcpss.calcSection = [0.4,0.5];
%dcpss.fs = Fs;
dcpss.isHp = 0;
dcpss.f_pass = 7;%通过频率5Hz
dcpss.f_stop = 5;%截止频率3Hz
dcpss.rp = 0.1;%边带区衰减DB数设置
dcpss.rs = 30;%截止区衰减DB数设置
for i = 1:length(para)
    pressure1 = [];
    pressure2 = [];
    pressure3 = [];
    [pressure1,pressure2,pressure3] = ...
        doubleVesselPulsationCalc(massFlowE,Fre,time,...
        para(i).L1,para(i).L2,para(i).L3,...
        para(i).Lv1,para(i).Lv2,para(i).l,para(i).Dpipe,para(i).Dv1,para(i).Dv2,...
        para(i).sectionL1,para(i).sectionL2,para(i).sectionL3,...
        'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,...
        'meanFlowVelocity',opt.meanFlowVelocity);
    plus1{i} = calcPuls(pressure1,dcpss);
    plus2{i} = calcPuls(pressure2,dcpss);
    plus3{i} = calcPuls(pressure3,dcpss);

    plus{i} = [plus1{i},plus2{i},plus3{i}];
    if isSaveXls
        tableHeader = {};
        for k=1:length(para(i).sectionL1)
            tableHeader = cellPush2Right(tableHeader,{sprintf('L1:%g',para(i).sectionL1(k))});
        end
        for k=1:length(para(i).sectionL2)
            tableHeader = cellPush2Right(tableHeader,{sprintf('L2:%g',para(i).sectionL2(k))});
        end
        for k=1:length(para(i).sectionL3)
            tableHeader = cellPush2Right(tableHeader,{sprintf('L3:%g',para(i).sectionL3(k))});
        end
        xlsDataToBeWrite = toOneCell(tableHeader,pressure1,pressure2,pressure3);
        xlswrite(fullfile(currentPath,[name{i},'.xls']),xlsDataToBeWrite);
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

    if isempty(plus3{i})
        maxPlus3(i) = nan;
    else
        maxPlus3(i) = max(plus3{i});
    end
    
    newSectionL2 = para(i).L1 + 2*para(i).l+para(i).Lv1 + para(i).sectionL2;
    newSectionL3 = para(i).L1+4*para(i).l+para(i).Lv1+para(i).Lv2+para(i).L2 + para(i).sectionL3;
    realLengthSection{i} = [para(i).sectionL1,newSectionL2,newSectionL3];

    L_straight(i) = para(i).L1 + para(i).L2 + para(i).L3 + para(i).Lv1 + para(i).Lv2 + 4*para(i).l;
    temp = find(realLengthSection{i}>para(i).L1);
    sepratorIndex(i) = temp(1);

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
if isShowStraightPipe
    maxPlus1 = [maxPlus1Straight(1,:);maxPlus1];
    maxPlus3 = [maxPlus3Straight(1,:);maxPlus3];
end

%% 绘图
marker_style = {'-o','-d','-<','-s','->','-<','-p','-*','-v','-^','-+','-x','-h'};
color_style = [...
    245,18,103;...
    36,100,196;...
    18,175,134;...
    237,144,10;... 
    131,54,229;...    
    255,99,56;...
    ]./255;
market_style_length = length(marker_style);
color_style_length = size(color_style,1);
marker_index = 0;
color_index = 0;

figure
marker_index = 0;
color_index = 0;
h = [];
textLegend = {};
hold on;
if isShowStraightPipe
    h(1) = plot(sectionL{i},plusStraight{i}./1000,'-r','LineWidth',1.5);
    textLegend = ['直管',name{i}];
else
    textLegend = name;
end
for i = 1:length(para)
    marker_index = marker_index + 1;
    if marker_index > market_style_length
        marker_index = 1;
    end
    color_index = color_index + 1;
    if color_index > color_style_length
        color_index = 1;
    end
    if isShowL2Plus
        Y = plus{i};
    else
        Y = [plus1{i},plus3{i}];
    end
    Y = Y./1000;
    Y=Y';
    if isXLength
        newSectionL2 = para(i).L1 + 2*para(i).l+para(i).Lv1 + para(i).sectionL2;
        newSectionL3 = para(i).L1+4*para(i).l+para(i).Lv1+para(i).Lv2+para(i).L2 + para(i).sectionL3;
        if isShowL2Plus
            X = realLengthSection{i};
        else
            X = [para(i).sectionL1,newSectionL3];
        end
    else
        X = 1:length(Y);
    end
    if isShowStraightPipe
        h(i+1) = plot(X,Y,marker_style{marker_index},'color',color_style(color_index,:),'LineWidth',1.5);    
    else
        h(i) = plot(X,Y,marker_style{marker_index},'color',color_style(color_index,:),'LineWidth',1.5);
        str = strrep(marker_style{marker_index},'-',':');
        plot(X,simulationData{i},str,'color',color_style(color_index,:),'LineWidth',1.5);
    end
end
legend(h,textLegend,0);
grid on;
if isXLength
    xlabel('距离(m)');
else
    xlabel('测点');
end
ylabel('脉动压力峰峰值(kPa)');
title(titleText);
set(gcf,'color','w');


h = [];
marker_index = 1;
color_index = 1;    
barMaxPlusBefore = maxPlus1./1000;
barMaxPlusMiddle = maxPlus2./1000;
barMaxPlusAfter = maxPlus3./1000;  
figure
hold on;
if isempty(xTick)
    X = 1:length(para); 
else
    X = xTick;
end
h(1) = plot(X,barMaxPlusBefore,marker_style{marker_index},'color',color_style(color_index,:),'LineWidth',1.5);
% h(2) = plot([X(1),X(end)],[maxL1./1000,maxL1./1000],'-','color',color_style(color_index,:));
marker_index = 2;
color_index = 2; 
h(3) = plot(X,barMaxPlusAfter,marker_style{marker_index},'color',color_style(color_index,:),'LineWidth',1.5);
% h(4) = plot([X(1),X(end)],[maxL2./1000,maxL2./1000],'-','color',color_style(color_index,:));
marker_index = 3;
color_index = 3; 
h(5) = plot(X,barMaxPlusMiddle,marker_style{marker_index},'color',color_style(color_index,:),'LineWidth',1.5);
title('缓冲罐间隔对气流脉动的影响');
ylabel('脉动压力最大值(kPa)');
xlabel('缓冲罐间距(m)');
if isempty(xTick)
    set(gca,'XTick',0:size(inputData,1));
    set(gca,'XTickLabel',xTickLabel);
end
legend(h,{'缓冲罐前管道脉动最大值','单缓冲罐前管道脉动最大值','缓冲罐后管道脉动最大值','单缓冲罐后管道脉动最大值','连接管道脉动最大值'});
grid on;
set(gcf,'color','w');

