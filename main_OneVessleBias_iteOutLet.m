%% 错位缓冲罐迭代
clc;
close all;
clear;
currentPath = fileparts(mfilename('fullpath'));
%错位缓冲罐的气流脉动计算
%   Detailed explanation goes here
%           |  L2
%        l  |     Lv    outlet
%   bias2___|_______________
%       |                   |
%       |lv2  V          lv1|  Dv
%       |___________________|
%                    l  |   bias1  
%                       |
%              inlet:   | L1 Dpipe 
%容器的传递矩阵
isOpening = 	0;%管道闭口
%rpm = 300;outDensity = 1.9167;multFre=[10,20,30];%环境25度绝热压缩到0.2MPaG的温度对应密度
rpm = 420;outDensity = 1.5608;multFre=[14,28,42];%环境25度绝热压缩到0.15MPaG的温度对应密度
Fs = 4096;
[massFlowRaw,time,~,opt.meanFlowVelocity] = massFlowMaker(0.25,0.098,rpm...
	,0.14,1.075,outDensity,'rcv',0.15,'k',1.4,'pr',0.15,'fs',Fs,'oneSecond',6);

[FreRaw,AmpRaw,PhRaw,massFlowERaw] = frequencySpectrum(detrend(massFlowRaw),Fs);
% 提取主要频率
[pks,locs] = findpeaks(AmpRaw,'SORTSTR','descend');
Fre = FreRaw(locs);
massFlowE = massFlowERaw(locs);
temp = 1:20;%(Fre<29) ;%| (Fre>30 & Fre < 100);
massFlowE = massFlowE(temp);
Fre = Fre(temp);

isDamping = 1;
%绘图参数
isXShowRealLength = 1;
isShowStraightPipe=1;%是否显示直管
isShowOnlyVessel=1;%是否显示无内件缓冲罐

calcDatas = {};%绘图数据列
opt.frequency = 10;%脉动频率
opt.acousticVelocity = 345;%声速
opt.isDamping = isDamping;%是否计算阻尼
opt.coeffDamping = nan;%阻尼
opt.coeffFriction = 0.04;%管道摩察系数
% opt.meanFlowVelocity = 14.5;%管道平均流速
opt.isUseStaightPipe = 0;%计算容器传递矩阵的方法
opt.mach = opt.meanFlowVelocity / opt.acousticVelocity;
opt.notMach = 1;

L1 = 1.5;%L1(m)
L2 = 6;%缓冲罐中间连接管道的长度（m） 
Lv = 1.1;%管道直径（m）
l = 0.01;
Dpipe = 0.098;%0.115;%缓冲罐前管道的长度(m)   
Dv = 0.372;%[[0.157,0.25,0.5,0.75],[1:0.25:5]];%第一个缓冲罐的直径（m）
lv1 = 0.45;%
lv2 = 0.45;
Dbias = 0;%


sectionL1 = 0:0.2:L1;
sectionL2 = 0:0.2:L2;


dcpss = getDefaultCalcPulsSetStruct();
dcpss.calcSection = [0.3,0.7];
dcpss.fs = Fs;
dcpss.isHp = 0;
dcpss.f_pass = 7;%通过频率5Hz
dcpss.f_stop = 5;%截止频率3Hz
dcpss.rp = 0.1;%边带区衰减DB数设置
dcpss.rs = 30;%截止区衰减DB数设置

dataCount = 1;
outletPressure = 1000+1000i;
for i = 1:length(outletPressure)
    %计算直管的长度
    straightPipeLength = L1 + 2*l+Lv + L2;
    %计算直管的分段
    straightPipeSection = [sectionL1,...
                            L1 + 2*l+Lv + sectionL2];
    if isXShowRealLength
        X = straightPipeSection;
    else
        X = 1:length(Y);
    end
    %先对直管进行计算
    if i==1
        %计算直管
        %直管总长
        
        newSectionL2 = L1 + 2*l+Lv + sectionL2;
        temp = find(straightPipeLength>L1);%找到缓冲罐所在的索引
        sepratorIndex = temp(1);
        temp = straightPipePulsationCalc(massFlowE,Fre,time,straightPipeLength,straightPipeSection...
        ,'d',Dpipe,'a',opt.acousticVelocity,'isDamping',opt.isDamping...
        ,'friction',opt.coeffFriction,'meanFlowVelocity',opt.meanFlowVelocity...
        ,'m',opt.mach,'notMach',opt.notMach...
        ,'isOpening',isOpening...
        );
        plusStraight = calcPuls(temp,dcpss);
        maxPlus1Straight(i) = max(plusStraight(1:sepratorIndex(i)));
        maxPlus2Straight(i) = max(plusStraight(sepratorIndex(i):end));

        calcDatas{dataCount,1} = X;
        calcDatas{dataCount,2} = plusStraight;
        calcDatas{dataCount,3} = sprintf('直管');%具体描述可以在这改
        dataCount = dataCount + 1;
    end

    pressure1 = [];
    pressure2 = [];

    % [pressure1,pressure2] = ...
    %     vesselBiasPulsationCalc(massFlowE,Fre,time,...
    %     L1,L2,...
    %     Lv,l,Dpipe,Dv,lv1,...
    %     lv2,Dbias,...
    %     sectionL1,sectionL2,...
    %     'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,...
    %     'meanFlowVelocity',opt.meanFlowVelocity,'isUseStaightPipe',1,...
    %     'm',opt.mach,'notMach',opt.notMach...
    %     ,'isOpening',isOpening...
    %     );%,'coeffDamping',opt.coeffDamping
    % plus1{i} = calcPuls(pressure1,dcpss);
    % plus2{i} = calcPuls(pressure2,dcpss);
    % plus{i} = [plus1{i},plus2{i}];
    
    % calcDatas{dataCount,1} = X;
    % calcDatas{dataCount,2} = plus{i};
    % calcDatas{dataCount,3} = sprintf('进出都错位缓冲罐-lv1:%g,lv2:%g',inputData(i,7),inputData(i,8));%具体描述可以在这改
    % dataCount = dataCount + 1;

%     %计算入口顺接出口偏置
%     [pressure1OSB,pressure2OSB] = ...
%         vesselStraightBiasPulsationCalc(massFlowE,Fre,time,...
%         L1,L2,...
%         Lv,l,Dpipe,Dv,...
%         lv2,Dbias,...
%         sectionL1,sectionL2,...
%         'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,...
%         'meanFlowVelocity',opt.meanFlowVelocity,'isUseStaightPipe',1,...
%         'm',opt.mach,'notMach',opt.notMach...
%         ,'isOpening',isOpening...
%         );%,'coeffDamping',opt.coeffDamping
%     plus1OSB{i} = calcPuls(pressure1OSB,dcpss);
%     plus2OSB{i} = calcPuls(pressure2OSB,dcpss);
%     plusOSB{i} = [plus1OSB{i},plus2OSB{i}];
% 
%     calcDatas{dataCount,1} = X;
%     calcDatas{dataCount,2} = plusOSB{i};
%     calcDatas{dataCount,3} = sprintf('进顺出错缓冲罐-lv1:%g,lv2:%g',inputData(i,7),inputData(i,8));%具体描述可以在这改
%     dataCount = dataCount + 1;
    
    %计算入口偏置出口错位
    [pressure1OBS,pressure2OBS] = ...
        vesselBiasStraightPulsationCalc(massFlowE,Fre,time,...
        L1,L2,...
        Lv,l,Dpipe,Dv,...
        lv2,Dbias,...
        sectionL1,sectionL2,...
        'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,...
        'meanFlowVelocity',opt.meanFlowVelocity,'isUseStaightPipe',1,...
        'm',opt.mach,'notMach',opt.notMach...
        ,'isOpening',isOpening...
        ,'outletpressure',outletPressure(i)...
        );%,'coeffDamping',opt.coeffDamping
    plus1OBS{i} = calcPuls(pressure1OBS,dcpss);
    plus2OBS{i} = calcPuls(pressure2OBS,dcpss);
    plusOBS{i} = [plus1OBS{i},plus2OBS{i}];

    calcDatas{dataCount,1} = X;
    calcDatas{dataCount,2} = plusOBS{i};
    calcDatas{dataCount,3} = sprintf('出顺进错缓冲罐-p:%g',outletPressure(i));%具体描述可以在这改
    dataCount = dataCount + 1;
    % %计算单一缓冲罐
    % if i == 1
    %     [pressure1OV,pressure2OV] = oneVesselPulsationCalc(massFlowE,Fre,time,...
    %         L1,L2,...
    %         Lv,l,Dpipe,Dv,...
    %         sectionL1,sectionL2,...
    %         'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,...
    %         'meanFlowVelocity',opt.meanFlowVelocity,'isUseStaightPipe',1,...
    %         'm',opt.mach,'notMach',opt.notMach...
    %         ,'isOpening',isOpening...
    %         );
    %     plus1OV = calcPuls(pressure1OV,dcpss);
    %     plus2OV = calcPuls(pressure2OV,dcpss);
    %     plusOV = [plus1OV,plus2OV];

    %     calcDatas{dataCount,1} = X;
    %     calcDatas{dataCount,2} = plusOV;
    %     calcDatas{dataCount,3} = sprintf('顺接缓冲罐');%具体描述可以在这改
    %     dataCount = dataCount + 1;
    % end
    
    %计算脉动抑制率
    % temp = plusStraight;
    % temp2 = plus{i};

    % temp(temp<1e-4) = 1;
    % temp2(temp<1e-4) = 1;%temp小于1e-4时，temp2也设置为1.
    % reduceRate{i} = (temp - temp2)./temp;

    % if isempty(plus1{i})
    %     maxPlus1(i) = nan;
    % else
    %     maxPlus1(i) = max(plus1{i});
    % end

    % if isempty(plus2{i})
    %     maxPlus2(i) = nan;
    % else
    %     maxPlus2(i) = max(plus2{i});
    % end  

end

figure
handleCur = plotDataCells(calcDatas);