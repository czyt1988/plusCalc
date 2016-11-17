%% 错位缓冲罐结构与单一顺接缓冲罐对比
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
temp = [1:20];%(Fre<29) ;%| (Fre>30 & Fre < 100);

massFlowE4Vessel = massFlowE;
massFlowE = massFlowE(temp);
Fre4Vessel = Fre;
Fre = Fre(temp);
isDamping = 1;
%绘图参数
isXShowRealLength = 1;
isShowStraightPipe=1;%是否显示直管
isShowOnlyVessel=1;%是否显示无内件缓冲罐
inputData = [...
    %1    2  3     4   5   6     7          8      9     
    %L1  ,L2,Lv   ,l   ,Dpipe,Dv ,lv1    ,lv2   ,Dbias
    1.5 ,6,1.1,0.01,0.098,0.372,0.55-0.232   ,0.55-0.232  , 0 
    1.5 ,6,1.1,0.01,0.098,0.372,0.45  ,0.45 , 0
];%LinLout不能为零
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
for i = 1:size(inputData,1)
    name{i} = sprintf('lv1:%g,lv2:%g',inputData(i,7),inputData(i,8));
%     desp{i} = sprintf('L1:%g,L2:%g,Dpipe:%g,Dv:%g,l:%g,Lv:%g,Linner:%g,Lin:%g，Lout:%g，Din:%g'...
%         ,inputData(i,1),inputData(i,2),inputData(i,3),inputData(i,4)...
%         ,inputData(i,5),inputData(i,6),inputData(i,7),inputData(i,8)...
%         ,inputData(i,9),inputData(i,10));
    
    para(i).opt = opt;
    para(i).L1 = inputData(i,1);%L1(m)
    para(i).L2 = inputData(i,2);%缓冲罐中间连接管道的长度（m） 
    para(i).Lv = inputData(i,3);%管道直径（m）
    para(i).l = inputData(i,4);
    para(i).Dpipe = inputData(i,5);%0.115;%缓冲罐前管道的长度(m)   
    para(i).Dv = inputData(i,6);%[[0.157,0.25,0.5,0.75],[1:0.25:5]];%第一个缓冲罐的直径（m）
    para(i).lv1 = inputData(i,7);%
    para(i).lv2 = inputData(i,8);
    para(i).Dbias = inputData(i,9);%

    
    para(i).sectionL1 = 0:0.2:para(i).L1;
    para(i).sectionL2 = 0:0.2:para(i).L2;
end

dcpss = getDefaultCalcPulsSetStruct();
dcpss.calcSection = [0.3,0.7];
dcpss.fs = Fs;
dcpss.isHp = 0;
dcpss.f_pass = 7;%通过频率5Hz
dcpss.f_stop = 5;%截止频率3Hz
dcpss.rp = 0.1;%边带区衰减DB数设置
dcpss.rs = 30;%截止区衰减DB数设置

dataCount = 1;
for i = 1:length(para)
    %计算直管的长度
    straightPipeLength = para(i).L1 + 2*para(i).l+para(i).Lv + para(i).L2;
    %计算直管的分段
    straightPipeSection = [para(i).sectionL1,...
                            para(i).L1 + 2*para(i).l+para(i).Lv + para(i).sectionL2];
    if isXShowRealLength
        X = straightPipeSection;
    else
        X = 1:length(Y);
    end
    %先对直管进行计算
    if i==1
        %计算直管
        %直管总长
        
        newSectionL2 = para(i).L1 + 2*para(i).l+para(i).Lv + para(i).sectionL2;
        temp = find(straightPipeLength>para(i).L1);%找到缓冲罐所在的索引
        sepratorIndex = temp(1);
        temp = straightPipePulsationCalc(massFlowE,Fre,time,straightPipeLength,straightPipeSection...
        ,'d',para(i).Dpipe,'a',opt.acousticVelocity,'isDamping',opt.isDamping...
        ,'friction',opt.coeffFriction,'meanFlowVelocity',opt.meanFlowVelocity...
        ,'m',para(i).opt.mach,'notMach',para(i).opt.notMach...
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

    [pressure1,pressure2] = ...
        vesselBiasPulsationCalc(massFlowE,Fre,time,...
        para(i).L1,para(i).L2,...
        para(i).Lv,para(i).l,para(i).Dpipe,para(i).Dv,para(i).lv1,...
        para(i).lv2,para(i).Dbias,...
        para(i).sectionL1,para(i).sectionL2,...
        'a',para(i).opt.acousticVelocity,'isDamping',para(i).opt.isDamping,'friction',para(i).opt.coeffFriction,...
        'meanFlowVelocity',para(i).opt.meanFlowVelocity,'isUseStaightPipe',1,...
        'm',para(i).opt.mach,'notMach',para(i).opt.notMach...
        ,'isOpening',isOpening...
        );%,'coeffDamping',opt.coeffDamping
    plus1{i} = calcPuls(pressure1,dcpss);
    plus2{i} = calcPuls(pressure2,dcpss);
    plus{i} = [plus1{i},plus2{i}];
    
    calcDatas{dataCount,1} = X;
    calcDatas{dataCount,2} = plus{i};
    calcDatas{dataCount,3} = sprintf('进出都错位缓冲罐-lv1:%g,lv2:%g',inputData(i,7),inputData(i,8));%具体描述可以在这改
    dataCount = dataCount + 1;

    %计算入口顺接出口偏置
    [pressure1OSB,pressure2OSB] = ...
        vesselStraightBiasPulsationCalc(massFlowE,Fre,time,...
        para(i).L1,para(i).L2,...
        para(i).Lv,para(i).l,para(i).Dpipe,para(i).Dv,...
        para(i).lv2,para(i).Dbias,...
        para(i).sectionL1,para(i).sectionL2,...
        'a',para(i).opt.acousticVelocity,'isDamping',para(i).opt.isDamping,'friction',para(i).opt.coeffFriction,...
        'meanFlowVelocity',para(i).opt.meanFlowVelocity,'isUseStaightPipe',1,...
        'm',para(i).opt.mach,'notMach',para(i).opt.notMach...
        ,'isOpening',isOpening...
        );%,'coeffDamping',opt.coeffDamping
    plus1OSB{i} = calcPuls(pressure1OSB,dcpss);
    plus2OSB{i} = calcPuls(pressure2OSB,dcpss);
    plusOSB{i} = [plus1OSB{i},plus2OSB{i}];

    calcDatas{dataCount,1} = X;
    calcDatas{dataCount,2} = plusOSB{i};
    calcDatas{dataCount,3} = sprintf('进顺出错缓冲罐-lv1:%g,lv2:%g',inputData(i,7),inputData(i,8));%具体描述可以在这改
    dataCount = dataCount + 1;
    
    %计算入口偏置出口错位
    [pressure1OBS,pressure2OBS] = ...
        vesselBiasStraightPulsationCalc(massFlowE,Fre,time,...
        para(i).L1,para(i).L2,...
        para(i).Lv,para(i).l,para(i).Dpipe,para(i).Dv,...
        para(i).lv2,para(i).Dbias,...
        para(i).sectionL1,para(i).sectionL2,...
        'a',para(i).opt.acousticVelocity,'isDamping',para(i).opt.isDamping,'friction',para(i).opt.coeffFriction,...
        'meanFlowVelocity',para(i).opt.meanFlowVelocity,'isUseStaightPipe',1,...
        'm',para(i).opt.mach,'notMach',para(i).opt.notMach...
        ,'isOpening',isOpening...
        );%,'coeffDamping',opt.coeffDamping
    plus1OBS{i} = calcPuls(pressure1OBS,dcpss);
    plus2OBS{i} = calcPuls(pressure2OBS,dcpss);
    plusOBS{i} = [plus1OBS{i},plus2OBS{i}];

    calcDatas{dataCount,1} = X;
    calcDatas{dataCount,2} = plusOBS{i};
    calcDatas{dataCount,3} = sprintf('出顺进错缓冲罐-lv1:%g,lv2:%g',inputData(i,7),inputData(i,8));%具体描述可以在这改
    dataCount = dataCount + 1;
    %计算单一缓冲罐
    if i == 1
        [pressure1OV,pressure2OV] = oneVesselPulsationCalc(massFlowE,Fre,time,...
            para(i).L1,para(i).L2,...
            para(i).Lv,para(i).l,para(i).Dpipe,para(i).Dv,...
            para(i).sectionL1,para(i).sectionL2,...
            'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,...
            'meanFlowVelocity',opt.meanFlowVelocity,'isUseStaightPipe',1,...
            'm',para(i).opt.mach,'notMach',para(i).opt.notMach...
            ,'isOpening',isOpening...
            );
        plus1OV = calcPuls(pressure1OV,dcpss);
        plus2OV = calcPuls(pressure2OV,dcpss);
        plusOV = [plus1OV,plus2OV];

        calcDatas{dataCount,1} = X;
        calcDatas{dataCount,2} = plusOV;
        calcDatas{dataCount,3} = sprintf('顺接缓冲罐');%具体描述可以在这改
        dataCount = dataCount + 1;
    end
    
    %计算脉动抑制率
    temp = plusStraight;
    temp2 = plus{i};

    temp(temp<1e-4) = 1;
    temp2(temp<1e-4) = 1;%temp小于1e-4时，temp2也设置为1.
    reduceRate{i} = (temp - temp2)./temp;

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

end

figure
handleCur = plotDataCells(calcDatas);