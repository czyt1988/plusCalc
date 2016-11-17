%% 缓冲罐内置孔管结构与单一顺接缓冲罐对比
clc;
close all;
clear;
currentPath = fileparts(mfilename('fullpath'));
%缓冲罐中间插入孔管,两端堵死，开孔个数不足以等效为亥姆霍兹共鸣器(缓冲罐入口偏置)
%                 L1
%                     |
%                     |
%           l         |          Lv              l    L2  
%              _______|_________________________        
%             |    dp(n1)    |    dp(n2)        |
%             |   ___ _ _ ___|___ _ _ ___ lc    |     
%             |  |___ _ _ ___ ___ _ _ ___|Din   |----------
%             |   la1 lp1 la2|lb1 lp2 lb2       |
%             |______________|__________________|       
%                  Lin             Lout
%    Dpipe                   Dv                     Dpipe 
%              
%
% Lin 内插孔管入口段长度 
% Lout内插孔管出口段长度
% lc  孔管壁厚
% dp  孔管每一个孔孔径
% n1  孔管入口段开孔个数；    n2  孔管出口段开孔个数
% la1 孔管入口段距入口长度 
% la2 孔管入口段距隔板长度
% lb1 孔管出口段距隔板长度
% lb2 孔管出口段距开孔长度
% lp1 孔管入口段开孔长度
% lp2 孔管出口段开孔长度
% Din 孔管管径；
% xSection1，xSection2 孔管每圈孔的间距，从0开始算，x的长度为孔管孔的圈数+1，x的值是当前一圈孔和上一圈孔的距离，如果间距一样，那么x里的值都一样
isOpening = 0;%管道闭口
%  rpm = 300;outDensity = 1.9167;multFre=[10,20,30];%环境25度绝热压缩到0.2MPaG的温度对应密度
rpm = 420;outDensity = 1.5608;multFre=[14,28,42];%环境25度绝热压缩到0.15MPaG的温度对应密度
Fs = 4096;
[massFlowRaw,time,~,opt.meanFlowVelocity] = massFlowMaker(0.25,0.098,rpm...
	,0.14,1.075,outDensity,'rcv',0.15,'k',1.4,'pr',0.15,'fs',Fs,'oneSecond',6);

%massFlow = load(fullfile(currentPath,'mass_flow_0.1478_NorthZone.txt'));

[FreRaw,AmpRaw,PhRaw,massFlowERaw] = frequencySpectrum(detrend(massFlowRaw,'constant'),Fs);

% 提取主要频率
% massFlowE = massFlowERaw;
% Fre = FreRaw;

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

opt.acousticVelocity = 345;%声速
opt.isDamping = isDamping;%是否计算阻尼
opt.coeffDamping = nan;%阻尼
opt.coeffFriction = 0.04;%管道摩察系数
% opt.meanFlowVelocity =14.5;%14.5;%管道平均流速
opt.isUseStaightPipe = 1;%计算容器传递矩阵的方法
opt.mach = opt.meanFlowVelocity / opt.acousticVelocity;
opt.notMach = 1;


calcDatas = {};

iteData = [0.05:0.01:0.2];

for i = 1:length(iteData)     
    para(i).opt = opt;
    para(i).L1 = 1.5;%L1(m)
    para(i).L2 = 6;%L2（m）长度
    para(i).Dpipe = 0.098;%管道直径（m）
    para(i).vhpicStruct.l = 0.01;
    para(i).vhpicStruct.Dv = 0.372;%缓冲罐的直径（m）
    para(i).vhpicStruct.Lv = 1.1;%缓冲罐总长
    para(i).vhpicStruct.lc = 0.005;%内插管壁厚
    para(i).vhpicStruct.dp = 0.013;%开孔径
    para(i).vhpicStruct.Lin = 0.25;%内插管入口段长度
    para(i).vhpicStruct.lp2 = iteData(i);%内插管入口段非孔管开孔长度
    para(i).vhpicStruct.lp1 = 0.16;%内插管出口段孔管开孔长度
    para(i).vhpicStruct.n1 = 48;%入口段孔数
    para(i).vhpicStruct.n2 = 48;%出口段孔数
    para(i).vhpicStruct.lb2 = 0.03;%孔管入口段靠近入口长度
    para(i).vhpicStruct.lb1 = para(i).vhpicStruct.Lin - para(i).vhpicStruct.lb2 - para(i).vhpicStruct.lp2;%孔管
    para(i).vhpicStruct.la2 = 0.06;
    para(i).vhpicStruct.la1 = 0.03;
    para(i).vhpicStruct.Din = 0.053;
    para(i).vhpicStruct.Lout = 0.25;
    sectionNum1 =[1];%对应孔1的组数
    sectionNum2 =[1];%对应孔2的组数
    para(i).vhpicStruct.bp1 = para(i).vhpicStruct.n1.*(para(i).vhpicStruct.dp)^2./(4.*para(i).vhpicStruct.Din.*para(i).vhpicStruct.lp1);%开孔率
    para(i).vhpicStruct.bp2 = para(i).vhpicStruct.n2.*(para(i).vhpicStruct.dp)^2./(4.*para(i).vhpicStruct.Din.*para(i).vhpicStruct.lp2);%开孔率
    lTmp = para(i).vhpicStruct.lp1;
    para(i).vhpicStruct.xSection1 = [0,ones(1,sectionNum1).*(lTmp/(sectionNum1))];
    lTmp = para(i).vhpicStruct.lp2;
    para(i).vhpicStruct.xSection2 = [0,ones(1,sectionNum2).*(lTmp/(sectionNum2))];
    para(i).sectionL1 = 0:0.25:para(i).L1;
    para(i).sectionL2 = 0:0.25:para(i).L2;
    para(i).vhpicStruct.lv1 = para(i).vhpicStruct.Lv./2-0.232;%232
    para(i).vhpicStruct.lv2 = 0;%出口不偏置


    holepipeLength1 = para(i).vhpicStruct.Lin - para(i).vhpicStruct.la1 - para(i).vhpicStruct.la2;
    hl1 = sum(para(i).vhpicStruct.xSection1);
    if(~cmpfloat(holepipeLength1,hl1))
        error('孔管参数设置错误：holepipeLength1=%.8f,hl1=%.8f;Lin:%g,la1:%g,la2:%g,sum(xSection1):%g,dp:%g'...
            ,holepipeLength1,hl1...
            ,para(i).vhpicStruct.Lin,para(i).vhpicStruct.la1,para(i).vhpicStruct.la2...
            ,sum(para(i).vhpicStruct.xSection1),para(i).vhpicStruct.dp);
    end
    name{i} = sprintf('lp1:%g',iteData(i));
end

dcpss = getDefaultCalcPulsSetStruct();
dcpss.calcSection = [0.3,0.7];
dcpss.fs = Fs;
dcpss.isHp = 0;
dcpss.f_pass = 7;%通过频率5Hz
dcpss.f_stop = 5;%截止频率3Hz
dcpss.rp = 0.1;%边带区衰减DB数设置
dcpss.rs = 30;%截止区衰减DB数设置

dataCount = 2;
calcDatas{1,2} = 'x值';
calcDatas{1,3} = '压力脉动';
calcDatas{1,4} = '1倍频';
calcDatas{1,5} = '2倍频';
calcDatas{1,6} = '3倍频';
for i = 1:length(para)
    
    if i==1
        %计算直管
        %直管总长
        straightPipeLength = para(i).L1 + 2*para(i).vhpicStruct.l+para(i).vhpicStruct.Lv + para(i).L2;
        straightPipeSection = [para(i).sectionL1,...
                                para(i).L1 + 2*para(i).vhpicStruct.l+para(i).vhpicStruct.Lv + para(i).sectionL2];
        newSectionL2 = para(i).L1 + 2*para(i).vhpicStruct.l+para(i).vhpicStruct.Lv + para(i).sectionL2;
        temp = find(straightPipeLength>para(i).L1);%找到缓冲罐所在的索引
        sepratorIndex = temp(1);
        temp = straightPipePulsationCalc(massFlowE,Fre,time,straightPipeLength,straightPipeSection...
        ,'d',para(i).Dpipe,'a',opt.acousticVelocity,'isDamping',opt.isDamping...
        ,'friction',opt.coeffFriction,'meanFlowVelocity',opt.meanFlowVelocity...
        ,'m',para(i).opt.mach,'notMach',para(i).opt.notMach,...
        'isOpening',isOpening);
        plusStraight = calcPuls(temp,dcpss);
        maxPlus1Straight(i) = max(plusStraight(1:sepratorIndex(i)));
        maxPlus2Straight(i) = max(plusStraight(sepratorIndex(i):end));
        multFreAmpValue_straightPipe{i} = calcWaveFreAmplitude(temp,Fs,multFre,'freErr',1);

        if isXShowRealLength
            X = straightPipeSection;
        else
            X = 1:length(plusStraight);
        end
        calcDatas{dataCount,1} = sprintf('直管');
        calcDatas{dataCount,2} = X;
        calcDatas{dataCount,3} = plusStraight;
        calcDatas{dataCount,4} = multFreAmpValue_straightPipe{i}(1,:);
        calcDatas{dataCount,5} = multFreAmpValue_straightPipe{i}(2,:);
        calcDatas{dataCount,6} = multFreAmpValue_straightPipe{i}(3,:);
        dataCount = dataCount + 1;
    end

    
    pressure1ClosedIB = [];
    pressure2ClosedIB = [];

    [pressure1ClosedIB,pressure2ClosedIB] = ...
        vesselInBiasHaveInnerPerfBothClosedCompCalc(massFlowE,Fre,time,...
        para(i).L1,para(i).L2,para(i).Dpipe...
        ,para(i).vhpicStruct,...
        para(i).sectionL1,para(i).sectionL2,...
        'a',para(i).opt.acousticVelocity,'isDamping',para(i).opt.isDamping,'friction',para(i).opt.coeffFriction,...
        'meanFlowVelocity',para(i).opt.meanFlowVelocity,...
        'm',para(i).opt.mach,'notMach',para(i).opt.notMach,...
        'isOpening',isOpening);%,'coeffDamping',para(i).opt.coeffDamping,
    plus1ClosedIB{i} = calcPuls(pressure1ClosedIB,dcpss);
    plus2ClosedIB{i} = calcPuls(pressure2ClosedIB,dcpss);
    plusClosedIB{i} = [plus1ClosedIB{i},plus2ClosedIB{i}];
    multFreAmpValueClosedIB{i} = calcWaveFreAmplitude([pressure1ClosedIB,pressure2ClosedIB],Fs,multFre,'freErr',1);
    
    calcDatas{dataCount,1} = sprintf('入口偏置内插孔管两端堵死缓冲罐,lp1:%g',iteData(i));
    calcDatas{dataCount,2} = X;
    calcDatas{dataCount,3} = plusClosedIB{i};
    calcDatas{dataCount,4} = multFreAmpValueClosedIB{i}(1,:);
    calcDatas{dataCount,5} = multFreAmpValueClosedIB{i}(2,:);
    calcDatas{dataCount,6} = multFreAmpValueClosedIB{i}(3,:);
    dataCount = dataCount + 1;

%     [pressure1OpenIB,pressure2OpenIB] = ...
%         vesselInBiasHaveInnerPerfInCloOutOpenCompCalc(massFlowE,Fre,time,...
%         para(i).L1,para(i).L2,para(i).Dpipe...
%         ,para(i).vhpicStruct,...
%         para(i).sectionL1,para(i).sectionL2,...
%         'a',para(i).opt.acousticVelocity,'isDamping',para(i).opt.isDamping,'friction',para(i).opt.coeffFriction,...
%         'meanFlowVelocity',para(i).opt.meanFlowVelocity,...
%         'm',para(i).opt.mach,'notMach',para(i).opt.notMach...
%         ,'isOpening',isOpening);%,'coeffDamping',opt.coeffDamping
%     plus1OpenIB{i} = calcPuls(pressure1OpenIB,dcpss);
%     plus2OpenIB{i} = calcPuls(pressure2OpenIB,dcpss);
%     plusOpenIB{i} = [plus1OpenIB{i},plus2OpenIB{i}];
%     multFreAmpValueOpendIB{i} = calcWaveFreAmplitude([pressure1OpenIB,pressure2OpenIB],Fs,multFre,'freErr',1);
%     calcDatas{dataCount,1} = sprintf('入口偏置内插孔管开口缓冲罐,lp1:%g',iteData(i));
%     calcDatas{dataCount,2} = X;
%     calcDatas{dataCount,3} = plusOpenIB{i};
%     calcDatas{dataCount,4} = multFreAmpValueOpendIB{i}(1,:);
%     calcDatas{dataCount,5} = multFreAmpValueOpendIB{i}(2,:);
%     calcDatas{dataCount,6} = multFreAmpValueOpendIB{i}(3,:);
%     dataCount = dataCount + 1;

     %计算单一缓冲罐入口偏置
     if i == 1
%         [pressure1Temp,pressure2Temp] = vesselBiasStraightPulsationCalc(massFlowE4Vessel,Fre4Vessel,time,...
%             para(i).L1,para(i).L2,...
%             para(i).vhpicStruct.Lv,para(i).vhpicStruct.l,para(i).Dpipe,para(i).vhpicStruct.Dv,...
%             para(i).vhpicStruct.lv1,0,...%用0是指？？
%             para(i).sectionL1,para(i).sectionL2,...
%             'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,...
%             'meanFlowVelocity',opt.meanFlowVelocity,'isUseStaightPipe',1,...
%             'm',para(i).opt.mach,'notMach',para(i).opt.notMach,...
%             'isOpening',isOpening);
%         plusOVIB = [calcPuls(pressure1Temp,dcpss),calcPuls(pressure2Temp,dcpss)];
%         multFreAmpValue_OVIB{i} = calcWaveFreAmplitude([pressure1Temp,pressure2Temp],Fs,multFre,'freErr',1);
% 
%         calcDatas{dataCount,1} = sprintf('无内件缓冲罐-进口错位出口顺接');
% 	    calcDatas{dataCount,2} = X;
% 	    calcDatas{dataCount,3} = plusOVIB;
% 	    calcDatas{dataCount,4} = multFreAmpValue_OVIB{i}(1,:);
% 	    calcDatas{dataCount,5} = multFreAmpValue_OVIB{i}(2,:);
% 	    calcDatas{dataCount,6} = multFreAmpValue_OVIB{i}(3,:);
%     	dataCount = dataCount + 1;
%         
%         [pressure1Temp,pressure2Temp] = vesselBiasPulsationCalc(massFlowE4Vessel,Fre4Vessel,time,...
%         para(i).L1,para(i).L2,...
%         para(i).vhpicStruct.Lv,para(i).vhpicStruct.l,para(i).Dpipe,para(i).vhpicStruct.Dv...
%         ,para(i).vhpicStruct.lv1,...
%         para(i).vhpicStruct.lv2,0,...
%         para(i).sectionL1,para(i).sectionL2,...
%         'a',para(i).opt.acousticVelocity,'isDamping',para(i).opt.isDamping,'friction',para(i).opt.coeffFriction,...
%         'meanFlowVelocity',para(i).opt.meanFlowVelocity,'isUseStaightPipe',1,...
%         'm',para(i).opt.mach,'notMach',para(i).opt.notMach...
%         ,'isOpening',isOpening...
%         );%,'coeffDamping',opt.coeffDamping
%         plusOVIB = [calcPuls(pressure1Temp,dcpss),calcPuls(pressure2Temp,dcpss)];
%         multFreAmpValue_OVIB{i} = calcWaveFreAmplitude([pressure1Temp,pressure2Temp],Fs,multFre,'freErr',1);
%         calcDatas{dataCount,1} = sprintf('无内件缓冲罐-进出口错位');
% 	    calcDatas{dataCount,2} = X;
% 	    calcDatas{dataCount,3} = plusOVIB;
% 	    calcDatas{dataCount,4} = multFreAmpValue_OVIB{i}(1,:);
% 	    calcDatas{dataCount,5} = multFreAmpValue_OVIB{i}(2,:);
% 	    calcDatas{dataCount,6} = multFreAmpValue_OVIB{i}(3,:);
%     	dataCount = dataCount + 1;
%         
        [pressure1Temp,pressure2Temp] = vesselStraightBiasPulsationCalc(massFlowE,Fre,time...
        ,para(i).L1,para(i).L2...
        ,para(i).vhpicStruct.Lv,para(i).vhpicStruct.l,para(i).Dpipe,para(i).vhpicStruct.Dv...
        ,para(i).vhpicStruct.lv2,0 ...
        ,para(i).sectionL1,para(i).sectionL2,...
        'a',para(i).opt.acousticVelocity,'isDamping',para(i).opt.isDamping,'friction',para(i).opt.coeffFriction,...
        'meanFlowVelocity',para(i).opt.meanFlowVelocity,'isUseStaightPipe',1,...
        'm',para(i).opt.mach,'notMach',para(i).opt.notMach...
        ,'isOpening',isOpening...
        );%,'coeffDamping',opt.coeffDamping
        plus1OVIB = calcPuls(pressure1Temp,dcpss);
        plus2OVIB = calcPuls(pressure2Temp,dcpss);
        plusOVIB = [plus1OVIB,plus2OVIB];
        multFreAmpValue_OVIB{i} = calcWaveFreAmplitude([pressure1Temp,pressure2Temp],Fs,multFre,'freErr',1);
        calcDatas{dataCount,1} = sprintf('无内件缓冲罐-进口顺接出口错位');
	    calcDatas{dataCount,2} = X;
	    calcDatas{dataCount,3} = plusOVIB;
	    calcDatas{dataCount,4} = multFreAmpValue_OVIB{i}(1,:);
	    calcDatas{dataCount,5} = multFreAmpValue_OVIB{i}(2,:);
	    calcDatas{dataCount,6} = multFreAmpValue_OVIB{i}(3,:);
    	dataCount = dataCount + 1;
        if isempty(plus1OVIB)
            maxOVPlus1(i) = nan;
        else
            maxOVPlus1(i) = max(plus1OVIB);
        end

        if isempty(plus2OVIB)
            maxOVPlus2(i) = nan;
        else
            maxOVPlus2(i) = max(plus2OVIB);
        end  
    end
    
    %计算脉动抑制率
%     temp = plusOVIB;
%     temp2 = plusOpenIB{i};
% 
%     temp(temp<1e-4) = 1;
%     temp2(temp<1e-4) = 1;%temp小于1e-4时，temp2也设置为1.
%     reduceRate{i} = (temp - temp2)./temp;

    if isempty(plus1ClosedIB{i})
        maxPlus1(i) = nan;
    else
        maxPlus1(i) = max(plus1ClosedIB{i});
    end

    if isempty(plus2ClosedIB{i})
        maxPlus2(i) = nan;
    else
        maxPlus2(i) = max(plus2ClosedIB{i});
    end  

end

figure
plot(iteData,maxPlus1);
set(gcf,'color','w');

figure
plot(iteData,maxPlus2.*3.9./1000,'LineWidth',1.5);
xlabel('孔管开孔段长度(m)');
ylabel('管系最大压力脉动峰峰值(kPa)');
set(gcf,'color','w');
set(gcf,'position',[100,100,300,250]);
xlim([0.05,0.2]);
ignoreHeader = 1;
%绘制压力脉动
figure 
plotDataCells(calcDatas,'xcol',2,'ycol',3,'legendcol',1,'ignoreHeader',ignoreHeader);
title('脉动压力峰峰值');
% %绘制1倍频
% figure
% plotDataCells(calcDatas,'xcol',2,'ycol',4,'legendcol',1,'ignoreHeader',ignoreHeader);
% title('压力1倍频');
% %绘制2倍频
% figure
% plotDataCells(calcDatas,'xcol',2,'ycol',5,'legendcol',1,'ignoreHeader',ignoreHeader);
% title('压力2倍频');
% %绘制3倍频
% figure
% plotDataCells(calcDatas,'xcol',2,'ycol',6,'legendcol',1,'ignoreHeader',ignoreHeader);
% title('压力3倍频');
