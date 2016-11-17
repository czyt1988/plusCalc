%% 入口连接孔管加孔板组合结构与入口连接孔管，孔板，单一缓冲罐对比
clc;
close all;
clear;
currentPath = fileparts(mfilename('fullpath'));
%缓冲罐内入口连接孔管，加孔板结构?
%      L1     l        Lv       l    L2  
%              _________________        
%             |   | dp(n)    |lo|
%             |_ _|_ _ lc    |  |      
%  -----------|_ _ _ _   Din    |----------
%             |Lin|Lout      |  |
%             |___|__________|__|       
%    Dpipe            Dv           Dpipe 
%              
%
% Lin 内插孔管入口段长度 亥姆霍兹共鸣器直径
% Lout内插孔管出口段长度
% lc 孔管壁厚
% dp 孔管每一个孔孔径
% n  开孔个数
% Dp 等效亥姆霍兹共鸣器小径dp*n
% V  亥姆霍兹共鸣器体积
% lv 亥姆霍兹共鸣器长度
% lo 孔板距离缓冲罐出口的距离
%
%       L1  l    Lv          l    L2  
%                   _____________        
%                  | dp(n)    |lo|
%                  |_ _ lc    |  |      
%  ----------------|_ _   Din    |----------
%          lc_| |_ |Lout      |  |
%           |     ||__________|__|       
%           |  Dp |
%        lv |  V  |
%           |     |
%           |_____|
%             Lin
%  Dpipe           Dv               Dpipe
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
temp = (Fre<29) | (Fre>30 & Fre < 100);
Fre = Fre(temp);
massFlowE = massFlowE(temp);
isDamping = 1;
%绘图参数
isXShowRealLength = 1;
isShowStraightPipe=1;%是否显示直管
isShowOnlyVessel=1;%是否显示无内件缓冲罐
inputData = [...
    %1    2  3     4   5    6     7       8      9     10    11    12    13
    %14    15
    %L1 ,L2,Dpipe,Dv , l  , Lv   ,lc    ,lv   ,  dp ,   n,   Lin， Lout, V ,     Din,   lo  
    3.35,10,0.106,0.45,0.01,1.18 ,0.005 ,0.591  ,0.02,  28,  0.2,  0.2,  0.0938, 0.106, 0.3  
    3.35,10,0.106,0.45,0.01,1.18 ,0.005 ,0.591  ,0.0265,16,  0.2,  0.2,  0.0938, 0.106, 0.3
    3.35,10,0.106,0.45,0.01,1.18 ,0.005 ,0.591  ,0.013, 68,  0.25, 0.25, 0.0938, 0.106/2,0.3
];%LrLc不能为零

opt.frequency = 10;%脉动频率
opt.acousticVelocity = 345;%声速
opt.isDamping = isDamping;%是否计算阻尼
opt.coeffDamping = 0.1;%阻尼
opt.coeffFriction = 0.03;%管道摩察系数
opt.meanFlowVelocity = 14.6;%管道平均流速
opt.isUseStaightPipe = 0;%计算容器传递矩阵的方法
opt.mach = opt.meanFlowVelocity / opt.acousticVelocity;
opt.notMach = 1;
for i = 1:size(inputData,1)
    name{i} = sprintf('dp:%g,n:%g,Lin:%g,Din:%g',inputData(i,9),inputData(i,10),inputData(i,11),inputData(i,14));
    % desp{i} = sprintf('L1:%g,L2:%g,Dpipe:%g,Dv:%g,l:%g,Lv:%g,Linner:%g,Lr:%g,Dr1:%g,Dr2:%g,Din:%g'...
    %     ,inputData(i,1),inputData(i,2),inputData(i,3),inputData(i,4)...
    %     ,inputData(i,5),inputData(i,6),inputData(i,7),inputData(i,8)...
    %     ,inputData(i,9),inputData(i,10),inputData(i,11));
    
    para(i).opt = opt;
    para(i).L1 = inputData(i,1);%L1(m)
    para(i).L2 = inputData(i,2);%缓冲罐中间连接管道的长度（m） 
    para(i).Dpipe = inputData(i,3);%管道直径（m）
    para(i).Dv = inputData(i,4);
    para(i).l = inputData(i,5);%0.115;%缓冲罐前管道的长度(m)   
    para(i).Lv = inputData(i,6);%[[0.157,0.25,0.5,0.75],[1:0.25:5]];%第一个缓冲罐的直径（m）
    para(i).lc = inputData(i,7);%
    para(i).lv = inputData(i,8);%
    para(i).dp = inputData(i,9);%
    para(i).n = inputData(i,10);%
    para(i).Lin = inputData(i,11);%
    para(i).Lout = inputData(i,12);
    para(i).V = inputData(i,13);
    para(i).Din = inputData(i,14);
    para(i).lo = inputData(i,15);

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

for i = 1:length(para)
    pressure1 = [];
    pressure2 = [];

    [pressure1,pressure2] = ...
        inPerfPipeInletOrifVesselPulsationCalc(massFlowE,Fre,time,...
        para(i).L1,para(i).L2,...
        para(i).Dpipe,para(i).Dv,para(i).l,para(i).Lv,...
        para(i).lc,para(i).lv,para(i).dp,para(i).n,...
        para(i).Lin,para(i).Lout,para(i).V,para(i).Din,para(i).lo,...
        para(i).sectionL1,para(i).sectionL2,...
        'a',para(i).opt.acousticVelocity,'isDamping',para(i).opt.isDamping,'friction',para(i).opt.coeffFriction,...
        'meanFlowVelocity',para(i).opt.meanFlowVelocity,'isUseStaightPipe',1,...
        'm',para(i).opt.mach,'notMach',para(i).opt.notMach);%,'coeffDamping',opt.coeffDamping
    plus1{i} = calcPuls(pressure1,dcpss);
    plus2{i} = calcPuls(pressure2,dcpss);
    plus{i} = [plus1{i},plus2{i}];
    
    %计算入口接孔管缓冲罐
    if i == 1
        [pressure1OrificeP,pressure2OrificeP] = ...
        innerPerforatedPipeInletVesselPulsationCalc(massFlowE,Fre,time,...
        para(i).L1,para(i).L2,...
        para(i).Dpipe,para(i).Dv,para(i).l,para(i).Lv,...
        para(i).lc,para(i).lv,para(i).dp,para(i).n,...
        para(i).Lin,para(i).Lout,para(i).V,para(i).Din,...
        para(i).sectionL1,para(i).sectionL2,...
        'a',para(i).opt.acousticVelocity,'isDamping',para(i).opt.isDamping,'friction',para(i).opt.coeffFriction,...
        'meanFlowVelocity',para(i).opt.meanFlowVelocity,'isUseStaightPipe',1,...
        'm',para(i).opt.mach,'notMach',para(i).opt.notMach);%,'coeffDamping',opt.coeffDamping
        plus1Orifice = calcPuls(pressure1OrificeP,dcpss);
        plus2Orifice = calcPuls(pressure2OrificeP,dcpss);
        plusOrifice = [plus1Orifice,plus2Orifice];
    end
    if i == 2
        [pressure1OrificeP2,pressure2OrificeP2] = ...
        innerPerforatedPipeInletVesselPulsationCalc(massFlowE,Fre,time,...
        para(i).L1,para(i).L2,...
        para(i).Dpipe,para(i).Dv,para(i).l,para(i).Lv,...
        para(i).lc,para(i).lv,para(i).dp,para(i).n,...
        para(i).Lin,para(i).Lout,para(i).V,para(i).Din,...
        para(i).sectionL1,para(i).sectionL2,...
        'a',para(i).opt.acousticVelocity,'isDamping',para(i).opt.isDamping,'friction',para(i).opt.coeffFriction,...
        'meanFlowVelocity',para(i).opt.meanFlowVelocity,'isUseStaightPipe',1,...
        'm',para(i).opt.mach,'notMach',para(i).opt.notMach);%,'coeffDamping',opt.coeffDamping
        plus1Orifice2 = calcPuls(pressure1OrificeP2,dcpss);
        plus2Orifice2 = calcPuls(pressure2OrificeP2,dcpss);
        plusOrifice2 = [plus1Orifice2,plus2Orifice2];
    end
    
    if i == 3
        [pressure1OrificeP3,pressure2OrificeP3] = ...
        innerPerforatedPipeInletVesselPulsationCalc(massFlowE,Fre,time,...
        para(i).L1,para(i).L2,...
        para(i).Dpipe,para(i).Dv,para(i).l,para(i).Lv,...
        para(i).lc,para(i).lv,para(i).dp,para(i).n,...
        para(i).Lin,para(i).Lout,para(i).V,para(i).Din,...
        para(i).sectionL1,para(i).sectionL2,...
        'a',para(i).opt.acousticVelocity,'isDamping',para(i).opt.isDamping,'friction',para(i).opt.coeffFriction,...
        'meanFlowVelocity',para(i).opt.meanFlowVelocity,'isUseStaightPipe',1,...
        'm',para(i).opt.mach,'notMach',para(i).opt.notMach);%,'coeffDamping',opt.coeffDamping
        plus1Orifice3 = calcPuls(pressure1OrificeP3,dcpss);
        plus2Orifice3 = calcPuls(pressure2OrificeP3,dcpss);
        plusOrifice3 = [plus1Orifice3,plus2Orifice3];
    end
%    %计算带孔板缓冲罐
%     if i == 1
%         [pressure1OrificeP,pressure2OrificeP] = ...
%         vesselHaveOrificePulsationCalc(massFlowE,Fre,time,...
%                 para(i).L1,para(i).L2,...
%                 para(i).Lin,para(i).Lv-para(i).Lin,...
%                 para(i).l,para(i).Dpipe,para(i).Dv,para(i).Dv,para(i).Din,...
%                 para(i).sectionL1,para(i).sectionL2,...
%                 'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,...
%                 'meanFlowVelocity',opt.meanFlowVelocity,'isUseStaightPipe',1,...
%                 'm',para(i).opt.mach,'notMach',para(i).opt.notMach);
%         plus1Orifice = calcPuls(pressure1OrificeP,dcpss);
%         plus2Orifice = calcPuls(pressure2OrificeP,dcpss);
%         plusOrifice = [plus1Orifice,plus2Orifice];
%     end
%     if i == 3
%         [pressure1OrificeP2,pressure2OrificeP2] = ...
%         vesselHaveOrificePulsationCalc(massFlowE,Fre,time,...
%                 para(i).L1,para(i).L2,...
%                 para(i).Lin,para(i).Lv-para(i).Lin,...
%                 para(i).l,para(i).Dpipe,para(i).Dv,para(i).Dv,para(i).Din,...
%                 para(i).sectionL1,para(i).sectionL2,...
%                 'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,...
%                 'meanFlowVelocity',opt.meanFlowVelocity,'isUseStaightPipe',1,...
%                 'm',para(i).opt.mach,'notMach',para(i).opt.notMach);
%         plus1Orifice2 = calcPuls(pressure1OrificeP2,dcpss);
%         plus2Orifice2 = calcPuls(pressure2OrificeP2,dcpss);
%         plusOrifice2 = [plus1Orifice2,plus2Orifice2];
%     end

    %计算单一缓冲罐
    if i == 1
        [pressure1OV,pressure2OV] = oneVesselPulsationCalc(massFlowE,Fre,time,...
            para(i).L1,para(i).L2,...
            para(i).Lv,para(i).l,para(i).Dpipe,para(i).Dv,...
            para(i).sectionL1,para(i).sectionL2,...
            'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,...
            'meanFlowVelocity',opt.meanFlowVelocity,'isUseStaightPipe',1,...
            'm',para(i).opt.mach,'notMach',para(i).opt.notMach);
        plus1OV = calcPuls(pressure1OV,dcpss);
        plus2OV = calcPuls(pressure2OV,dcpss);
        plusOV = [plus1OV,plus2OV];
    end
    if i==1
        %计算直管
        %直管总长
        straightPipeLength = para(i).L1 + 2*para(i).l+para(i).Lv + para(i).L2;
        straightPipeSection = [para(i).sectionL1,...
                                para(i).L1 + 2*para(i).l+para(i).Lv + para(i).sectionL2];
        newSectionL2 = para(i).L1 + 2*para(i).l+para(i).Lv + para(i).sectionL2;
        temp = find(straightPipeLength>para(i).L1);%找到缓冲罐所在的索引
        sepratorIndex = temp(1);
        temp = straightPipePulsationCalc(massFlowE,Fre,time,straightPipeLength,straightPipeSection...
        ,'d',para(i).Dpipe,'a',opt.acousticVelocity,'isDamping',opt.isDamping...
        ,'friction',opt.coeffFriction,'meanFlowVelocity',opt.meanFlowVelocity...
        ,'m',para(i).opt.mach,'notMach',para(i).opt.notMach...
        );
        plusStraight = calcPuls(temp,dcpss);
        maxPlus1Straight(i) = max(plusStraight(1:sepratorIndex(i)));
        maxPlus2Straight(i) = max(plusStraight(sepratorIndex(i):end));
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
plotCount = 1;

for i = 1:length(para)
    marker_index = marker_index + 1;
    if marker_index > market_style_length
        marker_index = 1;
    end
    color_index = color_index + 1;
    if color_index > color_style_length
        color_index = 1;
    end
    Y = plus{i};
    Y = Y./1000;
    Y=Y';
    if isXShowRealLength
        X = straightPipeSection;
    else
        X = 1:length(Y);
    end
    if i==1 
        h(plotCount) = plot(X,plusStraight./1000,'-r','LineWidth',1.5);
        textLegend{plotCount} = '直管';
        plotCount = plotCount + 1;
    end
    if i==1
        %显示缓冲罐
        h(plotCount) = plot(X,plusOV./1000,'-b','LineWidth',1.5);
        textLegend{plotCount} = '单一缓冲罐';
        plotCount = plotCount + 1;
    end
%     if i==1
%         %显示孔板
%         h(plotCount) = plot(X,plusOrifice./1000,'+r','LineWidth',1.5);
%         textLegend{plotCount} = '缓冲罐内孔板(D=0.106)';
%         plotCount = plotCount + 1;
%     end
%     if i==3
%         %显示孔板
%         h(plotCount) = plot(X,plusOrifice2./1000,':g','LineWidth',1.5);
%         textLegend{plotCount} = '缓冲罐内孔板(D=0.053)';
%         plotCount = plotCount + 1;
%     end
     if i==1
        %显示入口接孔管
        h(plotCount) = plot(X,plusOrifice./1000,'-y','LineWidth',1.5);
        textLegend{plotCount} = '缓冲罐入口接孔管n=28';
        plotCount = plotCount + 1;
     end
     if i==2
        %显示入口接孔管
        h(plotCount) = plot(X,plusOrifice2./1000,'-m','LineWidth',1.5);
        textLegend{plotCount} = '缓冲罐入口接孔管n=16';
        plotCount = plotCount + 1;
         end
     if i==3
        %显示入口接孔管
        h(plotCount) = plot(X,plusOrifice3./1000,'-c','LineWidth',1.5);
        textLegend{plotCount} = '缓冲罐入口接孔管n=68';
        plotCount = plotCount + 1;
    end
    h(plotCount) = plot(X,Y,marker_style{marker_index},'color',color_style(color_index,:),'LineWidth',1.5);
    textLegend{plotCount} = name{i};
    plotCount = plotCount + 1;
end
legend(h,textLegend,0);
set(gcf,'color','w');