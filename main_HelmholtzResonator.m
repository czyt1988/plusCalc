%% 当内插孔管可等效为亥姆霍兹共鸣器时，不同开孔形式等效对比
clc;
close all;
clear;
currentPath = fileparts(mfilename('fullpath'));
%亥姆霍兹共鸣器传递矩阵
% lv 共鸣器长
% lc 共鸣器连接管长
% Dp 共鸣器连接管直径 dp*n                         
%       __________                
%      |          |                   
%      |    V     | lv
%      |___    ___|     
%          |  | lc        
% _________|dp|__________                  
% _______________________  
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
Lin = 0.25;%亥姆霍兹共鸣器长
Dv = 0.45;%缓冲罐直径
l=0.01;
%以加工的第三种孔管为例输入参数V
inputData = [...
    %1    2  3     4   5    6       7    8  9     10   11   12    
    %L1  ,L2,Dpipe,Dp ,lc , V  ,    lv ，n, dp,   la1, la2, la     
    3.35,10,0.106,0.052,0.1,0.0392 ,0.799,4, 0.013,0.03,0.04,0.06   
];%LinLout不能为零

opt.frequency = 10;%脉动频率
opt.acousticVelocity = 345;%声速
opt.isDamping = isDamping;%是否计算阻尼
opt.coeffDamping = 0.1;%阻尼
opt.coeffFriction = 0.04;%管道摩察系数
opt.meanFlowVelocity = 14.5;%管道平均流速
opt.isUseStaightPipe = 0;%计算容器传递矩阵的方法
opt.mach = opt.meanFlowVelocity / opt.acousticVelocity;
opt.notMach = 1;
for i = 1:size(inputData,1)
    name{i} = sprintf('Dp:%g,lc:%g,V:%g,lv:%g',inputData(i,4),inputData(i,5),inputData(i,6),inputData(i,7));
    
    para(i).opt = opt;
    para(i).L1 = inputData(i,1);%L1(m)
    para(i).L2 = inputData(i,2);%缓冲罐中间连接管道的长度（m） 
    para(i).Dpipe = inputData(i,3);%管道直径（m）
    para(i).Dp = inputData(i,4);
    para(i).lc = inputData(i,5);%0.115;%缓冲罐前管道的长度(m)   
    para(i).V = inputData(i,6);%[[0.157,0.25,0.5,0.75],[1:0.25:5]];%第一个缓冲罐的直径（m）
    para(i).lv = inputData(i,7);
    para(i).n = inputData(i,8);
    para(i).dp = inputData(i,9);
    para(i).la1 = inputData(i,10);
    para(i).la2 = inputData(i,11);
    para(i).la = inputData(i,12);
    
    para(i).sectionL1 = 0:0.5:para(i).L1;
    para(i).sectionL2 = 0:0.5:para(i).L2;
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
        helmholtzResonatorCalc(massFlowE,Fre,time...
        ,para(i).L1,para(i).L2,para(i).Dpipe...
        ,para(i).dp,para(i).n,para(i).lc,para(i).V,para(i).lv...
        ,para(i).sectionL1,para(i).sectionL2,...
        'a',para(i).opt.acousticVelocity,'isDamping',para(i).opt.isDamping,'friction',para(i).opt.coeffFriction,...
        'meanFlowVelocity',para(i).opt.meanFlowVelocity,...
        'm',para(i).opt.mach,'notMach',para(i).opt.notMach,'sigma',para(i).dp*0.8);%,'coeffDamping',opt.coeffDamping
    plus1{i} = calcPuls(pressure1,dcpss);
    plus2{i} = calcPuls(pressure2,dcpss);
    plus{i} = [plus1{i},plus2{i}];
    
    %计算一列亥姆霍兹，不等效为n*dp
    if i == 1
        [pressure1OH,pressure2OH] = helmholtzResonator_nInParallelCalc(massFlowE,Fre,time,...
            para(i).L1,para(i).L2,...
            para(i).Dpipe,para(i).V,para(i).lv,para(i).lc,...
            para(i).dp,para(i).n,para(i).la1,para(i).la2,...
            para(i).sectionL1,para(i).sectionL2,...
            'a',para(i).opt.acousticVelocity,'isDamping',para(i).opt.isDamping,'friction',para(i).opt.coeffFriction,...
            'meanFlowVelocity',para(i).opt.meanFlowVelocity,...
            'm',para(i).opt.mach,'notMach',para(i).opt.notMach,'sigma',para(i).dp*0.8);%,'coeffDamping',opt.coeffDamping
        plus1OH = calcPuls(pressure1OH,dcpss);
        plus2OH = calcPuls(pressure2OH,dcpss);
        plusOH = [plus1OH,plus2OH];
    end
    
    %计算一行亥姆霍兹，不等效为n*dp
    if i == 1
        [pressure1H,pressure2H] = helmholtzResonator_nInSeriesCalc(massFlowE,Fre,time,...
            para(i).L1,para(i).L2,...
            para(i).Dpipe,para(i).V,para(i).lv,para(i).lc,...
            para(i).dp,para(i).la1,para(i).la2,para(i).la,...
            para(i).sectionL1,para(i).sectionL2,...
             'a',para(i).opt.acousticVelocity,'isDamping',para(i).opt.isDamping,'friction',para(i).opt.coeffFriction,...
           'meanFlowVelocity',para(i).opt.meanFlowVelocity,...
           'm',para(i).opt.mach,'notMach',para(i).opt.notMach,'sigma',para(i).dp*0.8);%,'coeffDamping',opt.coeffDamping
        plus1H = calcPuls(pressure1H,dcpss);
        plus2H = calcPuls(pressure2H,dcpss);
        plusH = [plus1H,plus2H];
    end
    
    %计算单一缓冲罐
    if i == 1
        [pressure1OV,pressure2OV] = oneVesselPulsationCalc(massFlowE,Fre,time,...
            para(i).L1,para(i).L2,...
            Lv,l,para(i).Dpipe,Dv,...
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
        straightPipeLength = para(i).L1 + para(i).L2;
        straightPipeSection = [para(i).sectionL1,...
                                 para(i).L1+para(i).sectionL2];
       
        temp = find(straightPipeLength>para(i).L1);%找到缓冲罐所在的索引
        sepratorIndex = temp(1);
        temp = straightPipePulsationCalc(massFlowE,Fre,time,straightPipeLength,straightPipeSection...
        ,'d',para(i).Dpipe,'a',opt.acousticVelocity,'isDamping',opt.isDamping...
        ,'friction',opt.coeffFriction,'meanFlowVelocity',opt.meanFlowVelocity...
        ,'m',para(i).opt.mach,'notMach',para(i).opt.notMach...
        );
        plusStraight = calcPuls(temp,dcpss);
        maxPlus1Straight(i) = max(plusStraight(1:sepratorIndex));
        maxPlus2Straight(i) = max(plusStraight(sepratorIndex:end));
   
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
    if isShowStraightPipe 
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
    if i==1
        %显示一列孔亥姆霍兹共鸣器
        h(plotCount) = plot(X,plusOH/1000,'-m','LineWidth',1.5);
        textLegend{plotCount} = '孔管开一圈4个孔';
        plotCount = plotCount + 1;
    end
        if i==1
        %显示一排孔亥姆霍兹共鸣器
        h(plotCount) = plot(X,plusH./1000,'-y','LineWidth',1.5);
        textLegend{plotCount} = '孔管开一排4个孔';
        plotCount = plotCount + 1;
    end
    h(plotCount) = plot(X,Y,marker_style{marker_index},'color',color_style(color_index,:),'LineWidth',1.5);
    textLegend{plotCount} = name{i};
    plotCount = plotCount + 1;
    
end
legend(h,textLegend,0);
set(gcf,'color','w');