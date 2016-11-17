%% 缓冲罐内置孔管结构与单一顺接缓冲罐对比
%内置孔管的入口部分环绕一周开孔，等效为亥姆霍兹共鸣器；出口部分不开孔，仅为内插管
clc;
close all;
clear;
currentPath = fileparts(mfilename('fullpath'));
%缓冲罐内入口连接孔管,Lin部分等效为亥姆霍兹共鸣器.Lout部分为内插开口管，不开孔
%      L1     l                 Lv              l    L2  
%              _________________________________        
%             | lc dp(n1) |                     |
%             |___ _ _ ___|___________          |     
%  -----------|___ _ _ ___ ___________ Din      |----------
%             |la1 lp1 la2|                     |
%             |___________|_____________________|       
%                  Lin         Lout
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
% V  亥姆霍兹共鸣器体积 V=(Sv-S)*Lin
% lv 等效共鸣器长 V./((pi*Lin^2)./4)（不一定等效准确）
% xSection1，xSection2 孔管每圈孔的间距，从0开始算，x的长度为孔管孔的圈数+1，x的值是当前一圈孔和上一圈孔的距离，如果间距一样，那么x里的值都一样

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

opt.frequency = 10;%脉动频率
opt.acousticVelocity = 345;%声速
opt.isDamping = isDamping;%是否计算阻尼
opt.coeffDamping = 0.1;%阻尼
opt.coeffFriction = 0.04;%管道摩察系数
opt.meanFlowVelocity = 14.5;%管道平均流速
opt.isUseStaightPipe = 1;%计算容器传递矩阵的方法
opt.mach = opt.meanFlowVelocity / opt.acousticVelocity;
opt.notMach = 0;

variant_n1 = [4];              %variant_n = [6,6];sectionNum1 =[1,6];%对应孔1的组数sectionNum2 =[1,1];%对应孔2的组数
sectionNum1 =[1];%对应孔1的组数
% sectionNum2 =[4];%对应孔2的组数
% variant_n2 = [4];
% variant_lp2 = [0.16];
%  variant_lv = [5,1.5,0.955,0.5,0.25];
% variant_V = [0.01];%等效体积为0.01时和普通孔管矩阵的计算结果较为接近
% variant_dp = [0.02];
variant_lc = [0.005,0.2,0.3,0.4,0.5,0.6];%相当于将孔管的孔改为内插小细管

for i = 1:length(variant_lc)
    
    para(i).opt = opt;
    para(i).L1 = 3.85;%L1(m)
    para(i).L2 = 3;%L2（m）长度
    para(i).Dpipe = 0.106;%管道直径（m）
    para(i).vhpicStruct.l = 0.01;
    para(i).vhpicStruct.Dv = 0.45;%缓冲罐的直径（m）
    para(i).vhpicStruct.Lv = 1.18;%缓冲罐总长
    para(i).vhpicStruct.lc = variant_lc(i);%内插管壁厚
    para(i).vhpicStruct.dp = 0.02;%开孔径
    para(i).vhpicStruct.Lin = 0.2;%内插管入口段长度
    para(i).vhpicStruct.lp1 = 0.04;%内插管入口段长度
%     para(i).vhpicStruct.lp2 = variant_lp2(i);%内插管出口段非孔管开孔长度
    para(i).vhpicStruct.n1 = 4;%入口段孔数
%     para(i).vhpicStruct.n2 = variant_n2(i);%出口段孔数
    para(i).vhpicStruct.la1 = 0.01;%孔管入口段靠近入口长度
    para(i).vhpicStruct.la2 = 0.15;%孔管
%     para(i).vhpicStruct.lb1 = 0.03;
%     para(i).vhpicStruct.lb2 = 0.01;
    para(i).vhpicStruct.Din = 0.106;
    para(i).vhpicStruct.Lout = 0.2;
    para(i).vhpicStruct.V = 0.01;%(pi.*para(i).vhpicStruct.Dv.^2./4-pi.*para(i).vhpicStruct.Din.^2./4)*para(i).vhpicStruct.Lin;  %0.03
    para(i).vhpicStruct.lv = 0.995;
    para(i).vhpicStruct.bp1 = para(i).vhpicStruct.n1.*(para(i).vhpicStruct.dp)^2./(4.*para(i).vhpicStruct.Din.*para(i).vhpicStruct.lp1);%开孔率
    l = para(i).vhpicStruct.lp1;
%     para(i).vhpicStruct.xSection1 = [0,ones(1,sectionNum1(i)).*(l/(sectionNum1(i)))];
    para(i).vhpicStruct.xSection1 = [0,ones(1,4).*(l/(4))];
%     l = para(i).vhpicStruct.lp2;
%     para(i).vhpicStruct.xSection2 = [0,ones(1,sectionNum2(i)).*(l/(sectionNum2(i)))];
%     para(i).vhpicStruct.xSection2
    para(i).sectionL1 = 0:0.5:para(i).L1;
    para(i).sectionL2 = 0:0.5:para(i).L2;
    
    holepipeLength1 = para(i).vhpicStruct.Lin - para(i).vhpicStruct.la1 - para(i).vhpicStruct.la2;
    hl1 = sum(para(i).vhpicStruct.xSection1);
    if(~cmpfloat(holepipeLength1,hl1))
        error('孔管参数设置错误：holepipeLength1=%.8f,hl1=%.8f;Lin:%g,la1:%g,la2:%g,sum(xSection1):%g,dp:%g'...
            ,holepipeLength1,hl1...
            ,para(i).vhpicStruct.Lin,para(i).vhpicStruct.la1,para(i).vhpicStruct.la2...
            ,sum(para(i).vhpicStruct.xSection1),para(i).vhpicStruct.dp);
    end
    
%     holepipeLength2 = para(i).vhpicStruct.Lout - para(i).vhpicStruct.lb1 - para(i).vhpicStruct.lb2;
%     hl2 = sum(para(i).vhpicStruct.xSection2);
%     if(holepipeLength2 ~= hl2)
%         error('孔管参数设置错误：holepipeLength2=%g,hl2=%g;Lout:%g,lb1:%g,lb2:%g,sum(xSection2):%g,dp:%g'...
%             ,holepipeLength2,hl2...
%             ,para(i).vhpicStruct.Lout,para(i).vhpicStruct.lb1,para(i).vhpicStruct.lb2...
%             ,sum(para(i).vhpicStruct.xSection2),para(i).vhpicStruct.dp);
%     end
    name{i} = sprintf('lc:%g',variant_lc(i));
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
        vesselHavePerfInletHelmParaNoneOutCalc(massFlowE,Fre,time,...
        para(i).L1,para(i).L2,para(i).Dpipe...
        ,para(i).vhpicStruct,...
        para(i).sectionL1,para(i).sectionL2,...
        'a',para(i).opt.acousticVelocity,'isDamping',para(i).opt.isDamping,'friction',para(i).opt.coeffFriction,...
        'meanFlowVelocity',para(i).opt.meanFlowVelocity,...
        'm',para(i).opt.mach,'notMach',para(i).opt.notMach);%,'coeffDamping',opt.coeffDamping
    plus1{i} = calcPuls(pressure1,dcpss);
    plus2{i} = calcPuls(pressure2,dcpss);
    plus{i} = [plus1{i},plus2{i}];
    

    %计算单一缓冲罐
    if i == 1
        [pressure1OV,pressure2OV] = oneVesselPulsationCalc(massFlowE,Fre,time,...
            para(i).L1,para(i).L2,...
            para(i).vhpicStruct.Lv,para(i).vhpicStruct.l,para(i).Dpipe,para(i).vhpicStruct.Dv,...
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
        straightPipeLength = para(i).L1 + 2*para(i).vhpicStruct.l+para(i).vhpicStruct.Lv + para(i).L2;
        straightPipeSection = [para(i).sectionL1,...
                                para(i).L1 + 2*para(i).vhpicStruct.l+para(i).vhpicStruct.Lv + para(i).sectionL2];
        newSectionL2 = para(i).L1 + 2*para(i).vhpicStruct.l+para(i).vhpicStruct.Lv + para(i).sectionL2;
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

    h(plotCount) = plot(X,Y,marker_style{marker_index},'color',color_style(color_index,:),'LineWidth',1.5);
    textLegend{plotCount} = name{i};
    plotCount = plotCount + 1;
    
end
legend(h,textLegend,0);
set(gcf,'color','w');

figure
Y = maxPlus1;
Y = Y./1000;
X = 1:length(Y);
hBar = bar(X,Y);
set(gca,'XTickLabel',name);
set(gcf,'color','w');

figure
Y = maxPlus2;
Y = Y./1000;
X = 1:length(Y);
hBar = bar(X,Y);
set(gca,'XTickLabel',name);
set(gcf,'color','w');