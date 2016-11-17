%% 缓冲罐内置孔管结构与单一顺接缓冲罐对比
clc;
close all;
clear;
currentPath = fileparts(mfilename('fullpath'));
%缓冲罐中间插入孔管,两端堵死，开孔个数不足以等效为亥姆霍兹共鸣器
%      L1     l                 Lv              l    L2  
%              _________________________________        
%             |    dp(n1)    |    dp(n2)        |
%             |   ___ _ _ ___|___ _ _ ___ lc    |     
%  -----------|  |___ _ _ ___ ___ _ _ ___|Din   |----------
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
% rpm = 300;outDensity = 1.9167;multFre=[10,20,30];%环境25度绝热压缩到0.2MPaG的温度对应密度
rpm = 420;outDensity = 1.5608;multFre=[14,28,42];%环境25度绝热压缩到0.15MPaG的温度对应密度
Fs = 4096;
[massFlowRaw,time,~,opt.meanFlowVelocity] = massFlowMaker(0.25,0.098,rpm...
	,0.14,1.075,outDensity,'rcv',0.15,'k',1.4,'pr',0.15,'fs',Fs,'oneSecond',1);

%massFlow = load(fullfile(currentPath,'mass_flow_0.1478_NorthZone.txt'));

[FreRaw,AmpRaw,PhRaw,massFlowE] = frequencySpectrum(detrend(massFlowRaw,'constant'),Fs);
Fre = FreRaw;
% 提取主要频率
[pks,locs] = findpeaks(AmpRaw,'SORTSTR','descend');
Fre = FreRaw(locs);
massFlowE = massFlowE(locs);
temp = [1:20];%(Fre<29) ;%| (Fre>30 & Fre < 100);
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
% opt.meanFlowVelocity =14.5;%14.5;%管道平均流速
opt.isUseStaightPipe = 1;%计算容器传递矩阵的方法
opt.mach = opt.meanFlowVelocity / opt.acousticVelocity;
opt.notMach = 0;

variant_n1 = [5:90];              %variant_n = [6,6];sectionNum1 =[1,6];%对应孔1的组数sectionNum2 =[1,1];%对应孔2的组数
sectionNum1 =5.*ones(length(variant_n1),1);%对应孔1的组数
sectionNum2 =5.*ones(length(variant_n1),1);%对应孔2的组数
variant_n2 = [5:90];

iterateSize = length(variant_n1);
for i = 1:iterateSize    
    para(i).opt = opt;
    para(i).L1 = 3;%L1(m)
    para(i).L2 = 6;%L2（m）长度
    para(i).Dpipe = 0.098;%管道直径（m）
    para(i).vhpicStruct.l = 0.01;
    para(i).vhpicStruct.Dv = 0.372;%缓冲罐的直径（m）
    para(i).vhpicStruct.Lv = 1.1;%缓冲罐总长
    para(i).vhpicStruct.lc = 0.005;%内插管壁厚
    para(i).vhpicStruct.dp = 0.013;%开孔径
    para(i).vhpicStruct.Lin = 0.25;%内插管入口段长度
    para(i).vhpicStruct.lp1 = 0.18;%内插管入口段非孔管开孔长度
    para(i).vhpicStruct.lp2 = 0.18;%内插管出口段孔管开孔长度
    para(i).vhpicStruct.n1 = variant_n1(i);%入口段孔数
    para(i).vhpicStruct.n2 = variant_n2(i);%出口段孔数
    para(i).vhpicStruct.la1 = 0.02;%孔管入口段靠近入口长度
    para(i).vhpicStruct.la2 = 0.05;%孔管
    para(i).vhpicStruct.lb1 = 0.05;
    para(i).vhpicStruct.lb2 = 0.02;
    para(i).vhpicStruct.Din = 0.053;
    para(i).vhpicStruct.Lout = 0.25;
    para(i).vhpicStruct.bp1 = para(i).vhpicStruct.n1.*(para(i).vhpicStruct.dp)^2./(4.*para(i).vhpicStruct.Din.*para(i).vhpicStruct.lp1);%开孔率
    para(i).vhpicStruct.bp2 = para(i).vhpicStruct.n2.*(para(i).vhpicStruct.dp)^2./(4.*para(i).vhpicStruct.Din.*para(i).vhpicStruct.lp2);%开孔率
    l = para(i).vhpicStruct.lp1;
    para(i).vhpicStruct.xSection1 = [0,ones(1,sectionNum1(i)).*(l/(sectionNum1(i)))];
    l = para(i).vhpicStruct.lp2;
    para(i).vhpicStruct.xSection2 = [0,ones(1,sectionNum2(i)).*(l/(sectionNum2(i)))];
    para(i).sectionL1 = 0:0.25:para(i).L1;
    para(i).sectionL2 = 0:0.25:para(i).L2;
%%
%     para(i).opt = opt;
%     para(i).L1 = 3.85;%L1(m)
%     para(i).L2 = 10;%L2（m）长度
%     para(i).Dpipe = 0.106;%管道直径（m）
%     para(i).vhpicStruct.l = 0.01;
%     para(i).vhpicStruct.Dv = 0.45;%缓冲罐的直径（m）
%     para(i).vhpicStruct.Lv = 1.18;%缓冲罐总长
%     para(i).vhpicStruct.lc = 0.005;%内插管壁厚
%     para(i).vhpicStruct.dp = 0.02;%开孔径
%     para(i).vhpicStruct.lp1 = 0.04;%内插管入口段非孔管开孔长度
%     para(i).vhpicStruct.Lin = 0.2;%内插管入口段长度
%     para(i).vhpicStruct.n1 = 20;%入口段孔数
%         para(i).vhpicStruct.Din = 0.106;
%     para(i).vhpicStruct.lp2 = variant_lp2(i);%内插管出口段非孔管开孔长度
%     para(i).vhpicStruct.bp1 = para(i).vhpicStruct.n1.*(para(i).vhpicStruct.dp)^2./(4.*para(i).vhpicStruct.Din.*para(i).vhpicStruct.lp1);%开孔率
%     para(i).vhpicStruct.n2 = variant_n2(i);%出口段孔数
%     para(i).vhpicStruct.la1 = 0.01;%孔管入口段靠近入口长度,当单独考虑lp1的影响时，必须保证la1始终不变，或la2始终不变
%     para(i).vhpicStruct.la2 = para(i).vhpicStruct.Lin-para(i).vhpicStruct.la1-para(i).vhpicStruct.lp1;%假定保证入口段la1始终不变，0.15
% %     para(i).vhpicStruct.lb1 = 0.04;
% %     para(i).vhpicStruct.lb2 = 0.02;
%     para(i).vhpicStruct.Lout = 0.2;
% %     para(i).vhpicStruct.bp1 = para(i).vhpicStruct.n1.*(para(i).vhpicStruct.dp)^2./(4.*para(i).vhpicStruct.Din.*para(i).vhpicStruct.lp1);%开孔率
%     l = para(i).vhpicStruct.lp1;
%      para(i).vhpicStruct.xSection1 = [0,ones(1,sectionNum1(i)).*(l/(sectionNum1(i)))];
% %     l = para(i).vhpicStruct.lp2;
% %     para(i).vhpicStruct.xSection2 = [0,ones(1,sectionNum1(i)).*(l/(sectionNum1(i)))];
%     para(i).sectionL1 = 0:0.5:para(i).L1;
%     para(i).sectionL2 = 0:0.5:para(i).L2;
    
    holepipeLength1 = para(i).vhpicStruct.Lin - para(i).vhpicStruct.la1 - para(i).vhpicStruct.la2;
    hl1 = sum(para(i).vhpicStruct.xSection1);
    if(~cmpfloat(holepipeLength1,hl1))
        error('孔管参数设置错误：holepipeLength1=%.8f,hl1=%.8f;Lin:%g,la1:%g,la2:%g,sum(xSection1):%g,dp:%g'...
            ,holepipeLength1,hl1...
            ,para(i).vhpicStruct.Lin,para(i).vhpicStruct.la1,para(i).vhpicStruct.la2...
            ,sum(para(i).vhpicStruct.xSection1),para(i).vhpicStruct.dp);
    end
    name{i} = sprintf('n1:%g',variant_n1(i));
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
    pressure1Closed = [];
    pressure2Closed = [];

    [pressure1Closed,pressure2Closed] = ...
        vesselHaveInnerPerfBothClosedCompCalc(massFlowE,Fre,time,...
        para(i).L1,para(i).L2,para(i).Dpipe...
        ,para(i).vhpicStruct,...
        para(i).sectionL1,para(i).sectionL2,...
        'a',para(i).opt.acousticVelocity,'isDamping',para(i).opt.isDamping,'friction',para(i).opt.coeffFriction,...
        'meanFlowVelocity',para(i).opt.meanFlowVelocity,...
        'm',para(i).opt.mach,'notMach',para(i).opt.notMach,...
        'isOpening',isOpening);%,'coeffDamping',para(i).opt.coeffDamping,
    plus1Closed{i} = calcPuls(pressure1Closed,dcpss);
    plus2Closed{i} = calcPuls(pressure2Closed,dcpss);
    plusClosed{i} = [plus1Closed{i},plus2Closed{i}];
    multFreAmpValueClosed{i} = calcWaveFreAmplitude([pressure1Closed,pressure2Closed],Fs,multFre,'freErr',1);
    
    [pressure1Open,pressure2Open] = ...
        vesselHaveInnerPerfInCloOutOpenCompCalc(massFlowE,Fre,time,...
        para(i).L1,para(i).L2,para(i).Dpipe...
        ,para(i).vhpicStruct,...
        para(i).sectionL1,para(i).sectionL2,...
        'a',para(i).opt.acousticVelocity,'isDamping',para(i).opt.isDamping,'friction',para(i).opt.coeffFriction,...
        'meanFlowVelocity',para(i).opt.meanFlowVelocity,...
        'm',para(i).opt.mach,'notMach',para(i).opt.notMach...
        ,'isOpening',isOpening);%,'coeffDamping',opt.coeffDamping
    plus1Open{i} = calcPuls(pressure1Open,dcpss);
    plus2Open{i} = calcPuls(pressure2Open,dcpss);
    plusOpen{i} = [plus1Open{i},plus2Open{i}];
    multFreAmpValueOpend{i} = calcWaveFreAmplitude([pressure1Open,pressure2Open],Fs,multFre,'freErr',1);
    %计算缓冲罐内置孔管隔板内插管
%     if i == 1
%         [pressure1PI,pressure2PI] = vesselHavePerfInletNoneOutCalc(massFlowE,Fre,time,...
%             para(i).L1,para(i).L2,para(i).Dpipe...
%         ,para(i).vhpicStruct,...
%         para(i).sectionL1,para(i).sectionL2,...
%         'a',para(i).opt.acousticVelocity,'isDamping',para(i).opt.isDamping,'friction',para(i).opt.coeffFriction,...
%         'meanFlowVelocity',para(i).opt.meanFlowVelocity,...
%         'm',para(i).opt.mach,'notMach',para(i).opt.notMach,...
%         'isOpening',isOpening);
%         plus1PI = calcPuls(pressure1PI,dcpss);
%         plus2PI = calcPuls(pressure2PI,dcpss);
%         plusPI = [plus1PI,plus2PI];
%         multFreAmpValue_PI{i} = calcWaveFreAmplitude([pressure1PI,pressure2PI],Fs,[10,20,30],'freErr',1);
%     end
    %计算单一缓冲罐
    if i == 1
        [pressure1OV,pressure2OV] = oneVesselPulsationCalc(massFlowE,Fre,time,...
            para(i).L1,para(i).L2,...
            para(i).vhpicStruct.Lv,para(i).vhpicStruct.l,para(i).Dpipe,para(i).vhpicStruct.Dv,...
            para(i).sectionL1,para(i).sectionL2,...
            'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,...
            'meanFlowVelocity',opt.meanFlowVelocity,'isUseStaightPipe',1,...
            'm',para(i).opt.mach,'notMach',para(i).opt.notMach,...
            'isOpening',isOpening);
        plus1OV = calcPuls(pressure1OV,dcpss);
        plus2OV = calcPuls(pressure2OV,dcpss);
        plusOV = [plus1OV,plus2OV];
        multFreAmpValue_OV{i} = calcWaveFreAmplitude([pressure1OV,pressure2OV],Fs,multFre,'freErr',1);
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
        ,'m',para(i).opt.mach,'notMach',para(i).opt.notMach,...
        'isOpening',isOpening);
        plusStraight = calcPuls(temp,dcpss);
        maxPlus1Straight(i) = max(plusStraight(1:sepratorIndex(i)));
        maxPlus2Straight(i) = max(plusStraight(sepratorIndex(i):end));
        multFreAmpValue_straightPipe{i} = calcWaveFreAmplitude(temp,Fs,multFre,'freErr',1);
    end
    %计算脉动抑制率
    temp = plusStraight;
    temp2 = plusClosed{i};

    temp(temp<1e-4) = 1;
    temp2(temp<1e-4) = 1;%temp小于1e-4时，temp2也设置为1.
    reduceRate{i} = (temp - temp2)./temp;

    if isempty(plus1Closed{i})
        maxPlus1(i) = nan;
    else
        maxPlus1(i) = max(plus1Closed{i});
    end

    if isempty(plus2Closed{i})
        maxPlus2(i) = nan;
    else
        maxPlus2(i) = max(plus2Closed{i});
    end  

end

%% 绘制闭口迭代三维图
figure
hold on;
for i=1:iterateSize
	Z = plusClosed{i};
	if isXShowRealLength
        X = straightPipeSection;
    else
        X = 1:length(Z);
    end
    Y = ones(length(Z),1) .* variant_n1(i);
    plot3(X,Y,Z);
end
xlabel('距离');
ylabel('孔数');
zlabel('脉动压力');
title('闭口');
view(17,44);
set(gcf,'color','w');
%% 绘制开口迭代三维图
figure
hold on;
for i=1:iterateSize
	Z = plusOpen{i};
	if isXShowRealLength
        X = straightPipeSection;
    else
        X = 1:length(Z);
    end
    Y = ones(length(Z),1) .* variant_n1(i);
    plot3(X,Y,Z);
end
xlabel('距离');
ylabel('孔数');
zlabel('脉动压力');
title('开口');
view(17,44);
set(gcf,'color','w');
% marker_style = {'-o','-d','-<','-s','->','-<','-p','-*','-v','-^','-+','-x','-h'};
% color_style = [...
%     245,18,103;...
%     36,100,196;...
%     18,175,134;...
%     237,144,10;... 
%     131,54,229;...
    
%     255,99,56;...
%     ]./255;
% market_style_length = length(marker_style);
% color_style_length = size(color_style,1);
% marker_index = 0;
% color_index = 0;

% figure
% marker_index = 0;
% color_index = 0;
% h = [];
% textLegend = {};
% hold on;
% plotCount = 1;

% for i = 1:length(para)
    
%     YClose = plusClosed{i};
%     YClose = YClose./1000;
%     YClose=YClose';
%     YOpen = plusOpen{i};
%     YOpen = YOpen./1000;
%     YOpen=YOpen';
%     if isXShowRealLength
%         X = straightPipeSection;
%     else
%         X = 1:length(YClose);
%     end
%     if i==1 
%         h(plotCount) = plot(X,plusStraight./1000,'-r','LineWidth',1.5);
%         textLegend{plotCount} = '直管';
%         plotCount = plotCount + 1;
%     end
%     if i==1
%         %显示缓冲罐
%         h(plotCount) = plot(X,plusOV./1000,'-b','LineWidth',1.5);
%         textLegend{plotCount} = '单一缓冲罐';
%         plotCount = plotCount + 1;
%     end
%     marker_index = marker_index + 1;
%     if marker_index > market_style_length
%         marker_index = 1;
%     end
%     color_index = color_index + 1;
%     if color_index > color_style_length
%         color_index = 1;
%     end
%     h(plotCount) = plot(X,YClose,marker_style{marker_index},'color',color_style(color_index,:),'LineWidth',1.5);
%     textLegend{plotCount} = [name{i} '-close'];
%     plotCount = plotCount + 1;
%     marker_index = marker_index + 1;
%     if marker_index > market_style_length
%         marker_index = 1;
%     end
%     color_index = color_index + 1;
%     if color_index > color_style_length
%         color_index = 1;
%     end
    
%     h(plotCount) = plot(X,YOpen,marker_style{marker_index},'color',color_style(color_index,:),'LineWidth',1.5);
%     textLegend{plotCount} = [name{i} '-open'];
%     plotCount = plotCount + 1;
% end
% legend(h,textLegend,0);
% set(gcf,'color','w');


% figure
% marker_index = 0;
% color_index = 0;
% h = [];
% textLegend = {};
% hold on;
% plotCount = 1;

% for i = 1:length(para)
%     marker_index = marker_index + 1;
%     if marker_index > market_style_length
%         marker_index = 1;
%     end
%     color_index = color_index + 1;
%     if color_index > color_style_length
%         color_index = 1;
%     end
%     YClose = multFreAmpValueClosed{i}(1,:);
%     YClose = YClose./1000;
%     YClose=YClose';
%     YOpen = multFreAmpValueOpend{i}(1,:);
%     YOpen = YOpen./1000;
%     YOpen=YOpen';
%     if isXShowRealLength
%         X = straightPipeSection;
%     else
%         X = 1:length(YClose);
%     end
%     if i==1 
%         h(plotCount) = plot(X,multFreAmpValue_straightPipe{i}(1,:)./1000,'-r','LineWidth',1.5);
%         textLegend{plotCount} = '直管';
%         plotCount = plotCount + 1;
%     end
%     if i==1
%         %显示缓冲罐
%         h(plotCount) = plot(X,multFreAmpValue_OV{i}(1,:)./1000,'-b','LineWidth',1.5);
%         textLegend{plotCount} = '单一缓冲罐';
%         plotCount = plotCount + 1;
%     end
%     if i==0
%         %显示内置孔管隔板内插管
%         h(plotCount) = plot(X, multFreAmpValue_PI{i}(1,:)./1000,'-g','LineWidth',1.5);
%         textLegend{plotCount} = '内置孔管隔板内插管';
%         plotCount = plotCount + 1;
%     end
    
%     marker_index = marker_index + 1;
%     if marker_index > market_style_length
%         marker_index = 1;
%     end
%     color_index = color_index + 1;
%     if color_index > color_style_length
%         color_index = 1;
%     end
%     h(plotCount) = plot(X,YClose,marker_style{marker_index},'color',color_style(color_index,:),'LineWidth',1.5);
%     textLegend{plotCount} = [name{i} '-close'];
%     plotCount = plotCount + 1;
    
%     marker_index = marker_index + 1;
%     if marker_index > market_style_length
%         marker_index = 1;
%     end
%     color_index = color_index + 1;
%     if color_index > color_style_length
%         color_index = 1;
%     end 
%     h(plotCount) = plot(X,YOpen,marker_style{marker_index},'color',color_style(color_index,:),'LineWidth',1.5);
%     textLegend{plotCount} = [name{i} '-open'];
%     plotCount = plotCount + 1;
% end
% legend(h,textLegend,0);
% set(gcf,'color','w');