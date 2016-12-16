%% 题目是Lin，其实是迭代Lv1-Lin部分
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
%rpm = 300;outDensity = 1.9167;multFre=[10,20,30];%环境25度绝热压缩到0.2MPaG的温度对应密度
rpm = 420;outDensity = 1.5608;multFre=[14,28,42];%环境25度绝热压缩到0.15MPaG的温度对应密度
Fs = 4096;
[massFlowRaw,time,~,opt.meanFlowVelocity] = massFlowMaker(0.25,0.098,rpm...
	,0.14,1.075,outDensity,'rcv',0.15,'k',1.4,'pr',0.15,'fs',Fs,'oneSecond',6);

%massFlow = load(fullfile(currentPath,'mass_flow_0.1478_NorthZone.txt'));

[FreRaw,AmpRaw,PhRaw,massFlowERaw] = frequencySpectrum(detrend(massFlowRaw,'constant'),Fs);
FreRaw = [7,14,21,28,14*3];
massFlowERaw = [0.02,0.2,0.03,0.003,0.007];

% 提取主要频率
massFlowE = massFlowERaw;
Fre = FreRaw;

% [pks,locs] = findpeaks(AmpRaw,'SORTSTR','descend');
% Fre = FreRaw(locs);
% massFlowE = massFlowERaw(locs);
% temp = [1:20];%(Fre<29) ;%| (Fre>30 & Fre < 100);
% 
% massFlowE4Vessel = massFlowE;
% massFlowE = massFlowE(temp);
% Fre4Vessel = Fre;
% Fre = Fre(temp);

isDamping = 1;
%绘图参数
isXShowRealLength = 1;
isShowStraightPipe=1;%是否显示直管
isShowOnlyVessel=1;%是否显示无内件缓冲罐

opt.acousticVelocity = 345;%声速
opt.isDamping = isDamping;%是否计算阻尼
opt.coeffDamping = nan;%阻尼
opt.coeffFriction = 0.04;%管道摩察系数
SreaightMeanFlowVelocity =20;%14.5;%管道平均流速
SreaightCoeffFriction = 0.03;
VesselMeanFlowVelocity =8;%14.5;%缓冲罐平均流速
VesselCoeffFriction = 0.003;
PerfClosedMeanFlowVelocity =9;%14.5;%堵死孔管平均流速
PerfClosedCoeffFriction = 0.04;
PerfOpenMeanFlowVelocity =15;%14.5;%开口孔管平均流速
PerfOpenCoeffFriction = 0.035;
% opt.meanFlowVelocity =14.5;%14.5;%管道平均流速
opt.isUseStaightPipe = 1;%计算容器传递矩阵的方法
opt.mach = opt.meanFlowVelocity / opt.acousticVelocity;
opt.notMach = 1;

variant_n1 = [32];              %variant_n = [6,6];sectionNum1 =[1,6];%对应孔1的组数sectionNum2 =[1,1];%对应孔2的组数
sectionNum1 =[1];%对应孔1的组数
sectionNum2 =[1];%对应孔2的组数
variant_n2 = [32];
variant_Lv1 = 0.26:0.1:0.84;
calcDatas = {};


for i = 1:length(variant_Lv1)     
    para(i).opt = opt;
    para(i).L1 = 3.5;%L1(m)
    para(i).L2 = 6;%L2（m）长度
    para(i).Dpipe = 0.098;%管道直径（m）
    para(i).vhpicStruct.l = 0.01;
    para(i).vhpicStruct.Dv = 0.372;%缓冲罐的直径（m）
    para(i).vhpicStruct.Lv = 1.1;%缓冲罐总长 
    para(i).vhpicStruct.Lv1 =variant_Lv1(i);%缓冲罐腔1总长
    para(i).vhpicStruct.Lv2 = para(i).vhpicStruct.Lv-para(i).vhpicStruct.Lv1;%缓冲罐腔2总长
    para(i).vhpicStruct.lc = 0.005;%内插管壁厚
    para(i).vhpicStruct.dp1 = 0.013;%开孔径
    para(i).vhpicStruct.dp2 = 0.013;%开孔径
    para(i).vhpicStruct.Lin = 0.25;%内插管入口段长度
    para(i).vhpicStruct.lp1 = 0.16;%内插管入口段非孔管开孔长度
    para(i).vhpicStruct.lp2 = 0.16;%内插管出口段孔管开孔长度
    para(i).vhpicStruct.n1 = variant_n1;%入口段孔数
    para(