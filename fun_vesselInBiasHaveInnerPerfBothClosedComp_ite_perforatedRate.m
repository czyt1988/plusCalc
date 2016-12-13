function [resN1,resLp1,resDin,resDp1,resOneVessel] = fun_vesselInBiasHaveInnerPerfBothClosedComp_ite_perforatedRate()

%迭代穿孔率的相关参数
currentPath = fileparts(mfilename('fullpath'));
isOpening = 0;%管道闭口
%rpm = 300;outDensity = 1.9167;multFre=[10,20,30];%环境25度绝热压缩到0.2MPaG的温度对应密度
rpm = 420;outDensity = 1.5608;%环境25度绝热压缩到0.15MPaG的温度对应密度
plusBaseFrequency = 2*(rpm/60);
multfre = [1,2,3] .* plusBaseFrequency;
Fs = 4096;
[massFlowRaw,time,~,opt.meanFlowVelocity] = massFlowMaker(0.25,0.098,rpm...
    ,0.14,1.075,outDensity,'rcv',0.15,'k',1.4,'pr',0.15,'fs',Fs,'oneSecond',6);
[FreRaw,AmpRaw,PhRaw,massFlowERaw] = frequencySpectrum(detrend(massFlowRaw,'constant'),Fs);
FreRaw = [7,14,21,28,14*3];
massFlowERaw = [0.02,0.2,0.03,0.003,0.007];
% 提取主要频率
massFlowE = massFlowERaw;
Fre = FreRaw;
%%参数

%计算脉动峰峰值的设置
dcpss = getDefaultCalcPulsSetStruct();
dcpss.calcSection = [0.3,0.7];
dcpss.fs = Fs;
dcpss.isHp = 0;
dcpss.f_pass = 7;%通过频率5Hz
dcpss.f_stop = 5;%截止频率3Hz
dcpss.rp = 0.1;%边带区衰减DB数设置
dcpss.rs = 30;%截止区衰减DB数设置




%

opt.frequency = plusBaseFrequency;%脉动频率
opt.acousticVelocity = 345;%声速
opt.isDamping = 1;%是否计算阻尼
opt.coeffDamping = nan;%阻尼
opt.coeffFriction = 0.04;%管道摩察系数
opt.meanFlowVelocity =9;%14.5;%管道平均流速
opt.isUseStaightPipe = 1;%计算容器传递矩阵的方法
opt.mach = opt.meanFlowVelocity / opt.acousticVelocity;
opt.notMach = 1;

opt.vesselCoeffFriction = 0.003;%缓冲罐管道摩察系数
opt.vesselMeanFlowVelocity =8;%14.5;%缓冲罐管道平均流速

%迭代
L1 = 3.5;%L1(m)
L2 = 6;%L2（m）长度
Dpipe = 0.098;%管道直径（m）
l = 0.01;
Dv = 0.372;%缓冲罐的直径（m）
Lv = 1.1;%缓冲罐总长 
Lv1 =Lv./2;%缓冲罐腔1总长
Lv2 = Lv - Lv1;
lc = 0.005;%内插管壁厚
dp1 = 0.013;%开孔径
dp2 = 0.013;%开孔径
%     Lin = 0.25;%内插管入口段长度
lp1 = 0.16;%内插管入口段非孔管开孔长度
lp2 = 0.16;%内插管出口段孔管开孔长度
n1 = 24;%入口段孔数
n2 = 24;%出口段孔数
la1 = 0.03;%孔管入口段靠近入口长度
la2 = 0.06;
lb2 = 0.03;
lb1 = 0.06;
Din = 0.049;
lv1 = Lv./2-0.232;%232
sectionNum1 = [1];%对应孔1的组数
sectionNum2 = [1];%对应孔2的组数
Dbias = 0;%无内插管
accuracy = 0.25;%计算精度
%迭代n1
n1Ite = [8,16,24,32,40,48,56,64,72,80,88];
temp = funIteratorVesselPipeLinePlusCalc(time,massFlowE,Fre,Fs,opt...
    ,L1,L2,l,Dpipe,Dv,Lv1,Lv2,lc,lp1,lp2,lv1...
    ,dp1,dp2,n1Ite,n2,la1,la2,lb1,lb2,Din...
    ,'n1'...
    ,'isOpening',0 ...
    ,'multfre',multfre...
    ,'dcpss',dcpss...
    ,'accuracy',accuracy...
    ,'sectionNum1',sectionNum1...
    ,'sectionNum1',sectionNum2...
    ,'Dbias',Dbias...
    ,'isCalcPureVessel',1 ...
);
resOneVessel = temp(1:2,:);
resN1 = {};
resN1 = cellPush2Bottom(resN1,temp(1,:),temp(3:end,:));

%迭代穿孔长度
lp1Ite = [0.05;0.06;0.07;0.08;0.10;0.12;0.16;0.24;0.32;0.40;0.48];
resLp1 = funIteratorVesselPipeLinePlusCalc(time,massFlowE,Fre,Fs,opt...
    ,L1,L2,l,Dpipe,Dv,Lv1,Lv2,lc,lp1Ite,lp2,lv1...
    ,dp1,dp2,n1,n2,la1,la2,lb1,lb2,Din...
    ,'lp1'...
    ,'isOpening',0 ...
    ,'multfre',multfre...
    ,'dcpss',dcpss...
    ,'accuracy',accuracy...
    ,'sectionNum1',sectionNum1...
    ,'sectionNum1',sectionNum2...
    ,'Dbias',Dbias...
    ,'isCalcPureVessel',0 ...
);

%迭代孔管管径
DinIte = [0.15*0.098,0.165*0.098,0.185*0.098,0.2*0.098,0.23*0.098,0.26*0.098,0.3*0.098,0.4*0.098,0.5*0.098,0.6*0.098,0.75*0.098,0.098,0.098*2];
resDin = funIteratorVesselPipeLinePlusCalc(time,massFlowE,Fre,Fs,opt...
    ,L1,L2,l,Dpipe,Dv,Lv1,Lv2,lc,lp1,lp2,lv1...
    ,dp1,dp2,n1,n2,la1,la2,lb1,lb2,DinIte...
    ,'Din'...
    ,'isOpening',0 ...
    ,'multfre',multfre...
    ,'dcpss',dcpss...
    ,'accuracy',accuracy...
    ,'sectionNum1',sectionNum1...
    ,'sectionNum1',sectionNum2...
    ,'Dbias',Dbias...
    ,'isCalcPureVessel',0 ...
);
%迭代孔径
Dp1Ite = [0.007,0.01,0.013,0.016,0.0175,0.019,0.0205,0.022,0.023,0.024,0.025];
resDp1 = funIteratorVesselPipeLinePlusCalc(time,massFlowE,Fre,Fs,opt...
    ,L1,L2,l,Dpipe,Dv,Lv1,Lv2,lc,lp1,lp2,lv1...
    ,Dp1Ite,dp2,n1,n2,la1,la2,lb1,lb2,Din...
    ,'dp1'...
    ,'isOpening',0 ...
    ,'multfre',multfre...
    ,'dcpss',dcpss...
    ,'accuracy',accuracy...
    ,'sectionNum1',sectionNum1...
    ,'sectionNum1',sectionNum2...
    ,'Dbias',Dbias...
    ,'isSavePureVessel',0 ...
);

end





function calcDatas = funIteratorVesselPipeLinePlusCalc(time,massFlowE,Fre,Fs,opt...
    ,L1,L2,l,Dpipe,Dv,Lv1,Lv2,lc,lp1,lp2,lv1 ...
    ,dp1,dp2,n1,n2,la1,la2,lb1,lb2,Din ...
    ,IteratorValueName,varargin)
%缓冲罐中间插入孔管,两端堵死，开孔个数不足以等效为亥姆霍兹共鸣器(缓冲罐入口偏置)
%                 L1
%                     |
%                     |
%           l         |          Lv              l    L2  
%              _______|_________________________        
%             |    dp1(n1)   |    dp2(n2)       |
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
% 缓冲罐内置孔管结构与单一顺接缓冲罐对比
% massFlow 为计算的质量流量
% meanFlowVelocity 入口管道的流速
% Fs 质量流量的采样率
% plusBaseFrequency 脉动基本频率
% acousticVelocity 声速
% L1 入口管道长度
% L2 出口管道长度
% l  缓冲罐两侧的延生长度
% Dpipe  管径
% Dv  缓冲罐直径
% Lv 缓冲罐长度
% Lv1 缓冲罐腔1总长
% lc 内插管壁厚
% dp1 开孔径1
% dp2 开孔径2
% n1 入口段孔数
% n2 出口段孔数
% la1 孔管入口段靠近入口长度
% la2
% lb1
% lb2
% Din
% Lin
% Lout
% IteratorValueName 参数长度不为1的变量名，就是要迭代的变量名
% varargin取值：
% 'multfre':需要计算的倍频(vector)([10,20,30])
% 'isOpening'：(bool)(1) 是否开口
% 'useCalcTopFreIndex':参与计算的最高频率数，如[1:20]代表前20个最大的频率进行计算，nan为全部参与
% 'isDamping' 是否有阻尼
% 'coeffDamping' 阻尼系数
% 'coeffFriction' 管道摩察系数
% 'dcpss' 计算脉动峰峰值的设置,如果为'default',会设置一个默认值,否则应该如下设置：
%       dcpss = getDefaultCalcPulsSetStruct();
%       dcpss.calcSection = [0.3,0.7];
%       dcpss.fs = Fs;
%       dcpss.isHp = 0;
%       dcpss.f_pass = 7;%通过频率5Hz
%       dcpss.f_stop = 5;%截止频率3Hz
%       dcpss.rp = 0.1;%边带区衰减DB数设置
%       dcpss.rs = 30;%截止区衰减DB数设置
%  'accuracy' 计算的精度默认为1，就是每隔1m取一个点
% 返回的cell
% {[]   ,'x值','压力脉动','1倍频','2倍频','3倍频','罐前压力脉动最大值','罐后压力脉动最大值'
% '直管',
% 'xx'

iteCount = eval(sprintf('length(%s)',IteratorValueName));

% 
% xSection1，xSection2 孔管每圈孔的间距，从0开始算，x的长度为孔管孔的圈数+1，x的值是当前一圈孔和上一圈孔的距离，如果间距一样，那么x里的值都一样
Lv = Lv2 + Lv1;%获取缓冲罐另一边长
lv2 = 0;
[multFre,varargin]= takeVararginProperty('multfre',varargin,[10,20,30]);
[isOpening,varargin]= takeVararginProperty('isOpening',varargin,0);

%计算脉动峰峰值的设置
[dcpss,varargin]= takeVararginProperty('dcpss',varargin,'default');
%计算的精度默认为1，就是每隔1m取一个点%L1 L2的间隔
[accuracy,varargin]= takeVararginProperty('accuracy',varargin,0.5);

[sectionNum1,varargin]= takeVararginProperty('sectionNum1',varargin,1);
[sectionNum2,varargin]= takeVararginProperty('sectionNum2',varargin,1);


% [lv1,varargin]= takeVararginProperty('lv1',varargin,0);%偏置长度
% [lv2,varargin]= takeVararginProperty('lv2',varargin,0);
[Dbias,varargin]= takeVararginProperty('Dbias',varargin,0);%内插管长
[isSavePureVessel,varargin]= takeVararginProperty('isSavePureVessel',varargin,0);%是否计算单一缓冲罐
if ~isstruct(dcpss)
    dcpss = getDefaultCalcPulsSetStruct();
    dcpss.calcSection = [0.3,0.7];
    dcpss.fs = Fs;
    dcpss.isHp = 0;
    dcpss.f_pass = 7;%通过频率5Hz
    dcpss.f_stop = 5;%截止频率3Hz
    dcpss.rp = 0.1;%边带区衰减DB数设置
    dcpss.rs = 30;%截止区衰减DB数设置
end




para = [];
para = setDataInStruct(para,opt,'opt',iteCount);
para = setDataInStruct(para,L1,'L1',iteCount);%L1(m)
para = setDataInStruct(para,L2,'L2',iteCount);%L2（m）长度
para = setDataInStruct(para,Dpipe,'Dpipe',iteCount);%管道直径（m）
para = setDataInStruct(para,l,'l',iteCount);
para = setDataInStruct(para,Dv,'Dv',iteCount);%缓冲罐的直径（m）
para = setDataInStruct(para,Lv,'Lv',iteCount);%缓冲罐总长
para = setDataInStruct(para,Lv1,'Lv1',iteCount);%缓冲罐腔1总长
para = setDataInStruct(para,Lv2,'Lv2',iteCount);%缓冲罐腔2总长
para = setDataInStruct(para,lc,'lc',iteCount);%内插管壁厚
para = setDataInStruct(para,dp1,'dp1',iteCount);%开孔径
para = setDataInStruct(para,dp2,'dp2',iteCount);%开孔径
para = setDataInStruct(para,lp1,'lp1',iteCount);%内插管入口段非孔管开孔长度
para = setDataInStruct(para,lp2,'lp2',iteCount);%内插管入口段非孔管开孔长度
para = setDataInStruct(para,n1,'n1',iteCount);%入口段孔数
para = setDataInStruct(para,n2,'n2',iteCount);%出口段孔数
para = setDataInStruct(para,n1,'n1',iteCount);
para = setDataInStruct(para,n2,'n2',iteCount);
para = setDataInStruct(para,la1,'la1',iteCount);%孔管入口段靠近入口长度
para = setDataInStruct(para,la2,'la2',iteCount);%孔管
para = setDataInStruct(para,lb1,'lb1',iteCount);
para = setDataInStruct(para,lb2,'lb2',iteCount);
para = setDataInStruct(para,Din,'Din',iteCount);

for i = 1:iteCount
    para(i).Lin =para(i).la1 + para(i).lp1+para(i).la2;%内插管入口段长度
    para(i).Lout = para(i).lb1 + para(i).lp2+para(i).lb2;%内插管入口段长度
end
para = setDataInStruct(para,calcPerforatingRatios(n1,dp1,Din,lp1),'bp1',iteCount);%开孔率
para = setDataInStruct(para,calcPerforatingRatios(n2,dp2,Din,lp2),'bp2',iteCount);%开孔率


for i=1:iteCount
    para(i).sectionL1 = 0:accuracy:para(i).L1;
    para(i).sectionL2 = 0:accuracy:para(i).L2;
    para(i).xSection1 = [0,ones(1,sectionNum1).*(para(i).lp1/(sectionNum1))];
    para(i).xSection2 = [0,ones(1,sectionNum2).*(para(i).lp2/(sectionNum2))];

    holepipeLength1 = para(i).Lin - para(i).la1 - para(i).la2;
    hl1 = sum(para(i).xSection1);
    if(~cmpfloat(holepipeLength1,hl1))
        error('孔管参数设置错误：holepipeLength1=%.8f,hl1=%.8f;Lin:%g,la1:%g,la2:%g,sum(xSection1):%g,dp:%g'...
            ,holepipeLength1,hl1...
            ,para(i).Lin,para(i).la1,para(i).la2...
            ,sum(para(i).xSection1),para(i).dp);
    end

end

dataCount = 2; index = 1;
calcDatas{1,1} = '说明\变量名';
indexName = index;
index = index + 1; rawIndex = index;
calcDatas{1,rawIndex} = 'rawData';

index = index + 1; xIndex = index;
calcDatas{1,xIndex} = 'x值';

index = index + 1; plusIndex = index;
calcDatas{1,plusIndex} = '压力脉动';

index = index + 1; fre1Index = index;
calcDatas{1,fre1Index} = '1倍频';

index = index + 1; fre2Index = index;
calcDatas{1,fre2Index} = '2倍频';

index = index + 1; fre3Index = index;
calcDatas{1,fre3Index} = '3倍频';

index = index + 1; 	restrainRateIndex = index;
calcDatas{1,restrainRateIndex} = '脉动抑制率';

index = index + 1; preMaxPlusIndex = index;
calcDatas{1,preMaxPlusIndex} = '罐前压力脉动最大值';

index = index + 1; backMaxPlusIndex = index;
calcDatas{1,backMaxPlusIndex} = '罐后压力脉动最大值';

index = index + 1; indexbp1 = index;
calcDatas{1,indexbp1} = 'bp1 - 开孔率';

index = index + 1; indexbp2 = index;
calcDatas{1,indexbp2} = 'bp2 - 开孔率';

index = index + 1; indexdp1 = index;
calcDatas{1,indexdp1} = 'dp1'; 

index = index + 1; indexdp2 = index;
calcDatas{1,indexdp2} = 'dp2'; 

index = index + 1; indexDin = index;
calcDatas{1,indexDin} = 'Din';

index = index + 1; indexlp1 = index;
calcDatas{1,indexlp1} = 'lp1';

index = index + 1; indexlp2 = index;
calcDatas{1,indexlp2} = 'lp2';

index = index + 1; indexn1 = index;
calcDatas{1,indexn1} = 'n1';

index = index + 1; indexlb1 = index;
calcDatas{1,indexlb1} = 'lb1';

index = index + 1; indexla1 = index;
calcDatas{1,indexla1} = 'la1';

index = index + 1; indexInput = index;
calcDatas{1,indexInput} = 'Input';

for i = 1:iteCount
    
     %计算单一缓冲罐
     X = [para(i).sectionL1,para(i).L1 + 2*para(i).l+para(i).Lv + para(i).sectionL2];
     if 1 == i
%             [pressure1OV,pressure2OV] = oneVesselPulsationCalc(massFlowE,Fre,time,...
%                 para(i).L1,para(i).L2,...
%                 para(i).Lv,para(i).l,para(i).Dpipe,para(i).Dv,...
%                 para(i).sectionL1,para(i).sectionL2,...
%                 'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,...
%                 'meanFlowVelocity',opt.meanFlowVelocity,'isUseStaightPipe',1,...
%                 'm',opt.mach,'notMach',opt.notMach,...
%                 'isOpening',isOpening);
           [pressure1OV,pressure2OV] = vesselBiasStraightPulsationCalc(massFlowE,Fre,time,...
                para(i).L1,para(i).L2...
                ,para(i).Lv,para(i).l,para(i).Dpipe,para(i).Dv...
                ,lv1,Dbias...
                ,para(i).sectionL1,para(i).sectionL2,...
                'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.vesselCoeffFriction,...
                'meanFlowVelocity',opt.vesselMeanFlowVelocity,'isUseStaightPipe',1,...
                'm',opt.mach,'notMach',opt.notMach,...
                'isOpening',isOpening);
            plus1OV = calcPuls(pressure1OV,dcpss);
            plus2OV = calcPuls(pressure2OV,dcpss);
            plusOV = [plus1OV,plus2OV];
            if isempty(plus1OV)
                maxPlus1 = nan;
            else
                maxPlus1= max(plus1OV);
            end
            if isempty(plus2OV)
                maxPlus2 = nan;
            else
                maxPlus2 = max(plus2OV);
            end  

            multFreAmpValue_OV = calcWaveFreAmplitude([pressure1OV,pressure2OV],Fs,multFre,'freErr',1);
         if isSavePureVessel
            calcDatas{dataCount,indexName} = '无内件缓冲罐';
            calcDatas{dataCount,rawIndex} = [pressure1OV,pressure2OV];
            calcDatas{dataCount,xIndex} = X;
            calcDatas{dataCount,plusIndex} = plusOV;
            calcDatas{dataCount,fre1Index} = multFreAmpValue_OV(1,:);
            calcDatas{dataCount,fre2Index} = multFreAmpValue_OV(2,:);
            calcDatas{dataCount,fre3Index} = multFreAmpValue_OV(3,:);
            calcDatas{dataCount,preMaxPlusIndex} = maxPlus1;
            calcDatas{dataCount,backMaxPlusIndex} = maxPlus2;

            dataCount = dataCount + 1;
         end
    end
    %计算孔管的数据
    vhpicStruct.l  = para(i).l;
    vhpicStruct.Dv = para(i).Dv;
    vhpicStruct.Lv = para(i).Lv;
    vhpicStruct.Lv1 = para(i).Lv1;
    vhpicStruct.Lv2 = para(i).Lv2;
    vhpicStruct.lc   = para(i).lc  ;
    vhpicStruct.dp1  = para(i).dp1 ;
    vhpicStruct.dp2  = para(i).dp2 ;
    vhpicStruct.Lin  = para(i).Lin ;
    vhpicStruct.lp1  = para(i).lp1 ;
    vhpicStruct.lp2  = para(i).lp2 ;
    vhpicStruct.n1   = para(i).n1  ;
    vhpicStruct.n2   = para(i).n2  ;
    vhpicStruct.la1  = para(i).la1 ;
    vhpicStruct.la2  = para(i).la2 ;
    vhpicStruct.lb1  = para(i).lb1 ;
    vhpicStruct.lb2  = para(i).lb2 ;
    vhpicStruct.Din  = para(i).Din ;
    vhpicStruct.Lout = para(i).Lout;
    vhpicStruct.bp1 = para(i).bp1;
    vhpicStruct.bp2 = para(i).bp2;
    vhpicStruct.xSection1 = para(i).xSection1;
    vhpicStruct.xSection2 = para(i).xSection2;
    vhpicStruct.lv1 = lv1;
    vhpicStruct.lv2 = lv2;
    vhpicStruct.Dbias = Dbias;
    [pressure1,pressure2] = ...
         vesselInBiasHaveInnerPerfBothClosedCompCalc(massFlowE,Fre,time,...
        para(i).L1,para(i).L2,para(i).Dpipe...
        ,vhpicStruct...
        ,para(i).sectionL1,para(i).sectionL2,...
        'a',para(i).opt.acousticVelocity,'isDamping',para(i).opt.isDamping,'friction',opt.coeffFriction,...
        'meanFlowVelocity',opt.meanFlowVelocity,...
        'm',opt.mach,'notMach',opt.notMach,...
        'isOpening',isOpening);%,'coeffDamping',para(i).opt.coeffDamping,
    pressure = [pressure1,pressure2];
    plus1 = calcPuls(pressure1,dcpss);
    plus2 = calcPuls(pressure2,dcpss);
    plus = [plus1,plus2];
    if isempty(plus1)
        maxPlus1 = nan;
    else
        maxPlus1= max(plus1);
    end
    if isempty(plus2)
        maxPlus2 = nan;
    else
        maxPlus2 = max(plus2);
    end  

    multFreAmpValue = calcWaveFreAmplitude(pressure,Fs,multFre,'freErr',1);
    temp = eval(sprintf('%s(i)',IteratorValueName));
    calcDatas{dataCount,indexName} = temp;
    calcDatas{dataCount,rawIndex} = pressure;
    calcDatas{dataCount,xIndex} = X;
    calcDatas{dataCount,plusIndex} = plus;
    calcDatas{dataCount,fre1Index} = multFreAmpValue(1,:);
    calcDatas{dataCount,fre2Index} = multFreAmpValue(2,:);
    calcDatas{dataCount,fre3Index} = multFreAmpValue(3,:);
    calcDatas{dataCount,restrainRateIndex} = (plusOV - plus) ./ plusOV;
    calcDatas{dataCount,preMaxPlusIndex} = maxPlus1;
    calcDatas{dataCount,backMaxPlusIndex} = maxPlus2;
    calcDatas{dataCount,indexbp1} = vhpicStruct.bp1;
    calcDatas{dataCount,indexbp2} = vhpicStruct.bp2;
    calcDatas{dataCount,indexdp1} = vhpicStruct.dp1;
    calcDatas{dataCount,indexdp2} = vhpicStruct.dp2;
    calcDatas{dataCount,indexDin} = vhpicStruct.Din;
    calcDatas{dataCount,indexlp1} = vhpicStruct.lp1;
    calcDatas{dataCount,indexlp2} = vhpicStruct.lp2;
    calcDatas{dataCount,indexn1} = vhpicStruct.n1;
    calcDatas{dataCount,indexlb1} = vhpicStruct.lb1;
    calcDatas{dataCount,indexla1} = vhpicStruct.la1;
    calcDatas{dataCount,indexInput} = vhpicStruct.la1;
    dataCount = dataCount + 1;
end

end

function param = setDataInStruct(param,data,dataName,iteNum)
    if length(data) > 1
        for i = 1:iteNum
            eval(sprintf('param(i).%s = data(i);',dataName));
        end
    else
        eval(sprintf('[param(1:iteNum).%s] = deal(data);',dataName));
    end
end

