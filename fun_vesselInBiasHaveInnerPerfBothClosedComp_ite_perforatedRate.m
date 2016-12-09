function resCell = fun_vesselInBiasHaveInnerPerfBothClosedComp_ite_perforatedRate(varargin)
[isDamping,varargin]= takeVararginProperty('isDamping',varargin,1);

%迭代穿孔率的相关参数
currentPath = fileparts(mfilename('fullpath'));
isOpening = 0;%管道闭口
%rpm = 300;outDensity = 1.9167;multFre=[10,20,30];%环境25度绝热压缩到0.2MPaG的温度对应密度
rpm = 420;outDensity = 1.5608;multFre=[14,28,42];%环境25度绝热压缩到0.15MPaG的温度对应密度
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


opt.acousticVelocity = 345;%声速
opt.isDamping = isDamping;%是否计算阻尼
opt.coeffDamping = nan;%阻尼
opt.coeffFriction = 0.04;%管道摩察系数
opt.SreaightMeanFlowVelocity =20;%14.5;%管道平均流速
opt.SreaightCoeffFriction = 0.03;
opt.VesselMeanFlowVelocity =8;%14.5;%缓冲罐平均流速
opt.VesselCoeffFriction = 0.003;
opt.PerfClosedMeanFlowVelocity =9;%14.5;%堵死孔管平均流速
opt.PerfClosedCoeffFriction = 0.04;
opt.PerfOpenMeanFlowVelocity =15;%14.5;%开口孔管平均流速
opt.PerfOpenCoeffFriction = 0.035;

opt.isUseStaightPipe = 1;%计算容器传递矩阵的方法
opt.mach = opt.meanFlowVelocity / opt.acousticVelocity;
opt.notMach = 1;

L1 = 3.5;%L1(m)
L2 = 6;%L2（m）长度
Dpipe = 0.098;%管道直径（m）
l = 0.01;
Dv = 0.372;%缓冲罐的直径（m）
Lv = 1.1;%缓冲罐总长 
Lv1 =Lv./2;%缓冲罐腔1总长
Lv2 = Lv-Lv1;%缓冲罐腔2总长
lc = 0.005;%内插管壁厚
dp1 = variant_dp1;%开孔径
dp2 = variant_dp2;%开孔径
%     Lin = 0.25;%内插管入口段长度
lp1 = variant_lp1(i);%内插管入口段非孔管开孔长度
lp2 = variant_lp2(i);%内插管出口段孔管开孔长度
n1 = variant_n1;%入口段孔数
n2 = variant_n2;%出口段孔数
la1 = 0.03;%孔管入口段靠近入口长度
lb2 = 0.06;
la2 = 0.06;
lb1 = 0.03;
Din = variant_Din;
%     Lout = 0.25;
Lin = la1+lp1+la2;
Lout = lb1+lp2+lb2;
bp1 = variant_n1.*(variant_dp1)^2./(4.*variant_Din.*variant_lp1(i));%开孔率
bp2 = variant_n2.*(variant_dp2)^2./(4.*variant_Din.*variant_lp2(i));%开孔率
nc1 = 8;%假设一圈有8个孔
nc2 = 8;%假设一圈有8个孔
Cloum1 = variant_n1./nc1;%计算一端固定开孔长度的孔管上能开多少圈孔
Cloum2 = variant_n2./nc2;
s1 = ((variant_lp1(i)./Cloum1)-variant_dp1)./2;%相邻两开孔之间间隔，默认等间隔
s2 = ((variant_lp2(i)./Cloum2)-variant_dp2)./2;
sc1 = (pi.*variant_Din - nc1.*dp1)./nc1;%一周开孔，相邻孔间距
sc2 = (pi.*variant_Din - nc2.*dp2)./nc2;
l = lp1;
xSection1 = [0,ones(1,sectionNum1).*(l/(sectionNum1))];
l = lp2;
xSection2 = [0,ones(1,sectionNum2).*(l/(sectionNum2))];
sectionL1 = 0:0.25:L1;
sectionL2 = 0:0.25:L2;
lv1 = Lv./2-0.232;%232
lv2 = 0;%出口不偏置
Dbias = 0;%无内插管
end


function data = calcPerforatingRatios(n,dp,Din,lp)
% 计算开孔率
% n 孔数
% dp 孔管每一个孔孔径
% Din 孔管管径
% lp 孔管开孔长度
data = (n.*(dp).^2)./(4.*Din.*lp);%开孔率
end


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


function calcDatas = funIteratorVesselPipeLinePlusCalc(time,massFlow,Fs,plusBaseFrequency,acousticVelocity...
    ,meanFlowVelocity ...
    ,L1,L2,l,Dpipe,Dv,Lv,Lv1,lc ...
    ,dp1,dp2,n1,n2,la1,la2,lb1,lb2,Din ...
    ,IteratorValueName,varargin)
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
Lv2 = Lv - Lv1;%获取缓冲罐另一边长
[multFre,varargin]= takeVararginProperty('multfre',varargin,[10,20,30]);
[isOpening,varargin]= takeVararginProperty('isOpening',varargin,1);
[useCalcTopFreIndex,varargin]= takeVararginProperty('isOpening',varargin,nan);

[FreRaw,AmpRaw,~,massFlowE] = frequencySpectrum(detrend(massFlow,'constant'),Fs);
Fre = FreRaw;
% 提取主要频率
[~,locs] = findpeaks(AmpRaw,'SORTSTR','descend');
Fre = FreRaw(locs);
massFlowE = massFlowE(locs);

if ~isnan(useCalcTopFreIndex)
    Fre = Fre(useCalcTopFreIndex);
    massFlowE = massFlowE(useCalcTopFreIndex);
end
%是否带阻尼
[isDamping,varargin]= takeVararginProperty('isDamping',varargin,1);
%阻尼系数
[coeffDamping,varargin]= takeVararginProperty('coeffDamping',varargin,1);
%管道摩察系数
[coeffFriction,varargin]= takeVararginProperty('coeffFriction',varargin,1);
%计算脉动峰峰值的设置
[dcpss,varargin]= takeVararginProperty('dcpss',varargin,'default');
%计算的精度默认为1，就是每隔1m取一个点
[accuracy,varargin]= takeVararginProperty('accuracy',varargin,1);

[sectionNum1,varargin]= takeVararginProperty('sectionNum1',varargin,1);
[sectionNum2,varargin]= takeVararginProperty('sectionNum2',varargin,1);
%L1 L2的间隔
[sectionL1Interval,varargin]= takeVararginProperty('sectionL1Interval',varargin,0.5);
[sectionL2Interval,varargin]= takeVararginProperty('sectionL2Interval',varargin,0.5);
[lv1,varargin]= takeVararginProperty('lv1',varargin,0);%偏置长度
[lv2,varargin]= takeVararginProperty('lv2',varargin,0);
[Dbias,varargin]= takeVararginProperty('Dbias',varargin,0);%内插管长
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


opt.frequency = plusBaseFrequency;%脉动频率
opt.acousticVelocity = acousticVelocity;%声速
opt.isDamping = isDamping;%是否计算阻尼
opt.coeffDamping = coeffDamping;%阻尼
opt.coeffFriction = coeffFriction;%管道摩察系数
opt.meanFlowVelocity =meanFlowVelocity;%14.5;%管道平均流速
opt.isUseStaightPipe = 1;%计算容器传递矩阵的方法
opt.mach = opt.meanFlowVelocity / opt.acousticVelocity;
opt.notMach = 0;

para(1:iteCount).opt = opt;
para(1:iteCount).L1 = L1;%L1(m)
para(1:iteCount).L2 = L2;%L2（m）长度
para(1:iteCount).Dpipe = Dpipe;%管道直径（m）
para(1:iteCount).l = l;
para(1:iteCount).Dv = Dv;%缓冲罐的直径（m）
para(1:iteCount).Lv = Lv;%缓冲罐总长

para(1:iteCount).Lv1 = Lv1;%缓冲罐腔1总长
para(1:iteCount).Lv2 = Lv2;%缓冲罐腔2总长
para(1:iteCount).lc = lc;%内插管壁厚
para(1:iteCount).dp1 = dp1;%开孔径
para(1:iteCount).dp2 = dp2;%开孔径
para(1:iteCount).lp1 = lp1;%内插管入口段非孔管开孔长度
para(1:iteCount).lp2 = lp2;%内插管入口段非孔管开孔长度
para(1:iteCount).n1 = n1;%入口段孔数
para(1:iteCount).n2 = n2;%出口段孔数
para(1:iteCount).la1 = la1;%孔管入口段靠近入口长度
para(1:iteCount).la2 = la2;%孔管
para(1:iteCount).lb1 = lb1;%孔管
para(1:iteCount).lb2 = lb2;%孔管
para(1:iteCount).Din = Din;%孔管
para(1:iteCount).Lin = para(1:iteCount).la1 + para(1:iteCount).lp1 +para(1:iteCount).la2;%内插管入口段长度
para(1:iteCount).Lout = para(1:iteCount).lb1 + para(1:iteCount).lp2 +para(1:iteCount).lb2;%内插管入口段长度
para(1:iteCount).bp1 = calcPerforatingRatios(n1,dp1,Din,lp1);%开孔率
para(1:iteCount).bp2 = calcPerforatingRatios(n2,dp2,Din,lp2);%开孔率



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

dataCount = 2; index = 2;

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

index = index + 1; preMaxPlusIndex = index;
calcDatas{1,preMaxPlusIndex} = '罐前压力脉动最大值';

index = index + 1; backMaxPlusIndex = index;
calcDatas{1,backMaxPlusIndex} = '罐后压力脉动最大值';

index = index + 1; indexDp = index;
calcDatas{1,indexDp} = 'Dp'; 

index = index + 1; indexDin = index;
calcDatas{1,indexDin} = 'Din';

index = index + 1; indexlp = index;
calcDatas{1,indexlp} = 'lp';

index = index + 1; indexn1 = index;
calcDatas{1,indexn1} = 'n1';

index = index + 1; indexlb1 = index;
calcDatas{1,indexlb1} = 'lb1';

index = index + 1; indexla1 = index;
calcDatas{1,indexla1} = 'la1';

for i = 1:iteCount
    
     %计算单一缓冲罐
     X = [para(i).sectionL1,para(i).L1 + 2*para(i).l+para(i).Lv + para(i).sectionL2];
     if 1 == i
        [pressure1OV,pressure2OV] = oneVesselPulsationCalc(massFlowE,Fre,time,...
            para(i).L1,para(i).L2,...
            para(i).Lv,para(i).l,para(i).Dpipe,para(i).Dv,...
            para(i).sectionL1,para(i).sectionL2,...
            'a',opt.acousticVelocity,'isDamping',opt.isDamping,'friction',opt.coeffFriction,...
            'meanFlowVelocity',opt.meanFlowVelocity,'isUseStaightPipe',1,...
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
        temp = eval(sprintf('%s(i)',IteratorValueName));
        calcDatas{dataCount,rawIndex} = [pressure1OV,pressure2OV];
        calcDatas{dataCount,xIndex} = X;
        calcDatas{dataCount,plusIndex} = [plus1OV,plus2OV];
        calcDatas{dataCount,fre1Index} = multFreAmpValue_OV(1,:);
        calcDatas{dataCount,fre2Index} = multFreAmpValue_OV(2,:);
        calcDatas{dataCount,fre3Index} = multFreAmpValue_OV(3,:);
        calcDatas{dataCount,preMaxPlusIndex} = maxPlus1;
        calcDatas{dataCount,backMaxPlusIndex} = maxPlus2;
        calcDatas{dataCount,backMaxPlusIndex} = maxPlus2;
        dataCount = dataCount + 1;
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
    vhpicStruct.nc1 = para(i).nc1;
    vhpicStruct.nc2 = para(i).nc2;
    vhpicStruct.xSection1 = para(i).xSection1;
    vhpicStruct.xSection2 = para(i).xSection2;
    vhpicStruct.lv1 = lv1;
    vhpicStruct.lv2 = lv2;
    vhpicStruct.Dbias = Dbias;
    [pressure1ClosedIB,pressure2ClosedIB] = ...
         vesselInBiasHaveInnerPerfBothClosedCompCalc(massFlowE,Fre,time,...
        para(i).L1,para(i).L2,para(i).Dpipe...
        ,vhpicStruct...
        ,para(i).sectionL1,para(i).sectionL2,...
        'a',para(i).opt.acousticVelocity,'isDamping',para(i).opt.isDamping,'friction',PerfClosedCoeffFriction,...
        'meanFlowVelocity',PerfClosedMeanFlowVelocity,...
        'm',para(i).opt.mach,'notMach',para(i).opt.notMach,...
        'isOpening',isOpening);%,'coeffDamping',para(i).opt.coeffDamping,
end
end
