function resCell = fun_vesselInBiasHaveInnerPerfBothClosedComp_ite_perforatedRate

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
    ,L1,L2,l,Dpipe,Dv,Lv,Lv1,lc ...
    ,dp1,dp2,n1,n2,la1,la2,lb1,lb2,Din,Lin,Lout
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
%% 缓冲罐内置孔管结构与单一顺接缓冲罐对比
% massFlow 为计算的质量流量
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
%{[]   ,'x值','压力脉动','1倍频','2倍频','3倍频','罐前压力脉动最大值','罐后压力脉动最大值'
% '直管',
% 'xx'

iteCount = eval(sprintf('length(%s)',IteratorValueName));

% 
% xSection1，xSection2 孔管每圈孔的间距，从0开始算，x的长度为孔管孔的圈数+1，x的值是当前一圈孔和上一圈孔的距离，如果间距一样，那么x里的值都一样
Lv2 = Lv - Lv1;%获取缓冲罐另一边长
[multFre,varargin]= takeVararginProperty('multfre',varargin,[10,20,30]);
[isOpening,varargin]= takeVararginProperty('isOpening',varargin,1);
[useCalcTopFreIndex,varargin]= takeVararginProperty('isOpening',varargin,nan);

[FreRaw,AmpRaw,PhRaw,massFlowE] = frequencySpectrum(detrend(massFlow,'constant'),Fs);
Fre = FreRaw;
% 提取主要频率
[pks,locs] = findpeaks(AmpRaw,'SORTSTR','descend');
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
opt.meanFlowVelocity =14.5;%14.5;%管道平均流速
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
for i=1:iteCount
    para(i).sectionL1 = 0:accuracy:para(i).L1;
    para(i).sectionL2 = 0:accuracy:para(i).L2;
end

dataCount = 2;
calcDatas{1,2} = 'x值';
calcDatas{1,3} = '压力脉动';
calcDatas{1,4} = '1倍频';
calcDatas{1,5} = '2倍频';
calcDatas{1,6} = '3倍频';
calcDatas{1,7} = '罐前压力脉动最大值';
calcDatas{1,8} = '罐后压力脉动最大值';
for i = 1:iteCount
    if i==1
        %计算直管
        %直管总长
        straightPipeLength = para(i).L1 + 2*para(i).l+para(i).Lv + para(i).L2;
        straightPipeSection = [para(i).sectionL1,...
                                para(i).L1 + 2*para(i).l+para(i).Lv + para(i).sectionL2];
    
        temp = find(para(i).L1>straightPipeSection);%找到缓冲罐所在的索引
        sepratorIndex = temp(end);
        temp = straightPipePulsationCalc(massFlowE,Fre,time,straightPipeLength,straightPipeSection...
        ,'d',para(i).Dpipe,'a',opt.acousticVelocity,'isDamping',opt.isDamping...
        ,'friction',opt.coeffFriction,'meanFlowVelocity',opt.meanFlowVelocity...
        ,'m',opt.mach,'notMach',opt.notMach,...
        'isOpening',isOpening);
        plusStraight = calcPuls(temp,dcpss);
        maxPlus1Straight(i) = max(plusStraight(1:sepratorIndex(i)));
        maxPlus2Straight(i) = max(plusStraight(sepratorIndex(i):end));
        multFreAmpValue_straightPipe{i} = calcWaveFreAmplitude(temp,Fs,multFre,'freErr',1);

        X = straightPipeSection;
        calcDatas{dataCount,1} = sprintf('直管');
        calcDatas{dataCount,2} = X;
        calcDatas{dataCount,3} = plusStraight;
        calcDatas{dataCount,4} = multFreAmpValue_straightPipe{i}(1,:);
        calcDatas{dataCount,5} = multFreAmpValue_straightPipe{i}(2,:);
        calcDatas{dataCount,6} = multFreAmpValue_straightPipe{i}(3,:);
        calcDatas{dataCount,7} = maxPlus1Straight(i);
        calcDatas{dataCount,8} = maxPlus2Straight(i);
        dataCount = dataCount + 1;
    end
    
     %计算单一缓冲罐

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
        plusOV{i} = [plus1OV,plus2OV];
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
        
        multFreAmpValue_OV{i} = calcWaveFreAmplitude([pressure1OV,pressure2OV],Fs,multFre,'freErr',1);
        temp = eval(sprintf('%s(i)',IteratorValueName));
        calcDatas{dataCount,1} = sprintf('单一缓冲罐,%s:%g',IteratorValueName,temp);
        calcDatas{dataCount,2} = X;
        calcDatas{dataCount,3} = plusOV{i};
        calcDatas{dataCount,4} = multFreAmpValue_OV{i}(1,:);
        calcDatas{dataCount,5} = multFreAmpValue_OV{i}(2,:);
        calcDatas{dataCount,6} = multFreAmpValue_OV{i}(3,:);
        calcDatas{dataCount,7} = maxPlus1;
        calcDatas{dataCount,8} = maxPlus2;
        dataCount = dataCount + 1;

end
