function calcDatas = funIteratorVesselPipeLinePlusCalc(time,massFlow,Fs,plusBaseFrequency,acousticVelocity...
    ,L1,L2,l,Dpipe,Dv,Lv,IteratorValueName,varargin)
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
% 'xx',
iteCount = eval(sprintf('length(%s)',IteratorValueName));

% 
% xSection1，xSection2 孔管每圈孔的间距，从0开始算，x的长度为孔管孔的圈数+1，x的值是当前一圈孔和上一圈孔的距离，如果间距一样，那么x里的值都一样
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
