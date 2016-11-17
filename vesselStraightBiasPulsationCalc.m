function [pressure1,pressure2] = vesselStraightBiasPulsationCalc(massFlowE,Frequency,time ...
,L1,L2,Lv,l,Dpipe,Dv,lv2,Dbias,sectionL1,sectionL2,varargin)
%缓冲罐入口顺接，出口错位的气流脉动计算
% Dbias 偏置管内插入缓冲罐的管径，如果偏置管没有内插如缓冲罐，Dbias为0
%   Detailed explanation goes here
%           |  L2
%        l  |     Lv    outlet
%   bias2___|_______________
%       |                   |  Dpipe
%       |lv2  V          lv1|———— L1  
%       |___________________| inlet
%           Dv              l      
%                       
%              
%容器的传递矩阵
pp=varargin;
a = nan;%声速


isDamping = 0;
isOpening = 1;
coeffFriction = nan;
meanFlowVelocity = nan;
isUseStaightPipe = 1;%使用直管理论代替缓冲罐，那么缓冲罐时相当于三个直管拼接
mach = nan;
notMach = 0;%强制不使用mach
while length(pp)>=2
    prop =pp{1};
    val=pp{2};
    pp=pp(3:end);
    switch lower(prop)
        case 'd' %管道直径
            Dpipe = val;
        case 'a'
        	a = val;
        case 'acousticvelocity'
        	a = val;
        case 'acoustic'
        	a = val;
        case 'isdamping' %是否包含阻尼
            isDamping = val;   
        case 'friction' %管道摩擦系数，计算阻尼系数时使用，如果输入是一个长度为2的向量，第一个代表直管的，第二个代表缓冲罐的
            coeffFriction = val;
        case 'coefffriction' %管道摩擦系数，计算阻尼系数时使用，如果输入是一个长度为2的向量，第一个代表直管的，第二个代表缓冲罐的
            coeffFriction = val;
        case 'meanflowvelocity' %平均流速，计算阻尼系数时使用，如果输入是一个长度为2的向量，第一个代表直管的，第二个代表缓冲罐的
            meanFlowVelocity = val;
        case 'flowvelocity' %平均流速，计算阻尼系数时使用,注意如果输入流速只有一个数值时，此流速代表缓冲罐的管道的流速，而不是缓冲罐里的流速
            meanFlowVelocity = val;
        case 'isusestaightpipe'
            isUseStaightPipe = val;%使用直管理论替代
        case 'usestaightpipe'
            isUseStaightPipe = val;
        case 'mach' %马赫数，加入马赫数将会使用带马赫数的公式计算
            mach = val;
        case 'm'
            mach = val;
        case 'isopening'
            isOpening = val;
        case 'notmach'
            notMach = val;
        otherwise
       		error('参数错误%s',prop);
    end
end
%如果用户没有定义k那么需要根据其他进行计算
if isnan(a)
    error('声速必须定义');
end
count = 1;
pressureE1 = [];
for i = 1:length(Frequency)
    f = Frequency(i);
    %最末端管道
    matrix_L2{count} = straightPipeTransferMatrix(L2,'f',f,'a',a,'d',Dpipe...
        ,'isDamping',isDamping,'coeffFriction',coeffFriction,'meanFlowVelocity',meanFlowVelocity...
        ,'m',mach,'notmach',notMach);
    matrix_Mv{count} = vesselStraightBiasTransferMatrix(Lv,l,lv2,Dbias ...
        ,'a',a,'d',Dpipe,'dv',Dv,'isDamping',isDamping,'coeffFriction',coeffFriction,'meanFlowVelocity',meanFlowVelocity,'f',f ...
        ,'isUseStaightPipe',isUseStaightPipe,'m',mach,'notmach',notMach);
    matrix_L1{count} = straightPipeTransferMatrix(L1,'f',f,'a',a,'D',Dpipe...
        ,'isDamping',isDamping,'coeffFriction',coeffFriction,'meanFlowVelocity',meanFlowVelocity...
        ,'m',mach,'notmach',notMach);
    matrix_total = matrix_L2{count} * matrix_Mv{count} * matrix_L1{count};
    A = matrix_total(1,1);
    B = matrix_total(1,2);
    C = matrix_total(2,1);
    D = matrix_total(2,2);

    if(isOpening)
        pressureE1(count) = ((-B/A)*massFlowE(count));
    else
        pressureE1(count) = ((-D/C)*massFlowE(count));
    end
    count = count + 1;
end

count = 1;
pressure1 = [];
if ~isempty(sectionL1)
    for len = sectionL1
        pressureEi = [];
        for i = 1:length(Frequency)
            f = Frequency(i);
            matrix_lx1 = straightPipeTransferMatrix(len,'f',f,'a',a,'D',Dpipe...
            ,'isDamping',isDamping,'coeffFriction',coeffFriction,'meanFlowVelocity',meanFlowVelocity...
            ,'m',mach,'notmach',notMach);
            pressureEi(i) = matrix_lx1(1,1)*pressureE1(i) + matrix_lx1(1,2)*massFlowE(i);
        end       
        pressure1(:,count) = changToWave(pressureEi,Frequency,time);
        count = count + 1;
    end
end

count = 1;
pressure2 = [];
if ~isempty(sectionL2)
    for len = sectionL2
        pressureEi = [];
        for i = 1:length(Frequency)
            f = Frequency(i);
            matrix_lx2 = straightPipeTransferMatrix(len,'f',f,'a',a,'D',Dpipe...
            ,'isDamping',isDamping,'coeffFriction',coeffFriction,'meanFlowVelocity',meanFlowVelocity...
            ,'m',mach,'notmach',notMach);
            matrix_Xl2_total = matrix_lx2  * matrix_Mv{i} * matrix_L1{i};
        
            pressureEi(i) = matrix_Xl2_total(1,1)*pressureE1(i) + matrix_Xl2_total(1,2)*massFlowE(i);
        end
        pressure2(:,count) = changToWave(pressureEi,Frequency,time);
        count = count + 1;
    end
end
end