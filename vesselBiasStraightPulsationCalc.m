function [pressure1,pressure2] = vesselBiasStraightPulsationCalc(massFlowE,Frequency,time ...
,L1,L2,Lv,l,Dpipe,Dv,lv1,Dbias,sectionL1,sectionL2,varargin)
%缓冲罐入口错位，出口顺接的气流脉动计算
% Dbias 偏置管内插入缓冲罐的管径，如果偏置管没有内插如缓冲罐，Dbias为0
%   Detailed explanation goes here
%   inlet   |  L1
%        l  |     Lv    
%   bias2___|_______________
%       |                   |  Dpipe
%       |lv1  V          lv2|———— L2  
%       |___________________| outlet
%           Dv              l     
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
outletPressure = nan;%指定出口压力
outletMassFlow = nan;%指定出口质量流量
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
        case 'notmach'
            notMach = val;
        case 'isopening'
            isOpening = val;
        case 'outletpressure'
            outletPressure = val;
        case 'outletmassflow'
            outletMassFlow = val;
        otherwise
       		error('参数错误%s',prop);
    end
end
if isnan(outletPressure)
    if isOpening
        outletPressure = 0;
    end
end
if isnan(outletMassFlow)
    if ~isOpening
        outletMassFlow = 0;
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
    matrix_Mv{count} = vesselBiasStraightTransferMatrix(Lv,l,lv1,Dbias ...
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

    
    if ~isnan(outletPressure)
        if length(outletPressure) > 1
            pressureE1(count) = (outletPressure(count)-(B*massFlowE(count)))/A;
        else
            pressureE1(count) = (outletPressure-(B*massFlowE(count)))/A;
        end
    elseif ~isnan(outletMassFlow)
        if length(outletMassFlow) > 1
            pressureE1(count) = (outletMassFlow(count)-(D*massFlowE(count)))/C;
        else
            pressureE1(count) = (outletMassFlow-(D*massFlowE(count)))/C;
        end
    end
%     if(isOpening)
%         pressureE1(count) = ((-B/A)*massFlowE(count));
%     else
% %         m2 = 50-1i.*50;%管道末端质量流量
% %         pressureE1(count) = m2+((-D/C)*massFlowE(count));
% %         p2 = 50+1000i;
% %         pressureE1(count) = p2+((-B/A)*massFlowE(count));
%             pressureE1(count) = ((-D/C)*massFlowE(count));    
%     end
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