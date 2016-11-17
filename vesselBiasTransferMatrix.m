function [M,k] = vesselBiasTransferMatrix(Lv,l,lv1,lv2,Dbias,varargin )
%缓冲罐错位的气流脉动计算
% Dbias 偏置管内插入缓冲罐的管径，如果偏置管没有内插如缓冲罐，Dbias为0
%   Detailed explanation goes here
%           |  L2
%        l  |     Lv    outlet
%   bias2___|_______________
%       |                   |
%       |lv2  V          lv1|  Dv
%       |___________________|
%                    l  |   bias1  
%                       |
%              inlet:   | L1 Dpipe 
%容器的传递矩阵
pp=varargin;
k = nan;
oumiga = nan;
f = nan;
a = nan;%声速
S = nan;
Sv = nan;

Dpipe = nan;
Dvessel = nan;
isDamping = 0;
coeffDamping = nan;
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
        case 's' %截面
            S = val;
        case 'd' %管道直径
            Dpipe = val;
            S = (pi.*Dpipe.^2)./4;
        case 'sv' %h缓冲罐截面
            Sv = val;
        case 'dv' %h缓冲罐截面
            Dvessel = val;
            Sv = (pi.*Dvessel.^2)./4;
        case 'k'
        	k = val;
        case 'oumiga'
        	oumiga = val;
        case 'f'
        	f = val;
        case 'a'
        	a = val;
        case 'acousticvelocity'
        	a = val;
        case 'acoustic'
        	a = val;
        case 'isdamping' %是否包含阻尼
            isDamping = val;   
        case 'coeffdamping' %阻尼系数，是一个长度为2的向量，第一个代表直管的，第二个代表缓冲罐的
            coeffDamping = val;
        case 'damping' %阻尼系数，是一个长度为2的向量，第一个代表直管的，第二个代表缓冲罐的
            coeffDamping = val;
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
        otherwise
       		error('参数错误%s',prop);
    end
end
%如果用户没有定义k那么需要根据其他进行计算
if isnan(a)
    error('声速必须定义');
end
if isnan(S)
    error('连接缓冲罐的管道截面积需要定义‘S’,或定义‘d’直径');
end
if isnan(Sv)
    error('缓冲罐的截面积需要定义‘Sv’,或定义‘dv’直径');
end
if isnan(k)
	if isnan(oumiga)
		if isnan(f)
			error('在没有输入k时，至少需要定义oumiga,f,acoustic中的两个');
		else
			oumiga = 2.*f.*pi;
		end
	end
	k = oumiga./a;
end
%流速修正
if ~isnan(meanFlowVelocity)
    if 1 == length(meanFlowVelocity)
        mfvVessel = meanFlowVelocity.*S./Sv;
        meanFlowVelocity = [meanFlowVelocity,mfvVessel];
    end
end
if isDamping
    if isnan(coeffDamping)
        if isnan(coeffFriction)
            error('若需要计算阻尼，且没有定义阻尼系数，需要定义“coeffFriction”管道摩擦系数');
        end
        if isnan(meanFlowVelocity)
            error('若需要计算阻尼，且没有定义阻尼系数，需要定义“meanFlowVelocity”平均流速');
        end

        if isnan(Dpipe)
            Dpipe = (4.*S./pi).^0.5;
        end
        if isnan(Dvessel)
            Dvessel = (4.*Sv./pi).^0.5;
        end
        Dtemp = [Dpipe,Dvessel];
        coeffDamping = (4.*coeffFriction.*meanFlowVelocity./Dtemp)./(2.*a);       
    end
    if isUseStaightPipe
        if isnan(meanFlowVelocity)
            error('使用直管等效，且有阻尼情况，需要定义平均流速');
        end
    end
end
if ~notMach%允许使用马赫
    if isnan(mach)
        if ~isnan(meanFlowVelocity)
            mach = meanFlowVelocity./a;
        end
    elseif(length(mach) == 1)
          mach(2) = meanFlowVelocity(2)/a;
    end
else
    mach = nan;
end
optMachStraight.notMach = notMach;
optMachStraight.mach = mach(1);
optMachVessel.notMach = notMach;
if(notMach)
    if(length(mach) == 1)
        optMachVessel.mach = mach(1);
    end
else
    optMachVessel.mach = mach(2);
end

%前管道传递矩阵
M1 = straightPipeTransferMatrix(l,'k',k,'S',S,'a',a,...
     'isDamping',isDamping,'coeffDamping',coeffDamping(1)...
     ,'mach',optMachStraight.mach,'notmach',optMachStraight.notMach);
%缓冲罐传递矩阵
%阻尼参数设置
optDamp.isDamping = 0;
optDamp.coeffDamping = coeffDamping(2);
optDamp.meanFlowVelocity = meanFlowVelocity(2);
Mv = vesselMatrix(isUseStaightPipe,Lv,lv1,lv2,k,Dvessel,Dbias,a,optDamp,optMachVessel);
%后管道传递矩阵
M2 = straightPipeTransferMatrix(l,'k',k,'S',S,'a',a,...
     'isDamping',isDamping,'coeffDamping',coeffDamping(1)...
     ,'mach',optMachStraight.mach,'notmach',optMachStraight.notMach);
    
M = M2*Mv*M1;
end

function Mv = vesselMatrix(isUseStaightPipe,Lv,lv1,lv2,k,Dv,Dbias,a,optDamping,optMach)
    if ~isstruct(optDamping)
        if isnan(optDamping)
            optDamping.isDamping = 0;
            optDamping.coeffDamping = 0;
            optDamping.meanFlowVelocity = 10;
        end
    end
    if ~isstruct(optMach)
        if isnan(optMach)
            optMach.notMach = 1;
            optMach.mach = 0;
        end
    end
    if isUseStaightPipe%使用直管理论
        ML = straightPipeTransferMatrix(Lv-lv1-lv2,'k',k,'D',Dv,'a',a,...
                'isDamping',optDamping.isDamping,'coeffDamping',optDamping.coeffDamping...
                ,'mach',optMach.mach,'notmach',optMach.notMach);%直管传递矩阵
        A = 0;
        B = 0;
        %考虑马赫数和不考虑马赫数对变径的传递矩阵有影响
        if optMach.notMach%不考虑马赫数
            Sv = pi.*Dv.^2./4;
            A = optDamping.coeffDamping*optDamping.meanFlowVelocity/Sv;%不考虑马赫数变径传递矩阵的右上角项
            B = A;%不考虑马赫数变径传递矩阵是相等的
            Mv = [1,B;0,1]*ML*[1,A;0,1];%直管加两个变径
            return;
        end
        %内插管腔体传递矩阵-左边
        innerLM = innerPipeCavityTransferMatrix(Dv,Dbias,lv2,'a',a,'k',k);
        %内插管腔体传递矩阵-右边
        innerRM = innerPipeCavityTransferMatrix(Dv,Dbias,lv1,'a',a,'k',k);
        Mv = innerLM*ML*innerRM;
        return;
    end
    %使用容积传递矩阵
    V = Sv*L;
    Mv = [1,0;-1i.*V.*k./a,1];
end