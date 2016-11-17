function M = vesselHaveInnerEnlargeTransferMatrix(Lv,l,Dv,Dpipe,Linner ...
    ,Le,Dr1,Dr2,varargin)

%计算内插扩渐扩张管的脉动
%  长度 L1     l    Lv      l    L2  
%              ______________        
%             |      |Le(Lc) |
%             |     /        |      
%  -----------|Dr1     Dr2   |----------
%             |     \        |
%             |______|_______|       
% 直径 Dpipe      Dv           Dpipe 
%             |Linner| 
pp=varargin;
k = nan;
oumiga = nan;
f = nan;
a = nan;%声速
isDamping = 1;%默认使用阻尼
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
S = pi .* Dpipe.^2 ./ 4;
Sv = pi .* Dv.^2 ./ 4;
mfvVessel = nan;
if ~isnan(meanFlowVelocity)
    if 1 == length(meanFlowVelocity)
        mfvVessel = meanFlowVelocity.*S./Sv;
        meanFlowVelocity = [meanFlowVelocity,mfvVessel];
    end
else 
    error(['需要指定流速，流速是管道进入缓冲罐时的流速，',...
    '若需要指定缓冲罐流速，可以使用一个含有两个元素的向量[pipe，vessel]']);
end
mfvVessel = meanFlowVelocity(2);
if isDamping
    if isnan(coeffDamping)
        if isnan(coeffFriction)
            error('若需要计算阻尼，且没有定义阻尼系数，需要定义“coeffFriction”管道摩擦系数');
        end
        if isnan(meanFlowVelocity)
            error('若需要计算阻尼，且没有定义阻尼系数，需要定义“meanFlowVelocity”平均流速');
        end
        Dtemp = [Dpipe,Dv];
        coeffDamping = (4.*coeffFriction.*meanFlowVelocity./Dtemp)./(2.*a);       
    end
    if length(coeffDamping)<2
        %必须考虑两个
        coeffDamping(2) = (4.*coeffFriction.*mfvVessel./Dv)./(2.*a);
    end
    if isUseStaightPipe
        if isnan(meanFlowVelocity)
            error('使用直管等效，且有阻尼情况，需要定义平均流速');
        end
    end
end
if length(meanFlowVelocity) < 2
    if isnan(coeffDamping) < 2
        error('参数不全，至少meanFlowVelocity，或coeffDamping需要定义一个');
    end
end
if ~notMach%允许使用马赫
    if isnan(mach)%如果没有设置马赫数
        if ~isnan(meanFlowVelocity)%如果设置了平均流速
            mach = meanFlowVelocity./a;
        end
    elseif(length(mach) == 1)%如果设置了马赫数，但马赫数的长度为1，说明只设置了一个马赫数，缓冲罐的马赫数没有设置，那么计算缓冲罐的马赫数，之前程序把缓冲罐的平均流速设置了，所以直接用meanFlowVelocity(2)
          mach(2) = meanFlowVelocity(2)/a;
    end
else
    mach = nan;%不允许使用马赫，那么马赫数设置为nan
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

M1 = straightPipeTransferMatrix(l,'k',k,'d',Dpipe,'a',a...
      ,'isDamping',isDamping,'coeffDamping',coeffDamping(1) ...
        ,'mach',optMachStraight.mach,'notmach',optMachStraight.notMach);
M2 = straightPipeTransferMatrix(l,'k',k,'d',Dpipe,'a',a...
      ,'isDamping',isDamping,'coeffDamping',coeffDamping(1) ...
        ,'mach',optMachStraight.mach,'notmach',optMachStraight.notMach);
optDamping.isDamping = isDamping;
optDamping.coeffDamping = coeffDamping(2);%缓冲罐的阻尼系数
optDamping.meanFlowVelocity = meanFlowVelocity(2);%缓冲罐平均流速
Mv = haveInnerReduceTransferMatrix(a,k,Lv,Dv,Dpipe,Linner ...
    ,Le,Dr1,Dr2,optDamping,optMachVessel,mfvVessel);
M = M2 * Mv * M1;
end
%这里都是用直管等效
function M = haveInnerReduceTransferMatrix(a,k,Lv,Dv,Dpipe,Linner ...
    ,Le,Dr1,Dr2,optDamping,optMach,mfvVessel)
%  长度 L1     l    Lv      l    L2  
%              _____Le(Lc)___       
%             |      |       |
%             |     /        |      
%  -----------|Dr1     Dr2   |----------
%             |     \        |
%             |______|_______|       
% 直径 Dpipe      Dv           Dpipe 
%             |Linner| 
%section        1          2   缓冲罐分的两个区域
    if ~isstruct(optDamping)
        if isnan(optDamping)
            optDamping.isDamping = 0;
            optDamping.coeffDamping = 0;%注意，这个是缓冲罐的祝你系数
            optDamping.meanFlowVelocity = 10;
        end
    end
    if ~isstruct(optMach)
        if isnan(optMach)
            optMach.notMach = 1;
            optMach.mach = 0;
        end
    end
    
    Lv2 = Lv - Linner;
    if (Lv2 < 0)
        error('长度尺寸有误');
    end

    Mv1 = straightPipeTransferMatrix(Linner-Le,'k',k,'d',Dv,'a',a,...
                'isDamping',optDamping.isDamping,'coeffDamping',optDamping.coeffDamping...
                ,'mach',optMach.mach,'notmach',optMach.notMach);

    Mv2 = straightPipeTransferMatrix(Lv-Linner,'k',k,'d',Dv,'a',a,...
                'isDamping',optDamping.isDamping,'coeffDamping',optDamping.coeffDamping...
                ,'mach',optMach.mach,'notmach',optMach.notMach);
    %急速变径的传递矩阵
    Sv = pi.* Dv.^2 ./ 4;
    Spipe = pi.* Dpipe.^2 ./ 4;
    %SinnerPipe = pi .* Din.^2 ./ 4;
    %mfvInnerPipe = mfvVessel .* Sv ./ SinnerPipe;%内插管的平均流速
    %RM = sudReduceTransferMatrix(Spipe,Sv,1,a,'coeffdamping',optDamping.coeffDamping,'mach',optMach.mach);
    %LM = sudReduceTransferMatrix(Sv,Spipe,0,a,'coeffdamping',optDamping.coeffDamping,'mach',optMach.mach);
    %由渐缩管道引起的入口处急速变径的传递矩阵
    LM = sudEnlargeTransferMatrix(Spipe,Sv,a,'coeffdamping',optDamping.coeffDamping,'mach',optMach.mach,'notMach',optMach.notMach);
    RM = sudReduceTransferMatrix(Sv,Spipe,a,'coeffdamping',optDamping.coeffDamping,'mach',optMach.mach,'notMach',optMach.notMach);
    Sv = pi.* Dv.^2 ./ 4;
    Sr2 = pi.* Dr2.^2 ./ 4;
    RMr = sudEnlargeTransferMatrix(Sr2,Sv,a...
        ,'coeffdamping',optDamping.coeffDamping...
        ,'mach',optMach.mach,'notmach',optMach.notMach);
    %内插管腔体左边
    innerLMEnlarge = innerEnlargeCavityTransferMatrix(Dv,Dr1,Le,'a',a,'k',k);
    %内插逐渐缩颈管传递矩阵
    innerEnlarge = innerEnlargeTransferMatrix(Dr1,Dr2,Le,'a',a,'k',k);
    M = RM * Mv2 * RMr * innerEnlarge * innerLMEnlarge *  Mv1 * LM;
end