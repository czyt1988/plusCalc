function M = vesselHavePerfPipeInletOrifTransferMatrix(Dpipe,Dv,l,Lv,...
    lc,lv,Dp,Lin,Lout,V,Din,lo,varargin)

%缓冲罐内入口连接孔管，加孔板结构?
%      L1     l        Lv       l    L2  
%              _________________        
%             |   | dp(n)    |lo|
%             |_ _|_ _ lc    |  |      
%  -----------|_ _ _ _   Din    |----------
%             |Lin|Lout      |  |
%             |___|__________|__|       
%    Dpipe            Dv           Dpipe 
%              
%
% Lin 内插孔管入口段长度 亥姆霍兹共鸣器直径
% Lout内插孔管出口段长度
% lc 孔管壁厚
% dp 孔管每一个孔孔径
% n  开孔个数
% Dp 等效亥姆霍兹共鸣器小径dp*n
% V  亥姆霍兹共鸣器体积
% lv 亥姆霍兹共鸣器长度
% lo 孔板距离缓冲罐出口的距离
%
%       L1  l    Lv          l    L2  
%                   _____________        
%                  | dp(n)    |lo|
%                  |_ _ lc    |  |      
%  ----------------|_ _   Din    |----------
%          lc_| |_ |Lout      |  |
%           |     ||__________|__|       
%           |  Dp |
%        lv |  V  |
%           |     |
%           |_____|
%             Lin
%  Dpipe           Dv               Dpipe

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
        case 'coeffdamping' %阻尼系数，
            coeffDamping = val;
        case 'damping' %阻尼系数，
            coeffDamping = val;
        case 'friction' %管道摩擦系数，计算阻尼系数时使用
            coeffFriction = val;
        case 'coefffriction' %管道摩擦系数，计算阻尼系数时使用
            coeffFriction = val;
        case 'meanflowvelocity' %平均流
            meanFlowVelocity = val;
        case 'flowvelocity' %平均流
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

if isnan(a)
    error('声速必须定义');
end
if isnan(k)
	if isnan(oumiga)
		if isnan(f)
			error('在没有输入k时，至少定义oumiga,f,acoustic中的两个');
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
    error(['需指定流速，流速是管道进入缓冲罐时的流速，',...
    '若需要指定缓冲罐流速，可以使用一个含有两个元素的向量[pipe，vessel]']);
end
mfvVessel = meanFlowVelocity(2);
if isDamping
    if isnan(coeffDamping)
        if isnan(coeffFriction)
            error('若需要计算阻尼，且没有定义阻尼系数，需定义“coeffFriction”管道摩擦系数');
        end
        if isnan(meanFlowVelocity)
            error('若需要计算阻尼，且没有定义阻尼系数，需定义“meanFlowVelocity”平均流速');
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
            error('使用直管等效，且有阻尼情况，需定义平均流速');
        end
    end
end
if length(meanFlowVelocity) < 2
    if isnan(coeffDamping) < 2
        error('参数不全，至少meanFlowVelocity，或coeffDamping');
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

M1 = straightPipeTransferMatrix(l,'k',k,'d',Dpipe,'a',a...
      ,'isDamping',isDamping,'coeffDamping',coeffDamping(1) ...
        ,'mach',optMachStraight.mach,'notmach',optMachStraight.notMach);
M2 = straightPipeTransferMatrix(l,'k',k,'d',Dpipe,'a',a...
      ,'isDamping',isDamping,'coeffDamping',coeffDamping(1) ...
        ,'mach',optMachStraight.mach,'notmach',optMachStraight.notMach);
optDamping.isDamping = isDamping;
optDamping.coeffDamping = coeffDamping(2);%缓冲罐的阻尼系数
optDamping.meanFlowVelocity = meanFlowVelocity(2);%缓冲罐平均流???
Mv = havePerforatedPipeInletOrifTransferMatrix(a,k,Dpipe,Dv,Lv,...
    lc,lv,Dp,Lin,Lout,V,Din,lo,optDamping,optMachVessel,mfvVessel);
M = M2 * Mv * M1;
end
%这里都是用直管等???
function M = havePerforatedPipeInletOrifTransferMatrix(a,k,Dpipe,Dv,Lv ...
    ,lc,lv,Dp,Lin,Lout,V,Din,lo,optDamping,optMach,mfvVessel)
%缓冲罐内入口连接孔管，加孔板结构?
%      L1     l        Lv       l    L2  
%              _________________        
%             |   | dp(n)    |lo|
%             |_ _|_ _ lc    |  |      
%  -----------|_ _ _ _   Din    |----------
%             |Lin|Lout      |  |
%             |___|__________|__|       
%    Dpipe            Dv           Dpipe 
%              
%
% Lin 内插孔管入口段长度 亥姆霍兹共鸣器直径
% Lout内插孔管出口段长度
% lc 孔管壁厚
% dp 孔管每一个孔孔径
% n  开孔个数
% Dp 等效亥姆霍兹共鸣器小径dp*n
% V  亥姆霍兹共鸣器体积
% lv 亥姆霍兹共鸣器长度
% 
%
%       L1  l    Lv          l    L2  
%                   _____________        
%                  | dp(n)    |lo|
%                  |_ _ lc    |  |      
%  ----------------|_ _   Din    |----------
%          lc_| |_ |Lout      |  |
%           |     ||__________|__|       
%           |  Dp |
%        lv |  V  |
%           |     |
%           |_____|
%             Lin
%  Dpipe           Dv               Dpipe
%section  1                           2   缓冲罐分的两个区
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
    
    Lv1 = Lv - Lin - Lout-lo;%缓冲罐内插孔管无交接的区域长???
    if ((Lv1 < 0))
        error('长度尺寸有误');
    end
    Mv1 = straightPipeTransferMatrix(Lv1,'k',k,'d',Dv,'a',a,...
                'isDamping',optDamping.isDamping,'coeffDamping',optDamping.coeffDamping...
                ,'mach',optMach.mach,'notmach',optMach.notMach);
    Mv2 = straightPipeTransferMatrix(lo,'k',k,'d',Dv,'a',a,...
                'isDamping',optDamping.isDamping,'coeffDamping',optDamping.coeffDamping...
                ,'mach',optMach.mach,'notmach',optMach.notMach);
    
    %急??变径的传递矩???
    Sv = pi.* Dv.^2 ./ 4;
    Spipe = pi.* Dpipe.^2 ./ 4;
    %SinnerPipe = pi .* Din.^2 ./ 4;
    %mfvInnerPipe = mfvVessel .* Sv ./ SinnerPipe;%内插管的平均流??
    %LM = sudEnlargeTransferMatrix(Spipe,Sv,a,'coeffdamping',optDamping.coeffDamping,'mach',optMach.mach,'notMach',optMach.notMach);
    %LM = sudReduceTransferMatrix(Spipe,Sv,1,a,'coeffdamping',optDamping.coeffDamping,'mach',optMach.mach);
    RM = sudReduceTransferMatrix(Sv,Spipe,a,'coeffdamping',optDamping.coeffDamping,'mach',optMach.mach,'notMach',optMach.notMach);
    %RM = sudReduceTransferMatrix(Sv,Spipe,0,a,'coeffdamping',optDamping.coeffDamping,'mach',optMach.mach);
    %内插管腔体传递矩???左边
    innerLM = innerPipeCavityTransferMatrix(Dv,Din,Lout,'a',a,'k',k);
    %内插管的管道传??矩阵
    innerPipeDampingCoeff = Dv^3 / Din^3 * optDamping.coeffDamping;%阻尼系数的传???
    innerPM = straightPipeTransferMatrix(Lout,'k',k,'d',Din,'a',a,...
                'isDamping',optDamping.isDamping,'coeffDamping',innerPipeDampingCoeff...
                ,'mach',optMach.mach,'notmach',optMach.notMach);
    %innerPML = [1,0;0,1]; 
    %亥姆霍兹共鸣器传递矩???
    HM = HelmholtzResonatorTransferMatrix(V,lv,lc,Dp,'a',a,'k',k);
    %孔板矩阵
    OM = orificeTransferMatrix(Dv,Din,mfvVessel);
    M = RM * Mv2 * OM * Mv1 * innerLM * innerPM * HM;
end