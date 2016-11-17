function M = vesselHavePerforatedPipeInletTransferMatrix(Dpipe,Dv,l,Lv,...
    lc,lv,dp,n,Lin,Lout,V,Din,varargin)

%计算管容管容的脉�?
%  长度 L1     l    Lv      l    L2  
%              _____________        
%             |   | dp(n)   |
%             |_ _|_ _ lc   |      
%  -----------|_ _ _ _Din   |----------
%             |lin|lout     |
%             |___|_________|       
% 直径 Dpipe        Dv           Dpipe 
%              

%当孔管与入口管紧密连接，上述结构等效为亥姆霍兹共鸣器与膨�?��的串�?
% lin 孔管入口长（共鸣器直径）
% lout孔管出口段长
% lc 孔管壁厚
% dp 孔管�?��孔的孔径
% n  孔管�?��个数
% Dp 孔管总孔径dp*n
% V  亥姆霍兹共鸣器体�?
% lv 共鸣器长
% 
%
%  长度 L1       l    Lv      l    L2  
%                   _________        
%                  | dp(n)   |
%                  |_ _ lc   |      
%  ----------------|_ _Din   |----------
%          lc | |  |lout     |
%            —| |�?|_________|       
%           |  Dp |
%        lv |  V  |
%           |     |
%            —�?—�?�?
%             lin
% 直径 Dpipe           Dv         Dpipe 

pp=varargin;
k = nan;
oumiga = nan;
f = nan;
a = nan;%声�?
isDamping = 1;%默认使用阻尼
coeffDamping = nan;
coeffFriction = nan;
meanFlowVelocity = nan;
isUseStaightPipe = 1;%使用直管理论代替缓冲罐，那么缓冲罐时相当于三个直管拼�?
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
        case 'coeffdamping' %阻尼系数，是�?��长度�?的向量，第一个代表直管的，第二个代表缓冲罐的
            coeffDamping = val;
        case 'damping' %阻尼系数，是�?��长度�?的向量，第一个代表直管的，第二个代表缓冲罐的
            coeffDamping = val;
        case 'friction' %管道摩擦系数，计算阻尼系数时使用，如果输入是�?��长度�?的向量，第一个代表直管的，第二个代表缓冲罐的
            coeffFriction = val;
        case 'coefffriction' %管道摩擦系数，计算阻尼系数时使用，如果输入是�?��长度�?的向量，第一个代表直管的，第二个代表缓冲罐的
            coeffFriction = val;
        case 'meanflowvelocity' %平均流�?，计算阻尼系数时使用，如果输入是�?��长度�?的向量，第一个代表直管的，第二个代表缓冲罐的
            meanFlowVelocity = val;
        case 'flowvelocity' %平均流�?，计算阻尼系数时使用,注意如果输入流�?只有�?��数�?时，此流速代表缓冲罐的管道的流�?，�?不是缓冲罐里的流�?
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
%如果用户没有定义k那么�?��根据其他进行计算
if isnan(a)
    error('声�?必须定义');
end
if isnan(k)
	if isnan(oumiga)
		if isnan(f)
			error('在没有输入k时，至少�?��定义oumiga,f,acoustic中的两个');
		else
			oumiga = 2.*f.*pi;
		end
	end
	k = oumiga./a;
end
%流�?修正
S = pi .* Dpipe.^2 ./ 4;
Sv = pi .* Dv.^2 ./ 4;
mfvVessel = nan;
if ~isnan(meanFlowVelocity)
    if 1 == length(meanFlowVelocity)
        mfvVessel = meanFlowVelocity.*S./Sv;
        meanFlowVelocity = [meanFlowVelocity,mfvVessel];
    end
else 
    error(['�?��指定流�?，流速是管道进入缓冲罐时的流速，',...
    '若需要指定缓冲罐流�?，可以使用一个含有两个元素的向量[pipe，vessel]']);
end
mfvVessel = meanFlowVelocity(2);
if isDamping
    if isnan(coeffDamping)
        if isnan(coeffFriction)
            error('若需要计算阻尼，且没有定义阻尼系数，�?��定义“coeffFriction”管道摩擦系�?);
        end
        if isnan(meanFlowVelocity)
            error('若需要计算阻尼，且没有定义阻尼系数，�?��定义“meanFlowVelocity”平均流�?);
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
            error('使用直管等效，且有阻尼情况，�?��定义平均流�?');
        end
    end
end
if length(meanFlowVelocity) < 2
    if isnan(coeffDamping) < 2
        error('参数不全，至少meanFlowVelocity，或coeffDamping�?��定义�?��');
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
optDamping.meanFlowVelocity = meanFlowVelocity(2);%缓冲罐平均流�?
Mv = havePerforatedPipeInletTransferMatrix(a,k,Dpipe,Dv,l,Lv,...
    lc,lv,dp,n,Lin,Lout,V,Din,optDamping,optMachVessel,mfvVessel);
M = M2 * Mv * M1;
end
%这里都是用直管等�?
function M = havePerforatedPipeInletTransferMatrix(a,k,Dpipe,Dv,l,Lv ...
    ,lc,lv,dp,n,Lin,Lout,V,Din,optDamping,optMach,mfvVessel)
%计算管容管容的脉�?
%  长度 L1     l    Lv      l    L2  
%              _____________        
%             |   | dp(n)   |
%             |_ _|_ _ lc   |      
%  -----------|_ _ _ _Din   |----------
%             |lin|lout     |
%             |___|_____Lv1_|       
% 直径 Dpipe        Dv           Dpipe 
%              

%当孔管与入口管紧密连接，上述结构等效为亥姆霍兹共鸣器与膨�?��的串�?
% lin 孔管入口长（共鸣器直径）
% lout孔管出口段长
% lc 孔管壁厚
% dp 孔管�?��孔的孔径
% n  孔管�?��个数
% Dp 孔管总孔径dp*n
% V  亥姆霍兹共鸣器体�?
% lv 共鸣器长
% 
%
%  长度 L1 l       Lv        l    L2  
%                   _________        
%                  | dp(n)   |
%                  |_ _ lc   |      
%  ----------------|_ _Din   |----------
%          lc | |  |lout     |
%            —| |�?|_________|       
%           |  Dp |       Lv1
%        lv |  V  |
%           |     |
%            —�?—�?�?
%             lin
% 直径 Dpipe           Dv         Dpipe
%section  1                        2   缓冲罐分的两个区�?
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
    
    Lv1 = Lv - Lin - Lout;%缓冲罐内插孔管无交接的区域长�?
    if ((Lv1 < 0))
        error('长度尺寸有误');
    end
    Mv1 = straightPipeTransferMatrix(Lv1,'k',k,'d',Dv,'a',a,...
                'isDamping',optDamping.isDamping,'coeffDamping',optDamping.coeffDamping...
                ,'mach',optMach.mach,'notmach',optMach.notMach);
    
    %急�?变径的传递矩�?
    Sv = pi.* Dv.^2 ./ 4;
    Spipe = pi.* Dpipe.^2 ./ 4;
    %SinnerPipe = pi .* Din.^2 ./ 4;
    %mfvInnerPipe = mfvVessel .* Sv ./ SinnerPipe;%内插管的平均流�?
    LM = sudEnlargeTransferMatrix(Spipe,Sv,a,'coeffdamping',optDamping.coeffDamping,'mach',optMach.mach,'notMach',optMach.notMach);
    %LM = sudReduceTransferMatrix(Spipe,Sv,1,a,'coeffdamping',optDamping.coeffDamping,'mach',optMach.mach);
    RM = sudReduceTransferMatrix(Sv,Spipe,a,'coeffdamping',optDamping.coeffDamping,'mach',optMach.mach,'notMach',optMach.notMach);
    %RM = sudReduceTransferMatrix(Sv,Spipe,0,a,'coeffdamping',optDamping.coeffDamping,'mach',optMach.mach);
    %内插管腔体传递矩�?左边
    innerLM = innerPipeCavityTransferMatrix(Dv,Din,Lout,'a',a,'k',k);
    %内插管的管道传�?矩阵
    innerPipeDampingCoeff = Dv^3 / Din^3 * optDamping.coeffDamping;%阻尼系数的传�?
    innerPM = straightPipeTransferMatrix(Lout,'k',k,'d',Din,'a',a,...
                'isDamping',optDamping.isDamping,'coeffDamping',innerPipeDampingCoeff...
                ,'mach',optMach.mach,'notmach',optMach.notMach);
    %innerPML = [1,0;0,1]; 
    %亥姆霍兹共鸣器传递矩�?
    Dp = dp.*n;
    HM = HelmholtzResonatorTransferMatrix(V,lv,lc,dp,n,'a',a,'k',k);
    M = RM * Mv1 * innerLM * innerPM * HM;
end