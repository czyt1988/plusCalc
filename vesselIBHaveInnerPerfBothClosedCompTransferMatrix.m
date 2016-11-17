function M = vesselIBHaveInnerPerfBothClosedCompTransferMatrix(Dpipe,Dv,l,Lv1,Lv2,...
    lc,dp1,dp2,lp1,lp2,n1,n2,la1,la2,lb1,lb2,Din,Dbias1,lv1,xSection1,xSection2,varargin)
%缓冲罐中间插入孔管,两端堵死，开孔个数不足以等效为亥姆霍兹共鸣器,缓冲罐入口偏置
%                 L1
%                     |
%                     |
%           l         |                      Lv2        l    L2  
%              _______|_________________________________        
%             |    dp1(n1)           |    dp2(n2)       |
%             |           ___ _ _ ___|___ _ _ ___ lc    |     
%             |          |___ _ _ ___ ___ _ _ ___|Din   |----------
%             |           la1 lp1 la2|lb1 lp2 lb2       |
%             |______________________|__________________|       
%                             Lin         Lout
%                       Lv1
%    Dpipe                       Dv                     Dpipe  
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
% xSection1，xSection2 孔管每圈孔的间距，从0开始算，x的长度为孔管孔的圈数+1，x的值是当前一圈孔和上一圈孔的距离，如果间距一样，那么x里的值都一样
pp=varargin;
k = nan;
oumiga = nan;
f = nan;
a = nan;%声速
isDamping = 1;%默认使用阻尼
coeffDamping = nan;
coeffFriction = nan;
meanFlowVelocity = nan;
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
Sv = pi .* Dv.^2 ./ 4;%缓冲罐截面积
Sp = pi*Din.^2./4;%孔管管径截面积
Sv_p = Sv-Sp;%去除孔管的缓冲罐截面积
Dv_inner = (4*Sv_p/pi).^0.5;%计算名义直径
mfvVessel = nan;
mfvInnerPipe = nan;
mfvVessel_Inner = nan;
if ~isnan(meanFlowVelocity)
    if 1 == length(meanFlowVelocity)
        mfvVessel = meanFlowVelocity*S/Sv;
        mfvVessel_Inner = meanFlowVelocity*S/Sv_p;
        mfvInnerPipe = meanFlowVelocity*S/Sp;
        meanFlowVelocity = [meanFlowVelocity,mfvVessel,mfvVessel_Inner,mfvInnerPipe];
    elseif 2 == length(meanFlowVelocity)
        mfvVessel = meanFlowVelocity(2);
        mfvVessel_Inner = meanFlowVelocity*S/Sv_p;
        mfvInnerPipe = meanFlowVelocity*S/Sp;
        meanFlowVelocity = [meanFlowVelocity,mfvVessel_Inner,mfvInnerPipe];
    elseif 3 == length(meanFlowVelocity)
        mfvVessel = meanFlowVelocity(2);
        mfvVessel_Inner = meanFlowVelocity(3);
        mfvInnerPipe = meanFlowVelocity*S/Sp;
        mfvInnerPipe = [meanFlowVelocity,mfvInnerPipe];
    elseif 4 == length(meanFlowVelocity)
        mfvVessel = meanFlowVelocity(2);
        mfvVessel_Inner = meanFlowVelocity(3);
        mfvInnerPipe = meanFlowVelocity(4);
    end
else 
    error(['需指定流速，流速是管道进入缓冲罐时的流速，',...
    '若需要指定缓冲罐流速，可以使用一个含有4个元素的向量[pipe，vessel,vessel_Inner,InnerPipe]']);
end

if isDamping
    if isnan(coeffDamping)
        if isnan(coeffFriction)
            error('若需要计算阻尼，且没有定义阻尼系数，需定义“coeffFriction”管道摩擦系数');
        end
        if isnan(meanFlowVelocity)
            error('若需要计算阻尼，且没有定义阻尼系数，需定义“meanFlowVelocity”平均流速');
        end
        if length(meanFlowVelocity) < 4
            error('“meanFlowVelocity”平均流速的长度过小，必须为4');
        end
        Dtemp = [Dpipe,Dv,Dv_inner,Din];
        coeffDamping = (4.*coeffFriction.*meanFlowVelocity./Dtemp)./(2.*a);       
    end
    if length(coeffDamping)<4
        %必须考虑4个
        coeffDamping(2) = (4.*coeffFriction.*mfvVessel./Dv)./(2.*a);
        coeffDamping(3) = (4.*coeffFriction.*mfvVessel./Dv_inner)./(2.*a);
        coeffDamping(4) = (4.*coeffFriction.*mfvVessel./Din)./(2.*a);
    end
end

mach = meanFlowVelocity./a;%马赫,mach(1)直管mach，mach(2)缓冲罐mach:mach(3)带内插管的缓冲罐的mach:mach(4)内插管的mach

optMach.notMach = notMach;
optMach.machStraight = mach(1);
optMach.machVessel = mach(2);
optMach.machVesselWithInnerPipe = mach(3);
optMach.machInnerPipe = mach(4);

optDamping.isDamping = isDamping;
optDamping.coeffDampStraight = coeffDamping(1);
optDamping.mfvStraight = meanFlowVelocity(1);

optDamping.coeffDampVessel = coeffDamping(2);%缓冲罐的阻尼系数
optDamping.mfvVessel = meanFlowVelocity(2);

optDamping.coeffDampVesselWithInnerPipe = coeffDamping(3);%缓冲罐的阻尼系数
optDamping.mfvVesselWithInnerPipe = meanFlowVelocity(3);

optDamping.coeffDampInnerPipe = coeffDamping(4);%缓冲罐的阻尼系数
optDamping.mfvInnerPipe = meanFlowVelocity(4);



M1 = straightPipeTransferMatrix(l,'k',k,'d',Dpipe,'a',a...
      ,'isDamping',isDamping,'coeffDamping',coeffDamping(1) ...
        ,'mach',optMach.machStraight,'notmach',optMach.notMach);
M2 = straightPipeTransferMatrix(l,'k',k,'d',Dpipe,'a',a...
      ,'isDamping',isDamping,'coeffDamping',coeffDamping(1) ...
        ,'mach',optMach.machStraight,'notmach',optMach.notMach);
    
Mv = haveInnerPerforatedPipeBCCompTransferMatrix(a,k,Dv,Dv_inner,Lv1,Lv2,...
    lc,dp1,dp2,lp1,lp2,n1,n2,la1,la2,lb1,lb2,Din,Dbias1,lv1,xSection1,xSection2,optDamping,optMach);
M = M2 * Mv * M1 ;

end
function M = haveInnerPerforatedPipeBCCompTransferMatrix(a,k,Dv,Dv_inner,Lv1,Lv2 ...
    ,lc,dp1,dp2,lp1,lp2,n1,n2,la1,la2,lb1,lb2,Din,Dbias1,lv1,xSection1,xSection2,optDamping,optMach)
%缓冲罐中间插入孔管,两端堵死，开孔个数不足以等效为亥姆霍兹共鸣器
%                 L1
%                     |
%                     |
%           l         |                      Lv2        l    L2  
%              _______|_________________________________        
%             |    dp1(n1)           |    dp2(n2)       |
%             |           ___ _ _ ___|___ _ _ ___ lc    |     
%             |          |___ _ _ ___ ___ _ _ ___|Din   |----------
%             |           la1 lp1 la2|lb1 lp2 lb2       |
%             |______________________|__________________|       
%                             Lin         Lout
%                       Lv1
%    Dpipe                       Dv                     Dpipe  
%section  1                              2   缓冲罐分的两个区              
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
% xSection1，xSection2 孔管每圈孔的间距，从0开始算，x的长度为孔管孔的圈数+1，x的值是当前一圈孔和上一圈孔的距离，如果间距一样，那么x里的值都一样

   if ~isstruct(optDamping)
        if isnan(optDamping)
            error('optDamping不能为空');
        end
    end
    if ~isstruct(optMach)
        if isnan(optMach)
            error('optMach不能为空');
        end
    end
    % 缓冲罐内不含孔管的部分
    Lin = la1 + lp1 + la2;
    Lout = lb1 + lp2 + lb2;
%     Lv1 = Lv/2 - Lin;%缓冲罐内插孔管无交接的区域长
%     Lv2 = Lv/2 - Lout;%缓冲罐内插孔管无交接的区域长
%    
%     if ((Lv1 < 0))
%         error('长度尺寸有误');
%     end
    Cav1 = Lv1 - Lin;%腔1非环形部分
    Cav2 = Lv2 - Lout;%腔2非环形部分
    Mv2 = straightPipeTransferMatrix(Cav2,'k',k,'d',Dv,'a',a,...
                'isDamping',optDamping.isDamping,'coeffDamping',optDamping.coeffDampVessel...
                ,'mach',optMach.machVessel,'notmach',optMach.notMach);
    %lb2处对应的缓冲罐腔室
    Mstr2 = straightPipeTransferMatrix(lb2,'k',k,'d',Dv_inner,'a',a,...
                'isDamping',optDamping.isDamping,'coeffDamping',optDamping.coeffDampVesselWithInnerPipe...
                ,'mach',optMach.machVesselWithInnerPipe,'notmach',optMach.notMach);
    %内插孔管出口段对应lp2
    Mp2 = innerPerfPipeOpenTransferMatrix(n2,dp2,Din,Dv,lp2,lc,lb1,lb2,xSection2...
        ,'k',k,'a',a,'meanflowvelocity',optDamping.mfvInnerPipe...
        ,'M1',optMach.machInnerPipe,'M2',optMach.machVesselWithInnerPipe);
    %内插孔管隔板区对应la2+lb1
    Mstr = straightPipeTransferMatrix(la2+lb1,'k',k,'d',Din,'a',a,...
                'isDamping',optDamping.isDamping,'coeffDamping',optDamping.coeffDampInnerPipe...
                ,'mach',optMach.machInnerPipe,'notmach',optMach.notMach);
    %内插孔管入口段对应lp1
    Mp1 = innerPerfPipeOpenInletClosedTransferMatrix(n1,dp1,Din,Dv,lp1,lc,la1,la2,xSection1...
        ,'k',k,'a',a,'meanflowvelocity',optDamping.mfvInnerPipe...
        ,'M1',optMach.machInnerPipe,'M2',optMach.machVesselWithInnerPipe);
    %内插孔管与进口管连接部分对应la1
    Mstr1 = straightPipeTransferMatrix(la1,'k',k,'d',Dv_inner,'a',a,...
                'isDamping',optDamping.isDamping,'coeffDamping',optDamping.coeffDampInnerPipe...
                ,'mach',optMach.machInnerPipe,'notmach',optMach.notMach);
    Mv1 = straightPipeTransferMatrix(Cav1-lv1,'k',k,'d',Dv,'a',a,...
            'isDamping',optDamping.isDamping,'coeffDamping',optDamping.coeffDampVessel...
            ,'mach',optMach.machVessel,'notmach',optMach.notMach);
    innerLM = innerPipeCavityTransferMatrix(Dv,Dbias1,lv1,'a',a,'k',k);
%     %入口孔管前端开放空腔对应
%     Mv1 = straightPipeTransferMatrix(Lv1,'k',k,'d',Dv,'a',a,...
%                 'isDamping',optDamping.isDamping,'coeffDamping',optDamping.coeffDampVessel...
%                 ,'mach',optMach.machVessel,'notmach',optMach.notMach);
      
%     %内插管腔体传递矩阵-右边对应lb1+lp2
%     Mca = innerPipeCavityTransferMatrix(Dv,Din,lb1+lp2,'a',a,'k',k);

    
   % M = Mv2 * Mstr2 * Mp2 * Mstr * Mp1 * Mstr1* Mv1;
   M = Mv2 * Mstr2 * Mp2 * Mstr * Mp1 * Mstr1 * Mv1 * innerLM;
%     M = innerLM * Mv2 * Mstr2 * Mp2 * Mstr * Mp1 * Mstr1* Mv1;
end