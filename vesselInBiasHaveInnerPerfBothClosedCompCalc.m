function [pressure1,pressure2] = vesselInBiasHaveInnerPerfBothClosedCompCalc(massFlowE,Frequency,time ...
    ,L1,L2,Dpipe,vhpicStruct,sectionL1,sectionL2,varargin)
%缓冲罐中间插入孔管,两端堵死，开孔个数不足以等效为亥姆霍兹共鸣器,缓冲罐入口偏置
%                 L1
%                     |
%                     |
%           l         |                      Lv2        l    L2  
%              _______|_________________________________        
%             |    dp(n1)            |    dp(n2)        |
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
a = 345;%声速

isDamping = 1;
coeffDamping = nan;
coeffFriction = nan;
meanFlowVelocity = nan;
mach = nan;
notMach = 0;%强制不使用mach
isOpening = 1;
while length(pp)>=2
    prop =pp{1};
    val=pp{2};
    pp=pp(3:end);
    switch lower(prop)
        case 'a' %声速
            a = val; 
        case 'acousticvelocity' %声速
            a = val;
        case 'acoustic' %声速
            a = val;
        case 'isdamping' %是否包含阻尼
            isDamping = val;   
        case 'friction' %管道摩擦系数，计算阻尼系数时使用
            coeffFriction = val;
        case 'coefffriction' %管道摩擦系数，计算阻尼系数时使用
            coeffFriction = val;
        case 'meanflowvelocity' %平均流速，计算阻尼系数时使用
            meanFlowVelocity = val;
        case 'flowvelocity' %平均流速，计算阻尼系数时使用
            meanFlowVelocity = val;
        case 'mach' %马赫数，加入马赫数将会使用带马赫数的公式计算
            mach = val;
        case 'coeffdamping'
            coeffDamping = val;
        case 'm'
            mach = val;
        case 'notmach' %强制用马赫数计算设定
            notMach = val;
        case 'k' %波数
        	k = val;
        case 'oumiga' %圆频率
        	oumiga = val;
        case 'isopening'%管道末端是否为无反射端(开口)，如果为0，就是为闭口，无流量端
            isOpening = val;
        otherwise
            error('参数错误%s',prop);
    end
end
if isnan(a)
    error('声速必须定义');
end
% 
% L2
% Dpipe
% isDamping
% coeffFriction
% meanFlowVelocity
% mach
% notMach
% coeffDamping
% k
% oumiga


count = 1;
pressureE1 = [];
for i = 1:length(Frequency)
    f = Frequency(i);
    %最末端管道
    matrix_L2{count} = straightPipeTransferMatrix(L2,'f',f,'a',a,'D',Dpipe...
        ,'isDamping',isDamping,'coeffFriction',coeffFriction,'meanFlowVelocity',meanFlowVelocity...
        ,'m',mach,'notmach',notMach,'coeffDamping',coeffDamping,'k',k,'oumiga',oumiga);
    matrix_Mv{count} = vesselIBHaveInnerPerfBothClosedCompTransferMatrix(Dpipe,vhpicStruct.Dv,vhpicStruct.l,vhpicStruct.Lv1,vhpicStruct.Lv2...
        ,vhpicStruct.lc,vhpicStruct.dp1,vhpicStruct.dp2,vhpicStruct.lp1,vhpicStruct.lp2,vhpicStruct.n1,vhpicStruct.n2...
        ,vhpicStruct.la1,vhpicStruct.la2,vhpicStruct.lb1,vhpicStruct.lb2,vhpicStruct.Din...
        ,0,vhpicStruct.lv1...
        ,vhpicStruct.xSection1,vhpicStruct.xSection2...
        ,'f',f,'a',a,'k',k,'oumiga',oumiga...
        ,'coeffDamping',coeffDamping,'isDamping',isDamping,'coeffFriction',coeffFriction,'meanFlowVelocity',meanFlowVelocity...
        ,'notmach',notMach...
        );
    matrix_L1{count} = straightPipeTransferMatrix(L1,'f',f,'a',a,'D',Dpipe...
        ,'isDamping',isDamping,'coeffFriction',coeffFriction,'meanFlowVelocity',meanFlowVelocity...
        ,'m',mach,'notmach',notMach,'coeffDamping',coeffDamping,'k',k,'oumiga',oumiga);
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
        count2 = 1;
        pTemp = [];
        pressureEi = [];
        for i = 1:length(Frequency)
            f = Frequency(i);
            matrix_lx1 = straightPipeTransferMatrix(len,'f',f,'a',a,'D',Dpipe...
            ,'isDamping',isDamping,'coeffFriction',coeffFriction,'meanFlowVelocity',meanFlowVelocity...
            ,'m',mach,'notmach',notMach);
            pressureEi(count2) = matrix_lx1(1,1)*pressureE1(count2) + matrix_lx1(1,2)*massFlowE(count2);
            count2 = count2 + 1;
        end       
        pressure1(:,count) = changToWave(pressureEi,Frequency,time);
        count = count + 1;
    end
end

count = 1;
pressure2 = [];
if ~isempty(sectionL2)
    for len = sectionL2
        count2 = 1;
        pressureEi = [];
        for i = 1:length(Frequency)
            f = Frequency(i);
            matrix_lx2 = straightPipeTransferMatrix(len,'f',f,'a',a,'D',Dpipe...
            ,'isDamping',isDamping,'coeffFriction',coeffFriction,'meanFlowVelocity',meanFlowVelocity...
            ,'m',mach,'notmach',notMach);
            matrix_Xl2_total = matrix_lx2  * matrix_Mv{count2} * matrix_L1{count2};
        
            pressureEi(count2) = matrix_Xl2_total(1,1)*pressureE1(count2) + matrix_Xl2_total(1,2)*massFlowE(count2);
            count2 = count2 + 1;
        end
        pressure2(:,count) = changToWave(pressureEi,Frequency,time);
        count = count + 1;
    end
end

end