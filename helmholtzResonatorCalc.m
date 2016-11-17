function [pressure1,pressure2] = helmholtzResonatorCalc(massFlowE,Frequency,time ...
    ,L1,L2,Dpipe,dp,n,lc,V,lv ...
     ,sectionL1,sectionL2,varargin)
%亥姆霍兹共鸣器传递矩阵，将多个开孔的孔管等效为亥姆霍兹共鸣器时，此代码认为多个开孔为一个总的开口直径Dp，不单独进行亥姆霍兹共鸣器的并联计算
% lv 共鸣器长
% lc 共鸣器连接管长
% Dp 共鸣器连接管直径             
%       __________                
%      |          |                   
%      |    V     | lv
%      |___    ___|     
%  L1      |  | lc    L2    
% _________|Dp|__________                  
% _______________________Dpipe  
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
sigma = 0;
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
        case 'm'
            mach = val;
        case 'notmach' %强制用马赫数计算设定
            notMach = val;
        case 'sigma'
            sigma = val;%修正系数
        case 'k' %波数
        	k = val;
        case 'oumiga' %圆频率
        	oumiga = val;
        otherwise
            error('参数错误%s',prop);
    end
end
if isnan(a)
    error('声速必须定义');
end

count = 1;
pressureE1 = [];
for i = 1:length(Frequency)
    f = Frequency(i);
    %最末端管道
    matrix_L2{count} = straightPipeTransferMatrix(L2,'f',f,'a',a,'D',Dpipe...
        ,'isDamping',isDamping,'coeffFriction',coeffFriction,'meanFlowVelocity',meanFlowVelocity...
        ,'m',mach,'notmach',notMach,'coeffDamping',coeffDamping,'k',k,'oumiga',oumiga);
    matrix_Mv{count} = HelmholtzResonatorTransferMatrix(V,lv,lc,dp,n ...
        ,'a',a,'f',f,'sigma',sigma,'k',k,'oumiga',oumiga);
    matrix_L1{count} = straightPipeTransferMatrix(L1,'f',f,'a',a,'D',Dpipe...
        ,'isDamping',isDamping,'coeffFriction',coeffFriction,'meanFlowVelocity',meanFlowVelocity...
        ,'m',mach,'notmach',notMach,'coeffDamping',coeffDamping,'k',k,'oumiga',oumiga);
    matrix_total = matrix_L2{count} * matrix_Mv{count} * matrix_L1{count};
    A = matrix_total(1,1);
    B = matrix_total(1,2);
    pressureE1(count) = ((-B/A)*massFlowE(count));
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