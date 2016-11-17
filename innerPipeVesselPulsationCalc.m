function [pressure1,pressure2] = innerPipeVesselPulsationCalc(massFlowE,Frequency,time ...
    ,L1,L2,Lv,l,Dpipe,Dv,Din,Linner,Lin,Lout...
     ,sectionL1,sectionL2,varargin)
%计算管容管容的脉动
%  长度 L1     l    Lv      l    L2  
%              ______________        
%             |     _|_  lout|      
%  -----------| lin _ _Din   |----------
%             |______|_______|       
% 直径 Dpipe       Dv       Dpipe 
%             |Linner| 
pp=varargin;
k = nan;
oumiga = nan;
a = 345;%声速

isDamping = 1;
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
            
        % case 'sv' %h缓冲罐截面
        %     Sv = val;
        % case 'dv' %h缓冲罐截面
        %     Dvessel = val;
            
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
        case 'isusestaightpipe'
            isUseStaightPipe = val;%使用直管理论替代
        case 'usestaightpipe'
            isUseStaightPipe = val;
        case 'm'
            mach = val;
        case 'notmach' %强制用马赫数计算设定
            notMach = val;
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
        ,'m',mach,'notmach',notMach);
    matrix_Mv{count} = vesselHaveInnerPipeTransferMatrix(Lv,l,Dv,Dpipe,Linner,Lin,Lout,Din ...
        ,'a',a,'isDamping',isDamping,'coeffFriction',coeffFriction,'meanFlowVelocity',meanFlowVelocity,'f',f,'a',a ...
        ,'isUseStaightPipe',isUseStaightPipe,'m',mach,'notmach',notMach);
    matrix_L1{count} = straightPipeTransferMatrix(L1,'f',f,'a',a,'D',Dpipe...
        ,'isDamping',isDamping,'coeffFriction',coeffFriction,'meanFlowVelocity',meanFlowVelocity...
        ,'m',mach,'notmach',notMach);
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