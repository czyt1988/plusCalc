function M = innerPerfPipeEndSplitTransferMatrix(Din,Dv,lb2,varargin)
%内插孔管出口部分，孔管出口处开放，缓冲罐全开,理论自推。本代码只针对lb2部分管内和腔内分流最后汇合，参考书76页的干涉式消声器
% n2 孔管穿孔数
% dp 孔管穿孔孔径
% Din 穿管管径
% lp2 穿孔部分长度
% k 波数
% M1 管内的马赫数
% R0 系数默认0.0055
% lc 孔管壁厚
% xSection2 孔管每圈孔的间距，从0开始算，x的长度为孔管孔的圈数+1，x的值是当前一圈孔和上一圈孔的距离，如果间距一样，那么x里的值都一样
% xSection2 =[0,1,1,1,1,1,1,1,1,1]
%下述传递矩阵只表示带孔部分
%      -----------------------------------|
%     |----- -- -- --dp----               |lc
%-----|   Din                             |------
%     |----- -- --n2-- ----               |
%      lb1_______lp2____lb2_______________|
%
% lp 孔管长度（从第一个开孔算起，到最后一个开孔结束，见书79页图4.5.1）

pp = varargin;
k = nan;
oumiga = nan;
f = nan;
a = 345;%声速
meanFlowVelocity = nan;
M1 = nan;%管道马赫数
M2 = nan;%缓冲罐马赫数
while length(pp)>=2
    prop =pp{1};
    val=pp{2};
    pp=pp(3:end);
    switch lower(prop)
        case 'k' %波数
        	k = val;
        case 'oumiga' %圆频率
        	oumiga = val;
        case 'f' %脉动频率
        	f = val;
        case 'a'
        	a = val;
        case 'acousticvelocity'
        	a = val;
        case 'acoustic'
        	a = val;
        case 'meanflowvelocity'
            meanFlowVelocity = val;
        case 'm1'
            M1 = val;
        case 'm2'
            M2 = val;
        case 'r0'
            R0 = val;
        otherwise
       		error('参数错误%s',prop);
    end
end
%如果用户没有定义k那么需要根据其他进行计算
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
Sp = pi*Din.^2./4;%孔管管径截面积
Sv = pi*Dv.^2./4;%缓冲罐截面积
Sv_inner = Sv-Sp;%内插管环绕腔体
if isnan(M1)
    if isnan(meanFlowVelocity)
        error('流速必须给定');
    end
    if isnan(a)
        error('声速必须给定');
    end
    M1 = meanFlowVelocity./a;
    M2 = M1 .* Sp ./ (Sv-Sp);
end


kc1 = k./(1-M1^2);
kc2 = k./(1-M2^2);
coeffA = Sp./tan(kc1.*lb2);
coeffB = Sv_inner./tan(kc2.*lb2);
coeffC = Sp.*exp(1i.*M1.*kc1.*lb2)./sin(kc1.*lb2);
coeffD = Sv_inner.*exp(1i.*M2.*kc2.*lb2)./sin(kc2.*lb2);
coeffE = Sp.*exp(-1i.*M1.*kc1.*lb2)./sin(kc1.*lb2);
coeffF = Sv_inner.*exp(-1i.*M2.*kc2.*lb2)./sin(kc2.*lb2);
coeffT(1,1) = (coeffA + coeffB)./(coeffC + coeffD);
coeffT(1,2) = 1i.*a./(coeffC + coeffD);
coeffT(2,1) = 1i./a.*(coeffE + coeffF)-1i./a.*(coeffA + coeffB)^2./(coeffC + coeffD);
coeffT(2,2) = (coeffA + coeffB)./(coeffC + coeffD);


coeffFM = coeffT(2,2)*coeffT(1,1)-coeffT(1,2)*coeffT(2,1);

A = coeffT(2,2)./coeffFM;
B = (-a./Sp)*coeffT(1,2)./coeffFM;
C = (Sp./a*coeffT(2,1))./(-coeffFM);
D = -coeffT(1,1)./(-coeffFM);
M = [A,B;C,D];%[p2,m2]孔管出口端（开孔部分）所在缓冲罐内的压力=M*[p1,m1]
end
 
