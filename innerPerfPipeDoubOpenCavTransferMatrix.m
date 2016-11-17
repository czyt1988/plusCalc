function M = innerPerfPipeDoubOpenCavTransferMatrix(n2,dp,Din,Dv,lp2,lc,lb1,lb2,varargin)
%孔管末端开口处腔体部分矩阵对应在lp2处
%内插孔管，孔管开孔数不是非常少，无法等价为亥姆霍兹共鸣器，且孔管出口开口，缓冲罐全开,理论来自文献《A study on the transmission loss of straight鈥恡hrough type reactive mufflers》
% n2 孔管穿孔数
% dp 孔管穿孔孔径
% Din 穿管管径
% lp2 穿孔部分长度
% k 波数
% M1 管内的马赫数
% R0 系数默认0.0055
% lc 孔管壁厚
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
R0 = 0.0055;% ?
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
% coeffD = calcD(n2,dp,Din,Dv,lp2,k,M1,M2,R0,lc,xSection2);
% coeffR = calcR(coeffD);
% coeffT = calcT(k,lb1,lb2,coeffR);

[Ta2,Tb2,Tc2,Td2] = calcT1(n2,dp,Din,Dv,lb1,lb2,lp2,k,M1,M2,R0,lc);
coeffFM = Td2*Ta2-Tb2*Tc2;



A = Td2./coeffFM;
B = (-a./(Sv-Sp))*Tb2./coeffFM;
C = (((Sv-Sp)./a)*Tc2)./(-coeffFM);
D = -((Sv-Sp)./Sp).*Ta2./(-coeffFM);
M = [A,B;C,D];%[p2,m2]孔管出口端（开孔部分）所在缓冲罐内的压力=M*[p1,m1]
end
 
function [Ta2,Tb2,Tc2,Td2] = calcT1(n2,dp,Din,Dv,lb1,lb2,lp2,k,M1,M2,R0,lc)
% n2 孔管穿孔数
% dp 孔管穿孔孔径
% Din 穿管管径
% lp2 穿孔部分长度
% k 波数
% M1 管内的马赫数
% R0 系数默认0.0055
% lc 孔管壁厚
    [a1,a2,a3,a4,a5,a6,a7,a8] = calcCoeffA(n2,dp,Din,Dv,lp2,k,M1,M2,R0,lc);
    [r1,r2,r3,r4] = calcCoeffR(a1,a2,a3,a4,a5,a6,a7,a8);
    beita = calcCoeffBeita(r1,r2,r3,r4);
    Fai = calcCoeffFai(beita,a1,a2,a3,a4);

    %计算T
    areaRatio = (Dv./Din)^2;

    AperB = calcCoeffAperB(k,areaRatio,M2,lb2);
    CperD = calcCoeffCperD(k,areaRatio,M2,lb2);
    Fi = calcCoeffFi(k,lb1,lp2,beita,Fai,AperB,CperD);
    [C31,C32,C41,C42] = calcCoeffC(k,lb1,lp2,AperB,CperD,Fai,Fi,beita);

    T01 = 1./beita(1) + C31./beita(3) + C41./beita(4);
    T02 = 1./beita(2) + C32./beita(3) + C42./beita(4);
    T03 = exp(beita(1).*lp2)./beita(1) + C31.*exp(beita(3).*lp2)./beita(3) + C41.*exp(beita(4).*lp2)./beita(4);
    T04 = exp(beita(2).*lp2)./beita(2) + C32.*exp(beita(3).*lp2)./beita(3) + C42.*exp(beita(4).*lp2)./beita(4);
    T05 = Fai(2,1).*exp(beita(1).*lp2)./beita(1) + Fai(2,3).*C31.*exp(beita(3).*lp2)./beita(3) + Fai(2,4).*C41.*exp(beita(4).*lp2)./beita(4);
    T06 = Fai(2,2).*exp(beita(2).*lp2)./beita(2) + Fai(2,3).*C32.*exp(beita(3).*lp2)./beita(3) + Fai(2,4).*C42.*exp(beita(4).*lp2)./beita(4);
    T07 = -(1./(1i.*k+M2.*beita(1)) + C31./(1i.*k+M2.*beita(3)) + C41./(1i.*k+M2.*beita(4)));
    T08 = -(1./(1i.*k+M2.*beita(2)) + C32./(1i.*k+M2.*beita(3)) + C42./(1i.*k+M2.*beita(4)));
    T09 = -(exp(beita(1).*lp2)./(1i.*k+M2.*beita(1)) + C31.*exp(beita(3).*lp2)./(1i.*k+M2.*beita(3)) + C41.*exp(beita(4).*lp2)./(1i.*k+M2.*beita(4)));
    T10 = -(exp(beita(2).*lp2)./(1i.*k+M2.*beita(2)) + C32.*exp(beita(3).*lp2)./(1i.*k+M2.*beita(3)) + C42.*exp(beita(4).*lp2)./(1i.*k+M2.*beita(4)));
    T11 = -(Fai(2,1).*exp(beita(1).*lp2) + Fai(2,3).*C31.*exp(beita(3).*lp2) + Fai(2,4).*C41.*exp(beita(4).*lp2))./(1i.*k);
    T12 = -(Fai(2,2).*exp(beita(2).*lp2) + Fai(2,3).*C32.*exp(beita(3).*lp2) + Fai(2,4).*C42.*exp(beita(4).*lp2))./(1i.*k);
    G2 = T05*T12-T06*T11;
    Ta2 = (T01*T12-T02*T11)./G2;
    Tb2 = (T02*T05-T01*T06)./G2;
    Tc2 = (T07*T12-T08*T11)./G2;
    Td2 = (T05*T08-T06*T07)./G2;

end




function [a1,a2,a3,a4,a5,a6,a7,a8] = calcCoeffA(n2,dp,Din,Dv,lp2,k,M1,M2,R0,lc)
% n2 孔管穿孔数
% dp 孔管穿孔孔径
% Din 穿管管径
% lp2 穿孔部分长度
% k 波数
% M1 管内的马赫数
% R0 系数默认0.0055
% lc 孔管壁厚

    bp2 = n2.*dp^2./(4.*Din.*lp2);%开孔率
    rp = dp/2;
    kr = k.*rp;
%     if kr <= M1 %M1为管内马赫数    
%         Cg = 12^(k.*rp./(M1-1));
%     else
%         Cg = 1;
%     end  
    sig0 = calcSigma0(dp);
    if M1 == 0
         kp = (R0+1i.*k.*(lc+sig0.*dp))./bp2;%kp为穿孔声阻抗率
    elseif M1>0
        kp = (R0 + 2.48.*(1-bp2).*(M1.^(1.04+40.*dp)) + 1i.*k.*(lc+(1-(1+398.34.*dp).*((1-bp2).^1.44).*(M1.^0.72)).*sig0.*dp))./bp2;%通过流
        %kp = (R0+0.48.*abs(k.*rp-M1)+1i.*k.*(lc+Cg.*sig0.*dp))./bp2;%掠过流
    else
        error('马赫数M1不能小于0');
    end
    
    a1 = -2.*M1.*(1i.*k+2./(Din.*kp))./(1-M1^2);
    a2 = 1.*(k^2-4.*1i.*k./(Din.*kp))./(1-M1^2);
    a3 = 4.*M1./(Din.*kp)./(1-M1^2);
    a4 = 4.*1i.*k./(Din.*kp)./(1-M1^2);
    a5 = M2.*4.*Din./((Dv.^2-(Din+2.*lc).^2).*kp)./(1-M2^2);
    a6 = 4.*1i.*k.*Din./((Dv.^2-(Din+2.*lc).^2).*kp)./(1-M2^2);
    a7 = -2.*M2.*(1i.*k+2.*Din./((Dv.^2-(Din+2.*lc).^2).*kp))./(1-M2^2);
    a8 = (k.^2-4.*1i.*k.*Din./((Dv.^2-(Din+2.*lc).^2).*kp))./(1-M2^2);
end

function [r1,r2,r3,r4] = calcCoeffR(a1,a2,a3,a4,a5,a6,a7,a8)
    r1 = ((a1+a7)+sqrt((a1-a7)^2+4*a3*a5))./2;
    r3 = ((a1+a7)-sqrt((a1-a7)^2+4*a3*a5))./2;
    r2 = ((a2+a8)-sqrt((a2-a8)^2+4*a4*a6))./2;
    r4 = ((a2+a8)+sqrt((a2-a8)^2+4*a4*a6))./2;
end

function beita = calcCoeffBeita(r1,r2,r3,r4)
    beita(1) = -(r1+sqrt(r1^2-4*r2))./2;
    beita(2) = -(r1-sqrt(r1^2-4*r2))./2;
    beita(3) = -(r3-sqrt(r3^2-4*r4))./2;
    beita(4) = -(r3+sqrt(r3^2-4*r4))./2;
end

function Fai = calcCoeffFai(beita,a1,a2,a3,a4)
    m = 1:4;
    Fai(1,m) = 1;
    Fai(2,m) = -(beita(m).^2 + a1*beita(m) + a2)./(a3*beita(m) + a4);
    Fai(3,m) = 1./beita(m);
    Fai(4,m) = Fai(2,m)./beita(m);
end

function AperB = calcCoeffAperB(k,areaRatio,M2,lb2)
    AperB = k.*(cot(k.*lb2) - ((areaRatio+1).*(1-1./areaRatio)^2.*(1-exp(1i.*k.*lb2.*M2)) ...
             + 1i.*(1+areaRatio.*(1-1./areaRatio)^2).*M2.*sin(k.*lb2).*cos(k.*lb2))...
             ./(sin(k.*lb2).*((1+areaRatio.*(1-1./areaRatio)^2).*cos(k.*lb2) ...
             + 1i.*M2.*sin(k.*lb2))));
end

function CperD = calcCoeffCperD(k,areaRatio,M2,lb2)
    CperD = k.*(cot(k.*lb2) + ((1-(1-1/areaRatio)^2).*(exp(-1i.*k.*lb2.*M2)-1))...
             ./(sin(k.*lb2).*((1+areaRatio.*(1-1/areaRatio)^2).*cos(k.*lb2) + 1i.*M2.*sin(k.*lb2))));
end

function Fi = calcCoeffFi(k,lb1,lp2,beita,Fai,AperB,CperD)
    Fi = ((Fai(2,3)+k.*tan(k.*lb1).*Fai(4,3)).*(1+AperB.*Fai(3,4)-Fai(2,4)-CperD.*Fai(4,4)).*exp(beita(4).*lp2)...
      -(Fai(2,4)+k.*tan(k.*lb1).*Fai(4,4)).*(1+AperB.*Fai(3,3)-Fai(2,3)-CperD.*Fai(4,3)).*exp(beita(3).*lp2));
end

function [C31,C32,C41,C42] = calcCoeffC(k,lb1,lp2,AperB,CperD,Fai,Fi,beita)
    C31 = ((Fai(2,4)+k.*tan(k.*lb1).*Fai(4,4)).*(1+AperB.*Fai(3,1)-Fai(2,1)-CperD.*Fai(4,1)).*exp(beita(1).*lp2)...
           -(Fai(2,1)+k.*tan(k.*lb1).*Fai(4,1)).*(1+AperB.*Fai(3,4)-Fai(2,4)-CperD.*Fai(4,4)).*exp(beita(4).*lp2))./Fi;
    C32 = ((Fai(2,4)+k.*tan(k.*lb1).*Fai(4,4)).*(1+AperB.*Fai(3,2)-Fai(2,2)-CperD.*Fai(4,2)).*exp(beita(2).*lp2)...
           -(Fai(2,2)+k.*tan(k.*lb1).*Fai(4,2)).*(1+AperB.*Fai(3,4)-Fai(2,4)-CperD.*Fai(4,4)).*exp(beita(4).*lp2))./Fi;
    C41 = ((Fai(2,1)+k.*tan(k.*lb1).*Fai(4,1)).*(1+AperB.*Fai(3,3)-Fai(2,3)-CperD.*Fai(4,3)).*exp(beita(3).*lp2)...
           -(Fai(2,3)+k.*tan(k.*lb1).*Fai(4,3)).*(1+AperB.*Fai(3,1)-Fai(2,1)-CperD.*Fai(4,1)).*exp(beita(1).*lp2))./Fi;
    C42 = ((Fai(2,2)+k.*tan(k.*lb1).*Fai(4,2)).*(1+AperB.*Fai(3,3)-Fai(2,3)-CperD.*Fai(4,3)).*exp(beita(3).*lp2)...
           -(Fai(2,3)+k.*tan(k.*lb1).*Fai(4,3)).*(1+AperB.*Fai(3,2)-Fai(2,2)-CperD.*Fai(4,2)).*exp(beita(2).*lp2))./Fi;
end

function sig0 = calcSigma0(bp)
    sig0 = 0.8216*(1-1.5443*bp^0.5+0.3508*bp+0.1935*bp^1.5);%穿孔管的端部修正系数（穿孔率小于40%）
if bp >= 0.4
	error('穿孔率不能大于40%');
end
end