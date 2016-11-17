function M = innerReduceTransferMatrix(Dr1,Dr2,Lr,varargin)
%逐渐缩颈管传递矩阵
%                            
%                   |Lr     
%                    \           
%            Dr1  Dr2|\    
%                     /
%                    /      
%                   |Lc|    cone锥形   
%                    
pp = varargin;
k = nan;
a = nan;%声速
sigma = 0;%修正系数
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
        case 'sigma'
            sigma = val;%修正系数
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

r1 = Dr1./2;
S1 = pi.*Dr1.^2./4;
r2 = Dr2./2;
S2 = pi.*Dr2.^2./4;
Lc = (r1./(r1-r2)).*Lr;
x1 = Lc.*r1./(r2-r1);
A = (r2./r1).*cos(k.*Lc) - (1./(k.*x1)).*sin(k.*Lc);
B = 1i.*(a./S1).*(r1./r2).*sin(k.*Lc);
C = 1i.*(S1./a).*((r2./r1+1./(k.^2.*x1.^2)).*sin(k.*Lc) ...
    - (Lc./(k.*x1.^2)).*cos(k.*Lc));
D = (r1./r2).*(cos(k.*Lc) + (1./(k.*x1)).*sin(k.*Lc));
M = [D./(A.*D-B.*C),-B./(A.*D-B.*C);
    C./(B.*C-D.*A),-A./(B.*C-D.*A)];
end

