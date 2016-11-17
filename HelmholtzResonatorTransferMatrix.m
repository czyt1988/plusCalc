function M = HelmholtzResonatorTransferMatrix(V,lv,lc,dp,n,varargin)
%亥姆霍兹共鸣器传递矩阵
% lv 共鸣器长
% lc 共鸣器连接管长
% Dp 共鸣器连接管直径 dp*n                         
%       __________                
%      |          |                   
%      |    V     | lv
%      |___    ___|     
%          |  | lc        
% _________|Dp|__________                  
% _______________________                   
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

Sv = V./lv;
Dp = dp.*n;
Sc = pi.*Dp.^2./4;
% sigma = 0.8*Dp;
lcp = lc + sigma;
X = -1i.*(Sc./a).*((tan(k.*lcp)+(Sv./Sc).*tan(k.*lv))./(1-(Sv./Sc).*tan(k.*lcp).*tan(k.*lv)));
M = [1,0;
    X,1];
end

