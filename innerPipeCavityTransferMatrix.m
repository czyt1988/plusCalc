function M = innerPipeCavityTransferMatrix(D,d,l,varargin)
%内插管传递矩阵
%   D 缓冲罐内径
%   d 内径
%   
%    
% ________
% __l____|
% __d_____  D
% _______|
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

S = pi.*D.^2./4;
Spipe = pi.*d.^2./4;
Sa = S - Spipe;
lp = l + sigma;
M = [1,0;
    -1i.*(Sa./a).*tan(k.*lp),1];
end

