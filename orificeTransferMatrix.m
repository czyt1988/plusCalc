function M = orificeTransferMatrix(D,d,Velocity )
%孔板传递矩阵
%   D 外径
%   d 内径
%   Velocity 平均速度
S = pi.*D.^2./4;
ks = (1+0.707./((1-(d.^2)./(D.^2)).^0.5)).^2 ...
.*...
((D./d).^2-1).^2;
A = -ks.*Velocity./S;
if isnan(A)
    A = 0;
end
M = [1,A;...
	0,1];
end

