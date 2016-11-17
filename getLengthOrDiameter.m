function res  = getLengthOrDiameter( V,data,LOrD )
%根据体积获得长或直径
%   LOrD 为1时，输入的是长度。
%   LOrD 为2时输入的是直径
L = nan;
D = nan;
if 1 == LOrD
	L = data;
else
	D = data;
end
if isnan(L)
	res = V./(pi.*D.^2./4);
else
	res = (4.*V./(L.*pi)).^0.5;
end
end

