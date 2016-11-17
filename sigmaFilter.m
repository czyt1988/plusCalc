function filterData = sigmaFilter(dataRaw,sigma)
%西格玛奇异值滤波
% dataRaw原始数据，sigma的倍数
if nargin<2
	sigma = 3;
end
	oi = fun_sigma_outlier_detection(dataRaw,sigma);
	dataRaw(oi) = [];
	filterData = dataRaw;
end