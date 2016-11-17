function [out_index,meadUpStd,meadDownStd,meanValue,stdValue] = fun_sigma_outlier_detection( series , stdTimes )
%西格玛奇异值检测，将返回stdTimes以外的索引
% series 为输入的向量
% stdTimes为标准差倍数，stdTimes为3时，置信度达到99%
% out_index 超出这个sigma范围的数的索引
% meadUpStd ：meanValue+stdTimes*stdValue
% meadDownStd ：meanValue-stdTimes*stdValue;
% meanValue : 均值
% stdValue 标准差

stdValue = std(series);
meanValue = mean(series);
meadUpStd = meanValue+stdTimes*stdValue;
meadDownStd = meanValue-stdTimes*stdValue;
out_index = find((series > (meadUpStd)) | (series < (meadDownStd)) );

end

