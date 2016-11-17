function d = arrayPresentSection( data,sec )
%提取一段array的百分比的内容
%   sec是一个长度为2的数据，0代表前段从第一个点开始算，1代表算到最后一个点，如果第一个为0.1则从总体的10%开始计算
lowIndex = floor(sec(1)*length(data))+1;
upIndex = floor(sec(2)*length(data));
d = data(lowIndex:upIndex);
end

