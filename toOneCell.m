function oneCell = toOneCell( headerData,varargin )
%把多个矩阵合并为一个cell，用于保存为xls
%   headerData 表头，如果不要，设置为{}
%   varargin 数据矩阵
oneCell = {};
while length(varargin)>0
    dataMat = varargin{1};
    varargin = varargin(2:end);
    numCell = num2cell(dataMat);
    oneCell = cellPush2Right(oneCell,numCell);
end
if length(headerData)>0
    if size(headerData,1)>1
        headerData = headerData';
    end
	oneCell = cellPush2Bottom(headerData,oneCell);
end

end

