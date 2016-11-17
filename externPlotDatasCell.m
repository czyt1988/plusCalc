function res = externPlotDatasCell(datasCell,varargin)
% []     ,'x','title1','title2','title3'
% 'data1',[x],[y]     ,[y]     ,[y]
% 'data2',[x],[y]     ,[y]     ,[y]
% 'data3',[x],[y]     ,[y]     ,[y]
%
% 'data1','data1' ,'data1' ,'data1' ,[],'data2','data2' ,'data2' ,'data2' ...
% 'x'    ,'title1','title2','title3',[],'x'    ,'title1','title2','title3',[]
%  1     , xx     ,  xx    , xx     ,[], 1     , xx     ,  xx    , xx     ,
%  2     , xx     ,  xx    , xx     ,[], 1     , xx     ,  xx    , xx     ,
%  3     , xx     ,  xx    , xx     ,[], 1     , xx     ,  xx    , xx     ,
%  4     , xx     ,  xx    , xx     ,[], 1     , xx     ,  xx    , xx     ,
%  ......
%  上述,'x','title1','title2','title3'叫paramLegend,可以自定义cell,或者定义paramLegendIndex
%  上述，'data1','data1' ,'data1' ,'data1'叫dataLegend,可以自定义cell
%  dataRowsIndex 为数据对应的索引，数据是以行为分布，默认除了第一行外的第二行到最后一行为dataRowsIndex
%                若有特殊情况不导出哪一行可以重新定义dataRowsIndex
%  paramIndex 为每个数据系列对应的变量索引，变量以列分布，默认除了第一列外的第二列到最后一列为paramIndex
%                若有特殊情况，不导出某一列变量，可以重新定义paramIndex
%  paramLegend 为param对应的描述，默认第一行从第二列到最后一列为paramLegend
%  dataLegend 为data对应的描述，默认第一列从第二行到最后一行为dataLegend
%

dataRowsIndexs = [2:size(datasCell,1)];
dataColumnIndex = [2:size(datasCell,2)];


dataParamLegend = {};
dataNameLegend = {};
while length(varargin)>=2
    prop =varargin{1};
    val=varargin{2};
    varargin=varargin(3:end);
    switch lower(prop)
    	case 'datarowsindexs'
    		dataRowsIndexs = val;
    	case 'datacolumnindex'
    		dataColumnIndex = val;
        case 'dataparamlegend' %是否允许补0
            dataParamLegend = val;
        case 'datanamelegend'
            dataNameLegend = val;
    end
end
if isempty(dataNameLegend)
    dataNameLegend = cell(length(dataRowsIndexs),1);
end
if isempty(dataParamLegend)
    dataParamLegend = cell(1,length(dataColumnIndex));
end
res = {};
externCells = {};
count = 1;
for i = dataRowsIndexs
    rowCell = datasCell(i,dataColumnIndex);
    externCells = exOneDataCells(rowCell,dataParamLegend,dataNameLegend{count});
    if count == 1
        res = cellPush2Right(res,externCells);
    else
        res = cellPush2Right(res,{''},externCells);
    end
    count = count + 1;
end
end


function exCell = exOneDataCells(rowDataCell,paramLegend,dataName)
colCount = size(rowDataCell,2);

exCell(2,:) = paramLegend;

for i=1:colCount
    exCell{1,i} = dataName;
	rowIndex = 3;
	data = rowDataCell{i};
	if ~(size(data,1) > 1)
		data = data';
    end
    
	exCell(rowIndex:(length(data)+rowIndex-1),i) = mat2cell(data,ones(length(data),1),[1]);
end

end