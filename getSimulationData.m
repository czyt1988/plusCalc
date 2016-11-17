function [ fre,mag ] = getSimulationData( matPath,structSection )
%获取模拟的数据
%  structSection  结构体的字段，如果没有设置，默认为rawData
if nargin == 1
    structSection = 'rawData';
end
load(matPath);
dsData = getfield(dataStruct,structSection);
mag = dsData.Mag;
fre = dsData.Fre;
end

