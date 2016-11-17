function [ data ] = getAnalysisDataFun( type )
%获取分析的数据
    currentPath = fileparts(mfilename('fullpath'));
    switch(lower(type))
        case 'dvthrfre'%双容拼接0m理论
            [data.fre,data.mag]=getThrFreData();
        case 'dvsimfre'%双容拼接0m模拟    
            [data.fre,data.mag]=getSimulationData(fullfile(currentPath,'originData\26米双容间距0米多测点.mat'));
        case 'dvexpfre'%双容拼接0m实验
            [data.fre,data.mag]=getExpFreData();
        case 'strfre'%直管模拟
            [data.fre,data.mag]=getSimulationData(fullfile(currentPath,'originData\直管对应26米双容间距0.9米.mat'));
        case 'ovsimfre'%单容模拟
            [data.fre,data.mag]=getSimulationData(fullfile(currentPath,'originData\26米单容V=V1+V2.mat'));
    end
end

function [fre,mag] = getThrFreData()
    currentPath = fileparts(mfilename('fullpath'));
    load(fullfile(currentPath,'originData\理论-0m间距-f小于80.mat'));
    [fre,mag] = fft_byColumn(pppp,3600);
    mag = mag./1000;
end

function [fre,mag] = getExpFreData()
    currentPath = fileparts(mfilename('fullpath'));
    load(fullfile(currentPath,'originData\双容无间距开口4.mat'));
    fre = dataStruct.subSpectrumData.Fre;
    mag = dataStruct.subSpectrumData.Mag;
end
