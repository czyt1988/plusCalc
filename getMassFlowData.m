function [ time,massFlow,Fre,massFlowE ] = getMassFlowData(varargin )
%获取质量流量
%   根据情况获取质量流量的数据

pp=varargin;
N = 4096;
isFindPeaks = 1;
isSin = 0;
sinAmp = nan;
Fs = nan;
sinFre = nan;
while length(pp)>=2
    prop =pp{1};
    val=pp{2};
    pp=pp(3:end);
    switch lower(prop)
        case 'n' %数据数量
            N = val;
        case 'isfindpeaks'
            isFindPeaks = val;
        case 'sin'
            sinAmp = val;%生成一个单一sin函数
        case'sinfre'
            sinFre = val;
        case 'fs'
            Fs = val;
    end
end
if isSin
    time = 0:1/Fs:((1/Fs)*(N-1))
    massFlow = sinAmp.*sin(2*pi*sinFre*time)
else
    currentPath = fileparts(mfilename('fullpath'));
    massFlow = load(fullfile(currentPath,'mass_flow_0.1478_NorthZone.txt'));
    N = 4096;
    time = massFlow(1:N,1);
    massFlowRaw = massFlow(1:N,2);
    Fs = 1/(time(2)-time(1));
    [FreRaw,AmpRaw,~,massFlowE] = fun_fft(detrend(massFlowRaw),Fs);
    % 提取主要频率
    if isFindPeaks
        [~,locs] = findpeaks(AmpRaw);
        Fre = FreRaw(locs);
        massFlowE = massFlowE(locs);
    end
    temp = Fre<80;%Fre<20 | (Fre>22&Fre<80);
    Fre = Fre(temp);
    massFlowE = massFlowE(temp);
end
end

