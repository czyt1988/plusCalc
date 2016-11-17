function [ multFreMag,multFrePh ] = calcBaseFres( rawDatas,Fs,baseFrequency,multFreTimes,allowDeviation )
%计算一串列数据的基本频率的幅值
[rawFre,rawMag,rawPh] = fft_byColumn(rawDatas,Fs);%原始数据频谱分析
multFreMag = [];
multFrePh = [];
for i=1:multFreTimes%计算倍频幅值，相位
    [m,p] = fun_findBaseFres(rawFre,rawMag,rawPh,baseFrequency*i,allowDeviation);
    multFreMag(i,:) = m;
    multFrePh(i,:) = p;
end

end

