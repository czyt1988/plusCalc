function [plusData,filterData] = calcPuls( pressure,dcpss )
%计算气流脉动
%   pressure 输入的原始压力
%   其他设置，用getDefaultCalcPulsSetStruct可获取默认的设置结构体
if nargin < 2 
	dcpss = getDefaultCalcPulsSetStruct();
end
%确定数据个数
if (1 == size(pressure,2) || 1 == size(pressure,1))
	DATA_COUNT = 1;
end
if dcpss.dim
	DATA_COUNT = size(pressure,2);
else
	DATA_COUNT = size(pressure,1);
end
%先把数据转移到cell里
for ii = 1:DATA_COUNT
	if 1 == dcpss.dim
		filterData{ii} = pressure(ii,:);
	else
		filterData{ii} = pressure(:,ii);
	end
end
%是否进行sigma滤波
if ~isnan(dcpss.sigma)
	filterData = cellfun(@(xx) sigmaFilter(xx,dcpss.sigma),filterData,'UniformOutput',0);
end
%进行高通滤波
if dcpss.isHp
	if isnan(dcpss.fs)
		error('设置了高通滤波，需要设置采用频率，fs');
	end
	filterData = cellfun(@(xx) highp(xx,dcpss.f_pass,dcpss.f_stop,dcpss.rp,dcpss.rs,dcpss.fs) ...
		,filterData,'UniformOutput',0);
end
%最后计算
filterData = cellfun(@(xx) arrayPresentSection(xx,dcpss.calcSection),filterData,'UniformOutput',0);
plusData = cellfun(@(xx) max(xx) - min(xx),filterData);

end




