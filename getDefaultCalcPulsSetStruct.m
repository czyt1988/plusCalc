function dcpss = getDefaultCalcPulsSetStruct()
	dcpss.dim = 2;%计算维度 1代表pressure每一列代表一组数据
	dcpss.sigma = 2.8;%sigma去噪的sigma的倍数，3代表置信度97%，如果为nan代表不进行sigma滤波
	dcpss.calcSection = [0,1];%计算区域，0代表前段从第一个点开始算，1代表算到最后一个点，如果第一个为0.1则从总体的10%开始计算
	dcpss.isHp = 0;%是否进行高通滤波
	dcpss.f_pass = 4;%通过频率5Hz
	dcpss.f_stop = 2;%截止频率3Hz
	dcpss.rp = 0.1;%边带区衰减DB数设置
	dcpss.rs = 30;%截止区衰减DB数设置
	dcpss.fs = nan;
end