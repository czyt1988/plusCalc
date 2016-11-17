function fhCalc()
clc;
clear;
close all;
currentPath = fileparts(mfilename('fullpath'));
%  长度 L1         Lv1       Lc       Lv2         L3
%              __________         __________
%             |          |       |          |
%  -----------|          |-------|          |-------------
%             |__________|       |__________|  
%  直径 Dpipe       Dv1      Dc       Dv2          Dpipe
%面积                        Ac
%体积                V1                 V2
	Lv1 = 0.1:0.2:0.9;
	Lc = 0:0.01:1.2;
	for ii = 1:length(Lv1)
		sd(ii).a = 345;
		sd(ii).Lc = Lc;%连接管长
		sd(ii).Dc = 0.157;
		sd(ii).Lv1 = Lv1(ii);
		sd(ii).Lv2 = 1-Lv1(ii);
		sd(ii).Dv1 = 0.5;
		sd(ii).Dv2 = 0.5;
		sd(ii).Ac = pi*sd(ii).Dc^2./4;
		sd(ii).V1 = pi*sd(ii).Dv1^2/4*sd(ii).Lv1;
		sd(ii).V2 = pi*sd(ii).Dv2^2/4*sd(ii).Lv2;
	end
	X = Lc;
	figure
	hold on;
	H = arrayfun(@(s)(plotData(s,X)),sd);
    figure
	hold on;
    H = arrayfun(@(s)(plotDiff(s,X,1)),sd);
end

function h = plotData(sd,X)
	fh = filterfrequency_vcv(sd.a,sd.Ac,sd.Lc,sd.Dc,sd.V1,sd.V2);
	h = plot(X,fh);
    text(X(end),fh(end),sprintf('Lv1:%g,Lv2:%g',sd.Lv1,sd.Lv2));
end

function h = plotDiff(sd,X,N)
	fh = filterfrequency_vcv(sd.a,sd.Ac,sd.Lc,sd.Dc,sd.V1,sd.V2);
    y = diff(fh,N);
	h = plot(X(N+1:end),y);
    text(X(end),y(end),sprintf('Lv1:%g,Lv2:%g',sd.Lv1,sd.Lv2));
end