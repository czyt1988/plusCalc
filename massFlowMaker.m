function [massFlow,time,meanMassFlow,meanFlowVelocity] = massFlowMaker(DCylinder,dPipe,rpm...
	,crank,connectingRod,densityOutlet,varargin)
%质量流量生成
%   alpha 相对余隙容积
%   pressureRadio 压力比（排气压力/吸气压力）
alpha=0.25;%相对余隙容积
pressureRadio=2;%压力比（排气压力/吸气压力）
k=1.4;%绝热指数
% DCylinder=250/1000;%缸径m
% dPipe=98/1000;%管内径m
% rpm=300;%转速r/min
% crank=140/1000;%曲柄长
% connectingRod=1075;%连杆长度
fs = 200;
oneSecond = 0;
while length(varargin)>=2
    prop =varargin{1};
    val=varargin{2};
    varargin=varargin(3:end);
    switch lower(prop)
        case 'relativeclearancevolume' %
            alpha = val;
        case 'rcv' %
            alpha = val;
        case 'pressureradio'
            pressureRadio = val;
        case 'pr'
            pressureRadio = val;
        case 'k'
            k = val;
        case 'fs'
        	fs = val;
        case 'onesecond'
            oneSecond = val;
        otherwise
        	error('参数错误:%s',prop);
    end
end
rps = rpm / 60;
spr = 1 / rps;%一转需要的秒
secNumPerRound = spr*fs;
sectionRang=linspace(0,360,secNumPerRound);
beita=asind(1-2.*(1+alpha).*(pressureRadio.^(-(1/k)))+2.*alpha);%角度不是弧度
alphac=270+beita;%盖侧缸排气阀的开启角
alphacx=90+beita;%轴侧缸排气阀的开启角
%双作用气缸排气速度曲线（盖侧进排气）
ACylinder=(pi.*(DCylinder.^2))/4;%气缸通流面积
APipe=(pi.*(dPipe.^2))/4;%管道通流面积
b=ACylinder/APipe;
lameda=crank/connectingRod;%曲柄长/连杆长
oumiga=(pi.*rpm)/30;%曲柄的角速度rad/s

%计算流速
massFlow=zeros(size(sectionRang));
for i=1:length(sectionRang)%盖侧单作用
    alpha=sectionRang(i);
    if (alpha<alphacx&&alpha>=0)
        v11=0;
    elseif (alpha<=180&&alpha>=alphacx)
        v11=b.*crank.*oumiga.*abs((sind(alpha)+(lameda/2).*sind(2.*alpha)));
    elseif (alpha<=alphac&&alpha>180)
        v11=0;
    elseif (alpha<=360&&alpha>=alphac)  
        v11=b.*crank.*oumiga.*abs((sind(alpha)+(lameda/2).*sind(2.*alpha)));
    end
    massFlow(i)=v11;
end

massFlow = massFlow.*densityOutlet.*APipe;
time = 0:1/fs:((length(massFlow)-1)/fs);
meanMassFlow = trapz(time,massFlow);
meanMassFlow = meanMassFlow/time(end);
meanFlowVelocity = (meanMassFlow/densityOutlet/APipe);
if oneSecond
    while length(massFlow) < fs+1
        massFlow = [massFlow,massFlow];
    end
    massFlow = massFlow(1:fs+1);
    time = 0:1/fs:1;
end

% for i=1:length(sectionRang)%盖侧单作用
%     alpha=sectionRang(i);
%     if (alpha<=alphacx && alpha>=0)
%         v11=0;
%     elseif (alpha<=180 && alpha>alphacx)
%         v11=densityOutlet .* ACylinder .* crank .* oumiga .* abs((sind(alpha)+(lameda/2).*sind(2.*alpha)));
%     elseif (alpha<=alphac && alpha>180)
%         v11=0;
%     elseif (alpha<=360&&alpha>alphac)  
%         v11=densityOutlet .* ACylinder .* crank .* oumiga.*abs((sind(alpha)+(lameda/2).*sind(2.*alpha)));
%     end
%     massFlow(i)=v11;
% end

end