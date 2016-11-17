%% 双容与单容对比
clc;
close all;
clear;
currentPath = fileparts(mfilename('fullpath'));
%  长度 L
%  -----------------------         
%  直径 Dpipe       
%% 质量流量获取
massFlow = load(fullfile(currentPath,'mass_flow_0.1478_NorthZone.txt'));
time = massFlow(:,1);
massFlow = massFlow(:,2);
massFlowE1 = fft(detrend(massFlow));
%% 基本参数设置
isXLength = 1;%如果X轴是距离将会按照实际尺寸来绘制。
opt.frequency = 10;
opt.acousticVelocity = 345;
Dpipe = 0.157;%管道直径（m）

L = [25.5:0.1:26.5];%直管长度
for i=1:length(L)
    sectionL{i} = 1:L(i);
end
%% 计算不同长度直管的脉动压力
for i=1:length(L)    
    [plusStraight{i},pressureStraight{i}] = fun_straightPipe(massFlowE1,L(i),Dpipe,opt,sectionL{i});
end

%% 绘图
marker_style = {'-o','-d','-<','-s','->','-<','-p','-*','-v','-^','-+','-x','-h'};
color_style = [36,100,196;...
    18,175,134;...
    237,144,10;...
    131,54,229;...
    245,18,103;...
    255,99,56;...
    ]./255;
market_style_length = length(marker_style);
color_style_length = size(color_style,1);
marker_index = 0;
color_index = 0;

figure
hold on;
for i = 1:length(plusStraight)
    marker_index = marker_index + 1;
    if marker_index > market_style_length
        marker_index = 1;
    end
    color_index = color_index + 1;
    if color_index > color_style_length
        color_index = 1;
    end
    Y = plusStraight{i};
    if isXLength
        X = sectionL{i};
    else
        X = linspace(0,max(L),length(Y));
    end

    h(i) = plot(X,Y./1000,marker_style{marker_index},'color',color_style(color_index,:),'LineWidth',1.5);
    textLegend{i}=sprintf('%g m直管',L(i));
end
legend(h,textLegend,0);
grid on;
if isXLength
    xlabel('距离(m)');
else
    xlabel('测点');
end
ylabel('脉动压力峰峰值(kPa)');
set(gcf,'color','w');

