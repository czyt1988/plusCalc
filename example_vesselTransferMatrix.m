%% 容器传递矩阵的运用及计算差距
%计算管容管容的脉动
%   massFlowE1 经过fft后的质量流量，直接对质量流量进行去直流fft
%  长度 L1     l    Lv1   l    L2  
%              __________        
%             |          |      
%  -----------|          |----------
%             |__________|       
% 直径 Dpipe       Dv1       Dpipe       
%   
%   opt.frequency 频率 Hz
%   opt.acousticVelocity 声速 366m/s
massFlow = load(fullfile(currentPath,'mass_flow_0.1478_NorthZone.txt'));
time = massFlow(:,1);
massFlow = massFlow(:,2);
massFlowE1 = fft(detrend(massFlow));
isDamping = 1;%是否计算阻尼
coeffFriction = 0.02;%管道摩察系数
meanFlowVelocity = 14.6;%管道平均流速
totalLength = 26;
L1 = 1;
Lv = 1;%缓冲罐长度
l = 0.115;%缓冲罐封头长度
Dpipe = 0.157;%管道直径（m）
Dv = 0.5;%缓冲罐的直径（m）
opt.frequency = 10;
opt.acousticVelocity = 366;
opt.isDamping = isDamping;
opt.coeffFriction = coeffFriction;
opt.meanFlowVelocity = meanFlowVelocity;
opt.isUseStaightPipe = 0;
L2 = totalLength-L1;
sectionL1 = 0:L1;
sectionL2 = 0:L2;
% 使用缓冲罐传递矩阵有阻尼形式计算的传递矩阵
textLegend{1} = '缓冲罐传递矩阵有阻尼形式';
[plusL1{1},plusL2{1}] = fun_oneTank(massFlowE1,L1,L2,Lv,l,Dpipe,Dv,opt,sectionL1,sectionL2);

textLegend{2} = '直管替代缓冲罐传递矩阵有阻尼形式';
opt.isUseStaightPipe = 1;
[plusL1{2},plusL2{2}] = fun_oneTank(massFlowE1,L1,L2,Lv,l,Dpipe,Dv,opt,sectionL1,sectionL2);

textLegend{3} = '缓冲罐传递矩阵无阻尼形式';
opt.isDamping = 0;
opt.isUseStaightPipe = 0;
[plusL1{3},plusL2{3}] = fun_oneTank(massFlowE1,L1,L2,Lv,l,Dpipe,Dv,opt,sectionL1,sectionL2);

textLegend{4} = '直管替代缓冲罐传递矩阵无阻尼形式';
opt.isDamping = 0;
opt.isUseStaightPipe = 1;
[plusL1{4},plusL2{4}] = fun_oneTank(massFlowE1,L1,L2,Lv,l,Dpipe,Dv,opt,sectionL1,sectionL2);

figure
h = [];
hold on;
X = 1:(length(plusL1{1}) + length(plusL2{1}));
h(1) = plot(X,[plusL1{1},plusL2{1}],'color',[64,98,237]./255,'LineWidth',2);
h(2) = plot(X,[plusL1{2},plusL2{2}],'color',[108,194,78]./255,'LineWidth',2);
h(3) = plot(X,[plusL1{3},plusL2{3}],'color',[147,52,234]./255);
h(4) = plot(X,[plusL1{4},plusL2{4}],'color',[49,201,113]./255);
grid on;
ylabel('maximum pulsating pressure(kPa)');
xlabel('mea point');
legend(h,textLegend);
title('example vessel transfer matrix');
set(gcf,'color','w');