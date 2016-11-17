clc;
close all;
clear;
currentPath = fileparts(mfilename('fullpath'));
mf = load([currentPath,'\mass_flow_0.1478_NorthZone.txt']);
figure
subplot(1,2,1)
plot(mf(:,1),mf(:,2));
set(gcf,'color','w');
hlx = xlabel('time(s)');
hly = ylabel('mass flow rate(kg/s)');
xlim([0,1]);
ylim([0,0.8]);
set(gca,'FontUnits','points','FontSize',8,'FontName','Times New Roman');
set(hlx,'FontUnits','points','FontSize',12,'FontName','Times New Roman');
set(hly,'FontUnits','points','FontSize',12,'FontName','Times New Roman');
%set(gcf,'unit','pixels','position',[200,200,300,200]);
%set(gca,'LineWidth',1);

%set(gcf,'PaperUnits','centimeters');
savePath = fullfile(currentPath,'数值分析\论文用图\质量流量波形图.png');
% 其中 width， 9厘米为半栏图片宽度， 14厘米为整栏图片宽度。
% export_fig('质量流量波形图','width',9,'format','eps','Bounds','tight',...
%     'FontMode','scaled','DefaultFixedFontSize',7,'FontSizeMin',6,'FontSizeMax',7,...
%     'color','rgb','Preview','tiff');
subplot(1,2,2)
[Fre, Amp,Ph] = fun_fft( mf(:,2),1/(mf(2,1)-mf(1,1)));
plot(Fre,Amp);
set(gcf,'color','w');
hlx = xlabel('frequency(Hz)');
hly =ylabel('amplitude(kg/s)');
xlim([0,100]);
set(gca,'FontUnits','points','FontSize',8,'FontName','Times New Roman');
set(hlx,'FontUnits','points','FontSize',12,'FontName','Times New Roman');
set(hly,'FontUnits','points','FontSize',12,'FontName','Times New Roman');
set(gcf,'unit','centimeter','position',[1 1 9 4.5]);
savePath = fullfile(currentPath,'FIG4.EPS');
export_fig(savePath,'-eps','-r1000');%600dpi
