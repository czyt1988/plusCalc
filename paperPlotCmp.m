function [ h1,h2 ] = paperPlotCmp( X1,Y1,X2,Y2,markCell,lineStyleCell,clr)
%论文绘图
%   绘制两个对比图,用同样的颜色淡不同的标记
h1 = plot(X1,Y1,'LineWidth',1.5,'color',clr,'Marker',markCell{1},'LineStyle',lineStyleCell{1});
hold on;
h2 = plot(X2,Y2,'LineWidth',1.5,'color',clr,'Marker',markCell{2},'LineStyle',lineStyleCell{2});
end

