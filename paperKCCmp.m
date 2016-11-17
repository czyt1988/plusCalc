
clc;
clear;
close all;
%%
mc = 0.5^2/0.157^2;
Lv = 0.926;
LHalf = 0:0.1:2;

KC2V = sameDoubleVesselKC(mc,LHalf,Lv,'a',345,'f',10);

figure
h2 = plot(LHalf,abs(KC2V),':b','LineWidth',1.5);
grid on;

xlabel('length of L2(m)');
ylabel('KC');

set(gcf,'color','w');