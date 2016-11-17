clc;
clear all;
syms P1 M1 P2 M2 A
S = [1,0;-A,1]*[P1;M1];
pretty(S)
P1 = P2;
M1 = M2 + A*P2;
P2 = P1
M2 = M1- A*P1

% 
% dmM = simple(dmM)
% pretty(dmM)