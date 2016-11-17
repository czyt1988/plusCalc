clc;
clear all;
syms fL3 aps gL3%f(L3) = fL3.g(L3) = gL3 aps = -a/S
M_L3 = [fL3,aps*gL3;...
    1/aps*gL3,fL3];
syms Alv2_ Blv2_ Clv2_ Dlv2_;
M_V2 = [Alv2_,aps*Blv2_;...
    (1/aps)*Blv2_,Alv2_];
syms fL2 aps gL2
M_L2 = [fL2,aps*gL2;...
    1/aps*gL2,fL2];
syms Alv1l1_ Blv1l1_ Clv1l1_ Dlv1l1_;
M_V1L1 = [Alv1l1_,Blv1l1_;...
    Clv1l1_,Dlv1l1_];

M_total = M_L3 * M_V2 * M_L2 * M_V1L1; 

M_total = simple(M_total);
NpM = M_total(1,2)/M_total(1,1)
NpM = simple(NpM);
pretty(NpM)