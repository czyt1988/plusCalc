clc;
clear all;
syms AL3 BL3 CL3 DL3%»¯¼ò
syms RV2 SV2 TV2 UV2
syms AL2 BL2 CL2 DL2%»¯¼ò
syms RV1 SV1 TV1 UV1
syms AL1 BL1 CL1 DL1%»¯¼ò

M = [AL3,BL3;CL3,DL3] * [RV2,SV2;TV2,UV2] * [AL2,BL2;CL2,DL2] * [RV1,SV1;TV1,UV1] * [AL1,BL1;CL1,DL1];
coefficient = -M(1,2)/M(1,1);
pretty(coefficient)