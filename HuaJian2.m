%% 验证闭口孔管传递矩阵化简计算
syms T11 T12 T21 T22 P1 B P2 D m1 m2 c s1 s2
[ansP2,ansM2] = solve([P1;m1*c/s1]==[T11,T12;T21,T22]*[P2;m2*c/s2],P2,m2);
ansP2 = simple(ansP2);
ansM2 = simple(ansM2);
disp('P2 = ');
pretty(ansP2)
disp('m2 = ');
pretty(ansM2)
[ansP1,ansM1] = solve([P2;m2*c/s2]==[T11,T12;T21,T22]*[P1;m1*c/s1],P1,m1);
ansP1 = simple(ansP1);
ansM1 = simple(ansM1);
disp('P1 = ');
pretty(ansP1)
disp('m1 = ');
pretty(ansM1)
%%
% 
% syms P1(0) p0 c u1(0) P2(0) u2(0) P1(lp2) u1(lp2) P2(lp2) u2(lp2) lb1 lp2 k i
% 
% P2(0) = (p0.*c.*u2(0))./(-1.*i.*tan(k.*lb1));
% P2(lp2) = (p0.*c.*u2(lp2))./(-1.*i.*tan(k.*(lb1+lp2)));
% 
% eq1 = R11.*P1(lp2) + R12.*p0.*c.*u1(lp2) + R13.*P2(lp2) + R14.*p0.*c.*u2(lp2) - P1(0);
% eq2 = R21.*P1(lp2) + R22.*p0.*c.*u1(lp2) + R23.*P2(lp2) + R24.*p0.*c.*u2(lp2) - p0.*c.*u1(0);
% eq3 = R31.*P1(lp2) + R32.*p0.*c.*u1(lp2) + R33.*P2(lp2) + R34.*p0.*c.*u2(lp2) - P2(0);
% eq4 = R41.*P1(lp2) + R42.*p0.*c.*u1(lp2) + R43.*P2(lp2) + R44.*p0.*c.*u2(lp2) - p0.*c.*u2(0);
% 
% D = solve(eq1,eq2,eq3,eq4,P1(0),u2(0));
% 
% P1(0) = D.P1(0)
% u2(0) = D.u2(0)
