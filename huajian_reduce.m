syms tAa taA ao ma s a st
tAa = s/st;
taA = st/s;
A = [1,(1-tAa^2)*ao*ma*a/s;
    0,1];
B = [1,-(2*ao*ma*a*(1-1/taA))/(s*taA);
    0,1];
C = B*A;
C = simple(C)
DD = (s-st)^2-(s*st^2);
DD = simple(DD)
pretty(C)

clear;
clc
s=0.001:0.01:5;
st = 0.001:0.01:5;
for i = 1:length(s)
    for j=1:length(st)
        if st(j) < s(i)
            DD(i,j) = nan;
        else
            DD(i,j) = (s(i)-st(j)).^2-(s(i).*st(j).^2);
        end
        if(abs(DD(i,j)) > 1)
            DD(i,j) = nan;
        end
    end
end
figure
[X,Y] = meshgrid(s,st);
contourf(X,Y,DD);
colorbar