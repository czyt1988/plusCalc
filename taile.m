%傅里叶级数拟合
N = 8;
Fre = 10;
[A0,An,Ph] = fourierSeriesFitFF(time,massFlowRaw,N,Fre);
figure
hold on;
plot(time,massFlowRaw,'-r');
dN = length(massFlowRaw);
wave = ones(dN,1).*A0;
for ii = 1:N
    wave = wave + An(ii).*cos((ii * Fre) .* time + Ph(ii).*ones(dN,1));
end
plot(time,wave,'-b');