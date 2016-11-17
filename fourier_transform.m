%fourier_transform
order = '8';
ft = fittype( ['fourier',order]);%n阶傅里叶拟合
fopts = fitoptions( ft );
fopts.Display = 'Off';
[fitresult, gof] = fit( time,massFlowRaw, ft, fopts );
yFit = feval(fitresult,time);
figure
plot(time,massFlowRaw,'-r');
hold on;
plot(time,yFit,'-b');
a0 = fitresult.a0;
w = fitresult.w;
% General model Fourier8:
%      f(x) = 
%                a0 + a1*cos(x*w) + b1*sin(x*w) + 
%                a2*cos(2*x*w) + b2*sin(2*x*w) + a3*cos(3*x*w) + b3*sin(3*x*w) + 
%                a4*cos(4*x*w) + b4*sin(4*x*w) + a5*cos(5*x*w) + b5*sin(5*x*w) + 
%                a6*cos(6*x*w) + b6*sin(6*x*w) + a7*cos(7*x*w) + b7*sin(7*x*w) + 
%                a8*cos(8*x*w) + b8*sin(8*x*w)
% Coefficients (with 95% confidence bounds):
%        a0 =       0.148  (0.1455, 0.1504)
%        a1 =   0.0001121  (-0.003385, 0.003609)
%        b1 =  -0.0003344  (-0.003716, 0.003047)
%        a2 =     0.02126  (0.01784, 0.02468)
%        b2 =    -0.02074  (-0.02421, -0.01728)
%        a3 =  -0.0005929  (-0.004021, 0.002836)
%        b3 =  -0.0002264  (-0.003689, 0.003236)
%        a4 =    -0.05018  (-0.05585, -0.04452)
%        b4 =     -0.2323  (-0.2359, -0.2288)
%        a5 =  -0.0004436  (-0.003885, 0.002998)
%        b5 =   0.0004223  (-0.00303, 0.003875)
%        a6 =    -0.01227  (-0.01572, -0.008824)
%        b6 =    -0.01381  (-0.01727, -0.01034)
%        a7 =   0.0005987  (-0.002859, 0.004057)
%        b7 =  -4.584e-05  (-0.003472, 0.00338)
%        a8 =    -0.09223  (-0.09653, -0.08794)
%        b8 =     0.06834  (0.06332, 0.07336)
%        w =       15.71  (15.7, 15.72)
% 
% Goodness of fit:
%   SSE: 25.2
%   R-square: 0.8484
%   Adjusted R-square: 0.8477
%   RMSE: 0.07861
        %cos             sin
coeff = [fitresult.a1,fitresult.b1 ...
		;fitresult.a2,fitresult.b2 ...
		;fitresult.a3,fitresult.b3 ...
		;fitresult.a4,fitresult.b4 ...
		;fitresult.a5,fitresult.b5 ...
		;fitresult.a6,fitresult.b6 ...
		;fitresult.a7,fitresult.b7 ...
		;fitresult.a8,fitresult.b8 ...
    ];
for i=1:size(coeff,1)
    An(i) = coeff(i,1);%An对应cos
	Bn(i) = coeff(i,2);%Bn对应sin
	Phn(i) = atan(Bn(i)/An(i));
	Cn(i) = (An(i)^2 + Bn(i)^2)^0.5;% Cn = sqrt(an^2+bn^2)
end
Phn([3:6,8]) = Phn([3:6,8]) + pi;


Wave = zeros(length(time),1);%初始化
for i=1:length(Cn)
	Wave = Wave + ( Cn(i) .* cos( i*w .* time - Phn(i)*ones(length(time),1) ) );
end
Wave = a0.*ones(length(time),1) + Wave;

% Wave2 = zeros(length(time),1);%初始化
% for i=1:length(Cn)
% 	Wave2 = Wave2 + ( Cn(i) .* cos( i*w .* time + Phn(i)*ones(length(time),1) ) );
% end
% Wave2 = a0.*ones(length(time),1) + Wave2;
% Wave = zeros(length(time),1);%初始化
% for i=1:length(Cn)
% 	Wave = Wave + ( An(i) .* cos( i*w .* time )  + Bn(i) .* sin(i*w .* time) );
% end
% Wave = a0.*ones(length(time),1) + Wave;

figure
plot(time,massFlowRaw,'-r');
hold on;
plot(time,Wave,'-b');
% plot(time,Wave2,'-','color',[180,62,173]./255);

% figure
% i = 8;
% plot(( An(i) .* cos( i*w .* time )  + Bn(i) .* sin(i*w .* time) ) , '-r');
% hold on;
% plot(Cn(i) .* cos( i*w .* time - (Phn(i) )*ones(length(time),1) ) , '-b');
