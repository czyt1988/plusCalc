clc;
close all;
clear;
currentPath = fileparts(mfilename('fullpath'));
Y = [5.43093333333333,5.77926666666667,5.88726666666667,5.46006666666667,5.34813333333333,4.82420000000000,3.27413333333333,2.17446666666667,2.86106666666667,4.16066666666667,5.20326666666667,5.46626666666667,5.77926666666667,5.42706666666667,2.09033333333333];
X = [5,7,8,9,10.5,12,14,16,18,20,22,24,25,26,32];
ft = fittype( ['fourier',num2str(3)] );%3阶傅里叶拟合
fopts = fitoptions( ft );
fopts.Display = 'Off';
[fitresult, gof] = fit( X',Y', ft, fopts );
save(fullfile(currentPath,'fitresultCfit.mat'),'fitresult');
X = 5:32;
yFit = feval(fitresult,X);
fh = fopen(fullfile(currentPath,'fitresultCfit.txt'),'w');
fprintf(fh,'长度，值\r\n');
for i=1:length(X)
    fprintf(fh,'%g\t%g\r\n',X(i),yFit(i));
end
fprintf(fh,'-------------------------------------------\r\n');
fprintf(fh,'Fourier-order:%d',3);
fprintf(fh,'形如:fitresult(x) =  a0 + a1*cos(x*w) + b1*sin(x*w) + a2*cos(2*x*w) + b2*sin(2*x*w) + a3*cos(3*x*w) + b3*sin(3*x*w)');
fprintf(fh,'a0:%g\r\n',fitresult.a0);
fprintf(fh,'a1:%g\r\n',fitresult.a1);
fprintf(fh,'a2:%g\r\n',fitresult.a2);
fprintf(fh,'a3:%g\r\n',fitresult.a3);
fprintf(fh,'b1:%g\r\n',fitresult.b1);
fprintf(fh,'b2:%g\r\n',fitresult.b2);
fprintf(fh,'b3:%g\r\n',fitresult.b3);
fprintf(fh,'w:%g\r\n',fitresult.w);
fclose(fh);