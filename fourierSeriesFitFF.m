function [A0,An,Ph ] = fourierSeriesFitFF( X,Y,N,Fre )
%傅里叶级数拟合
%   根据给定的x，y和设定的频率Fre，用N阶傅里叶级数进行拟合
%  最后拟合成 A0 + ΣAn cos(n Fre t + Ph)
if size(X,1) == 1
    error('X需要列向量');
end
hFun = fsFuncHandleMaker(N,Fre);
ft = fittype(hFun);
fitResult = fit(X,Y,ft);
A0 = fitResult.A0;
for ii=1:N
    An(ii) = eval(sprintf('fitResult.An%d',ii));
    Ph(ii) = eval(sprintf('fitResult.Ph%d',ii));
end
end

% function fp = fsFuncHandleMaker(N,Fre)
%     
% %     str = sprintf(...
% %     ['@(A0,An,Ph,X) '...
% %         ,'(Y = A0 .* ones(length(X),1);\n'...
% %         ,'for ii = 1:%d\n'...
% %         ,'   Y = Y + (An(ii) .* cos( ii*%g .* X + Ph(ii).*ones(length(X),1)));\n'...
% %         ,'end)'...
% %     ],N,Fre);
% 
%     strHead = '@(A0,An,Ph,x)';
%     str = '(A0 .* ones(length(x),1)';
%     for ii = 1:N
%         str = [str,'+',...
%             sprintf('(An(%d) .* cos( (%d*%g) .* x + Ph(%d).*ones(length(x),1)))',ii,ii,Fre,ii)];
%     end
%     str = [strHead,str,')'];
%     disp(str);
%     fp = str2func(str);
% end

function fp = fsFuncHandleMaker(N,Fre)   
    strHead = '@(A0'
    for ii = 1:N
        strHead = [strHead,','];
        strHead = [strHead,sprintf('An%d',ii)];
    end
    for ii = 1:N
        strHead = [strHead,','];
        strHead = [strHead,sprintf('Ph%d',ii)];
    end
    strHead = [strHead,',x)'];
    str = '(A0 .* ones(length(x),1)';
    for ii = 1:N
        str = [str,'+',...
            sprintf('(An%d .* cos( (%d*%g) .* x + Ph%d.*ones(length(x),1)))',ii,ii,Fre,ii)];
    end
    str = [strHead,str,')'];
    fp = str2func(str);
end
