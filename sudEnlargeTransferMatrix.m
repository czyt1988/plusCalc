function Matrix = sudEnlargeTransferMatrix( S_in,S_out,a,varargin)
%急速变径传递矩阵
%   S1进口截面积
%   S2出口截面积
%   S1BigThanS2 s1比s2大，说明是缩颈，否则为扩径
%   a声速
if(S_in > S_out)
    error('S_in 需要小于 S_out');
end

mach = 10/345;
coeffDamping = 0.03;
pp = varargin;
notMach = 1;
while length(pp)>=2
    prop =pp{1};
    val=pp{2};
    pp=pp(3:end);
    switch lower(prop)
        case 'coeffdamping' %阻尼系数
            coeffDamping = val;
        case 'damping' %阻尼系数
            coeffDamping = val;
        case 'mach' %马赫数，加入马赫数将会使用带马赫数的公式计算
            mach = val;
        case 'm'
            mach = val;
        case 'notmach'
            notMach = val;
        otherwise
       		error('参数错误%s',prop);
    end
end

M=0;
Tin_out = S_in/S_out;
if notMach
    mach = 0;
end
M=(1-Tin_out)*(2*coeffDamping*mach*a)/S_out;
if isnan(M)
    M = 0;
end

Matrix = [1,M;...
          0,1];

end

