function d = pressureDrop_kPa(u,p,K)
%% 计算压力降
% u 前管道流速 m/s
% p 密度 kg/m3
% 阻力系数
    d = K .* u.^2.*p./2000;
end