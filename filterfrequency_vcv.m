function fH = filterfrequency_vcv( a,Ac,Lc,Dc,V1,V2 )
%计算截止频率
%  长度 L1         Lv1       Lc       Lv2         L3
%              __________         __________
%             |          |       |          |
%  -----------|          |-------|          |-------------
%             |__________|       |__________|  
%  直径 Dpipe       Dv1      Dc       Dv2          Dpipe
%面积                        Ac
%体积                V1                 V2
% a 声速
% Ac 中间插管的截面积
Lc1 = Lc + 0.6*Dc;
fH = a./(2*pi)*((Ac./Lc1)*((1./V1)+(1./V2))).^0.5;
end

