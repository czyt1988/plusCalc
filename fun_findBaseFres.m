function [mags,phs,fres] = fun_findBaseFres(rawFre,rawMag,rawPh,findFre,allowDeviation)
%UNTITLED2 Summary of this function goes here
    for i=1:size(rawMag,2)
        f = rawFre(:,i);
        mag = rawMag(:,i);
        p = rawPh(:,i);
        
        k = find(f>(findFre-allowDeviation) & f<(findFre+allowDeviation));
        [mags(1,i),index] = max(mag(k));%找到范围内最大的值
        fres(1,i) = f(k(index));
        phs(1,i) = p(k(index));
    end
end

