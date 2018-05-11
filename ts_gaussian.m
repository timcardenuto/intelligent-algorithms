function [out] = ts_gaussian(x,c1,sigma1,c2,sigma2)
    out = [];
    for i = 1:length(x)
        if(x(i) <= c1)
            out = [out; exp(-0.5 * (((x(i)-c1)/sigma1).^2))];
        elseif(x(i) < c2)
            out = [out; 1];
        else
            out = [out; exp(-0.5 * (((x(i)-c2)/sigma2).^2))];
        end
    end
end
