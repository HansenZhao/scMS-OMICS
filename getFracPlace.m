function [ res ] = getFracPlace( r )
    res = find(abs(r*power(10,0:20) - round(r*power(10,0:20)))<1e-9);
    res = res(1)-1;  
end

