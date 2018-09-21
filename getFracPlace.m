function [ res ] = getFracPlace( r )
    res = find(mod(r*power(10,0:20),1)==0);
    res = res(1)-1;  
end

