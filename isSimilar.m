function [index,I] = isSimilar(vec,pool,thres)
    vec = vec(:);
    pool = pool(:)';
    diffMat = abs(vec-pool);
    [minDiff,index] = min(diffMat,[],2);
    I = minDiff <= thres;
    index(~I) = 0;
end