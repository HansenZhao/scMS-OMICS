function [ SNcell ] = getSCSN(I,maxStep)
    if ~exist('maxStep','var')
        maxStep = 3;
    end
    L = length(I);
    SNcell = cell(L,1);
    cellCounter = 0;
    indexCounter = 0;
    curIndex = zeros(maxStep,1);
    for m = 1:L
        indexCounter = indexCounter + 1;
        curIndex(indexCounter) = I(m);
        if m < L
            diff = I(m+1) - I(m);
            if diff > 1 || indexCounter == maxStep
                cellCounter = cellCounter + 1;
                SNcell{cellCounter} = curIndex(1:indexCounter);
                curIndex = zeros(maxStep,1);
                indexCounter = 0;
            end               
        else
            cellCounter = cellCounter + 1;
            SNcell{cellCounter} = curIndex(1:indexCounter);
        end
    end
    SNcell((cellCounter+1):end) = [];
end

