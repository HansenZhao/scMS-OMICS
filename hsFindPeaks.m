function [ pLoc,pInt,v ] = hsFindPeaks( intens,halfWidth,thres,isShow)
    % [ pLoc,pInt,v ] = hsFindPeaks( intens,halfWidth = 10,thres = 3.5,isShow=0)
    if ~exist('halfWidth','var')
        halfWidth = 10;
    end
    if ~exist('thres','var')
        thres = 3.5;
    end
    if ~exist('isShow','var')
        isShow = 0;
    end
    if sum(intens==0)/length(intens) > 0.6
        thres = thres*0.4;
    end
    [pkInts,pkLoc] = findpeaks(intens);
    totL = length(intens);
    L = length(pkLoc);
    v = zeros(L,1);
    nPeak = sum(v>thres);
    inc = 1;
    newIntens = intens;
    tmpMean = zeros(L,1);
%     plot(intens);hold on;
    while(inc>0)
        for m = 1:L
            index = max([1,pkLoc(m)-halfWidth]):1:min([totL,pkLoc(m)+halfWidth]);
            index(index==pkLoc(m)) = [];
            v(m) = (pkInts(m) - mean(newIntens(index)))/(eps+std(newIntens(index)));
            tmpMean(m) = mean(newIntens(index));
        end
        newIntens(pkLoc(v>thres)) = tmpMean(v>thres);
        inc = sum(v>thres) - nPeak;
        nPeak = sum(v>thres);
%         plot(newIntens);
    end
    I = v>thres;
    pLoc = pkLoc(I);
    pInt = pkInts(I); 
    if isShow
        figure;
        plot(subplot(211),intens); hold on; 
        if sum(I)>0
            scatter(pLoc,pInt,'filled');
        end
        plot(subplot(212),pkLoc,v); hold on; 
        if sum(I)>0
            scatter(pLoc,v(I),'filled');
        end
        linkaxes([subplot(211),subplot(212)],'x');
    end
end

