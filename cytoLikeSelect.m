function [I,t,isIn] = cytoLikeSelect(setObj,r1,r2,isShow)
    if ~exist('isShow','var')
        isShow = 0;
    end
    t = cell2mat(setObj.getFieldByName('sampleTime'));
    c1 = setObj.getMZRange(r1)+10;
    c2 = setObj.getMZRange(r2)+10;
    C1 = log10(c1(:)); C2 = log10(c2(:));
    
    hf = figure;
    [f,~] = ksdensity([C1,C2],[C1,C2]);
    scatter(C1,C2,10,f,'filled');
    a = impoly;
    pXpY = a.getPosition();
    isIn = inpolygon(C1,C2,pXpY(:,1),pXpY(:,2));
    close(hf);
    
    if isShow
        figure('Position',[0,0,600,600]);
        scatter(C1(~isIn),C2(~isIn),15,'filled'); hold on;
        scatter(C1(isIn),C2(isIn),15,'filled');box on;

        figure('Position',[0,0,1200,500]);
        plot(subplot(211),t,c1);
        hold on;
        scatter(t(isIn),c1(isIn),10,'filled');

        plot(subplot(212),t,c2);
        hold on;
        scatter(t(isIn),c2(isIn),10,'filled');
        linkaxes([subplot(211),subplot(212)],'x');
    end
    t = t(isIn);
    I = 1:length(c1);
    I = I(isIn);
end
