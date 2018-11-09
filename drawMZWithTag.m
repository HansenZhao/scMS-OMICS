function drawMZWithTag( ha,mz,intens,num )
    TARGET_WIDTH = 0.2;
    minDis = min(abs(mz(2:end)-mz(1:(end-1))));
    factor = TARGET_WIDTH/minDis;
    bar(ha,mz,intens,factor);
    ha.NextPlot = 'add';
    ha.XLim = [min(mz),max(mz)];
    [~,Is] = sort(intens,'descend');
    for m = 1:min([num,length(mz)])
        text(ha,mz(Is(m)),intens(Is(m)),sprintf('%.4f',mz(Is(m))),'FontSize',8);
    end
    ha.NextPlot = 'replace';
end

