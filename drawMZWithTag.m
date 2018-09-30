function drawMZWithTag( ha,mz,intens,num )
    bar(ha,mz,intens,100);
    ha.NextPlot = 'add';
    [~,Is] = sort(intens,'descend');
    for m = 1:min([num,length(mz)])
        text(ha,mz(Is(m)),intens(Is(m)),sprintf('%.4f',mz(Is(m))),'FontSize',8);
    end
    ha.NextPlot = 'replace';
end

