function [ alpha,beta,gam,pairs ] = msSetAlign( spLoc,rpLoc,thres )
    %[ alpha,beta,gam ] = msSetAlign( spLoc,rpLoc,thres )
    %   alpha = A - B
    %   beta = AxB
    %   gam = B-A
    spLoc = spLoc(:);
    rpLoc = rpLoc(:);
    Ls = length(spLoc);
    Lr = length(rpLoc);
    beta = 0;
    isAligned_s = zeros(Ls,1);
    pairs = zeros(min([Ls,Lr]),2);
    for m = 1:Lr
        I = find(and(abs(spLoc-rpLoc(m))<=thres,isAligned_s==0));
        if isempty(I)
            continue;
        else
            isAligned_s(I(1)) = 1;
            beta = beta+1;
            pairs(beta,:) = [spLoc(I(1)),rpLoc(m)];
        end
    end
    alpha = Ls-beta;
    gam = Lr-beta;
    pairs((beta+1:end),:)=[];   
end

