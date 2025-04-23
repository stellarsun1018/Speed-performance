% mmsigma and hitrates needs to be n-by-1, i.e. in single columns, and of
% the same size
function pxtSizes = tSizeGen(mmsigma,hitrates,pixellength)
    pxsigma = mmsigma ./ pixellength;
    p_edpt = @(x,sig) normpdf(x,0,sig);
    q = NaN(length(pxsigma),300);
    pxtSizes = NaN(length(hitrates),length(mmsigma));
    for i = 1:length(pxsigma)
        for j = 1:300
            q(i,j) = integral(@(x) p_edpt(x,pxsigma(i)),-j,j);
        end
        [~,pxtSizes(:,i)] = min(abs(repmat(q(i,:),length(hitrates),1) - hitrates),[],2);
    end

    return 
end