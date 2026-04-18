function [loss,fval] = lossCal(T,z)
        fval = z{1}{1}*mtkronprod(z{1}{end},z{1}(1:end-1),1,'H');
        for r = 2:length(z)
            fval = fval+z{r}{1}*mtkronprod(z{r}{end},z{r}(1:end-1),1,'H');
        end
        fval = fval-reshape(T,size(fval));
        fval = 0.5*(fval(:)'*fval(:));
        loss = frob(btdres(T,z))/frob(T);
    
end