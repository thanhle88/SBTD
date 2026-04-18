function T = getpartTensor(U,st,ed)
 tmpfac=U;
 
 for r = 1:size(U,2)
     tmpfac{r}{3}=U{r}{3}(st:ed,:);
 end
 T=btdgen(tmpfac);
end