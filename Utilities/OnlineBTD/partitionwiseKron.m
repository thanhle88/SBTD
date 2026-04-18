function AB = partitionwiseKron(U,L,d)
R=size(L,1);
AB=cell(1,R);
for r=1:R
    if d==1
       AB{r}= mtkronprod(U{r}{end},U{r}(1:3),d,'H') ;
    elseif d==2
       % AB{r}= kron(U{r}{1},U{r}{3});
        AB{r}= mtkronprod(U{r}{end},U{r}(1:3),d,'H');
    elseif d==3
       % AB{r}= kron(U{r}{1},U{r}{2});
        AB{r}= mtkronprod(U{r}{end},U{r}(1:3),d,'H') ;
    else
        error('Not valid input');
    end
end



end