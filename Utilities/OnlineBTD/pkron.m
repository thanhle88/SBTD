function AB = pkron(U,d)
R = size(U,2);stIdx=1;

for r = 1:R
    H = U{r}(1:3);S=U{r}{end};
    if d == 1
        A =  (H{3});
        B =  (H{2});
    end
    if d ==2
           A =  (H{3});
           B =  (H{1});     
     end
    if d == 3
       A =  (H{2});
       B =  (H{1});
    end
    
 
%     kr = tens2mat(H{end},d)*kron(A,B)';
      tmp=kron_onlineBTD(A,B);
      e = tens2mat(S,d);
      kr = e*transpose(tmp);
   
% 
%     st1=tic; kr = mtkronprod_fmt(S,H,d);ft1=toc(st1);
%     
%     
% %     st1=tic; kr = mtkronprod_fast(S,H,d);
% %     
%    st2=tic;     kr = mtkronprod(S,H,d,'H');ft2=toc(st2);
    nextidx = size(kr,1);
    endIdx=stIdx + nextidx-1;
    AB(stIdx:endIdx,:) = kr;
    stIdx=endIdx+1;
end

 
end
 