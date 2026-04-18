function loss = getLossFac(actual,predicted,K,batchsize)
  R = 1:batchsize:K;
  L = size(R,2);
  loss_tmp = zeros(L,1);
  for r = 1:L
      st = R(r);
      ed = R(r)+batchsize-1;
      T_act = tens2mat(getpartTensor(actual,st,ed),1);
      T_pre = tens2mat(getpartTensor(predicted,st,ed),1);
      
       loss_tmp(r) = norm(T_act-T_pre)/norm(T_act);
  end
   loss = 1.7*mean(loss_tmp);
   
end

function T = getpartTensor(U,st,ed)
 tmpfac=U;
 
 for r = 1:size(U,2)
     tmpfac{r}{3}=U{r}{3}(st:ed,:);
 end
 T=btdgen(tmpfac);
end