function [U,U_sub] = getnoisyFactorandSubset(fac,L,SNR,initial_batch_size,initilization_type)
    if initilization_type ==1
        U = noisy(fac,SNR);
        U_sub = U;
      for r = 1:size(L,1)
         U_sub{r}{3}=U_sub{r}{3}(1:initial_batch_size,:);
      end
    else
        I = size(fac{1}{1},1);J = size(fac{1}{2},1);K = size(fac{1}{3},1);
        U = btd_rnd([I J K], L)'; 
        U_sub = btd_rnd([I J initial_batch_size], L)'; 
    end
end