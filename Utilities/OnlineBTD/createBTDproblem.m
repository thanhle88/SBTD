function [X,fac]= createBTDproblem(I,J,K,L)
size_tens = [I J K];
size_core=L';
fac = btd_rnd(size_tens,size_core);
X   = btdgen(fac);
end
