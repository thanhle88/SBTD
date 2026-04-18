%"onlineBTD: Online Block Term Tensor Decomposition"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BTD_ALS BTD by Alternate least squares.
%   [Facs_final,MemUsed ] = onlineBTD(Xnew,facOld,L)  
% Inputs:
%        Xnew: New imcoming slice
%		 facOld: Old factors available
%        L: block rank
%        fac_type: 1 for LU factorization and 2 for QR factorization
% Output:
%		Facs_final: updated factors
%       MemUsed: total memory used in process 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Facs_final,MemUsed ] = onlineBTD(Xnew,facOld,L,fac_type)
        R = max(size(L));
        size_core = size(facOld{1}{end});
        facNew = facOld;
        % update temporal mode
          T = tens2mat(Xnew,3);% T_JKxI
          G = pkron(facNew,3);
          tmp = T*pinv(G);
          K = size(tmp,2);
          idx = [1:size_core(3):K K+1];
          for r = 1:R
           facNew{r}{3}= tmp(:,idx(r):idx(r+1)-1);  
          end
          
          % update non-temporal mode A
          
          for n =1:2
              T = tens2mat(Xnew,n);% T_JKxI
              Gn = pkron(facNew,n);
              tmp = T*pinv(Gn);
              K = size(tmp,2);
              idx = [1:size_core(n):K K+1];
              for r = 1:R
               facNew{r}{n}= tmp(:,idx(r):idx(r+1)-1) + facOld{r}{n}; 
              end
          end
          
          % update core
          facNew = updateCoreTensor(Xnew,facNew,L,fac_type);
          Facs_final = facNew;
           for r =1:R
               % Facs_final{r}{3}=[facOld{r}{3};facNew{r}{3}];  % edited by
               % Thanh
               Facs_final{r}{3}=[facNew{r}{3}]; 
           end
     mem_elements = eval('whos');
     memory_array=zeros(size(mem_elements,1),1);
     if size(mem_elements,1) > 0 
            for i = 1:size(mem_elements,1) 
                    memory_array(i) = mem_elements(i).bytes;
            end
            MemUsed = sum(memory_array);
    else
            MemUsed = 0;
    end 
        
end
 
function fac=updateCoreTensor(X,fac,L,fac_type)
      % unfold the tensors in all three modes
       Ra=size(L,1);
       X2=X(:);
           
       
       for r = 1:Ra
             if r ==1
                J = kron(kron(fac{r}{3},fac{r}{2}),fac{r}{1});
             
             else
                J = [J, kron(kron(fac{r}{3},fac{r}{2}),fac{r}{1})];
        
             end
       end
       M1 = geninv(J);
         if fac_type == 1
           Jt=J';
           [L1,U1] =lu(Jt*J);
           x = U1\(L1\Jt);
           M=  x*X2 ;
        else
            [Q,R]=qr(J,0);           
            M =  R\(Q'*X2); 
        end
         K = size(M,1);
         ll = size(fac{1}{end}(:),1);
         q = [1:ll:K K+1];
         for r = 1:Ra
          fac{r}{end}= reshape(M(q(r):q(r+1)-1,1),L{r}) ;  
         end
          
end