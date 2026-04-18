%"onlineBTD: Online Block Term Tensor Decomposition"
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BTD_ALS BTD by Alternate least squares.
%   [U,output] = btd_als(T,U0,L,varargin) computes R terms U{r} corresponding to a
%   block term decomposition of the N-th order tensor T by minimizing
%   0.5*frob(T-btdgen(U))^2. Each term U{r} is a cell array of N factor
%   matrices U{r}{n}, followed by a core tensor U{r}{N+1}. The algorithm is
%   initialized with the terms matrices U0{r}. The structure output returns
%   additional information:
%
%      output.loss  - The loss of the approximation algorithm.
%      output.fval -  BTD objective function value.
%      output.totTime  - Total CPU time of algorithm.
%      output.MemUsed  - Total memory used during algorithm.
%   btd_als(T,U0,L,options) may be used to set the following options:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [U,loss,fval,totTime,MemUsed] = btd_als(T,U0,L,varargin)
%ALS algorithm for decomposition in rank-(L,M,N) terms
    %% Extract number of dimensions and norm of X.
    N = size(L{1},2)+1;
    %% Set algorithm parameters from input or by using defaults
    params = inputParser;
    params.addParameter('tol',1e-4,@isscalar);
    params.addParameter('MaxIter',10000,@(x) isscalar(x) & x > 0);
    params.addParameter('printitn',2,@isscalar);
    params.addParameter('fac_type',2,@isscalar);
    params.parse(varargin{:});
   
    %% Copy from params object
    reslosschangetol = params.Results.tol;
    maxiters = params.Results.MaxIter;
    printitn = params.Results.printitn;
    fac_type = params.Results.fac_type;
    totTime=0;
    %% Set up and error checking on initial guess for U.
    U = U0;
    [loss, fval] = lossCal(T,U);
    % fprintf(' Loss = %e, Fval = %e\n', loss,fval);
    for iter = 1:maxiters
        st = tic;
        reslossold = loss(end) ;
     
        for d=1:N
            
            if d==N
             U = updateCoreTensor(T,U,L,fac_type);
 
            else
             U = updateFactor(T,U,d);
            end
        end
        totTime =totTime + toc(st);
        [loss(iter+1), fval(iter+1)] =    lossCal(T,U); 
           
        reslosschange = abs(reslossold - loss(iter+1));
          if mod(iter,printitn) == 0
             % fprintf(' Iter %2d: loss = %e loss-delta = %7.1e, fval = %e\n', iter,  loss(iter+1), reslosschange,fval(iter+1));
          end
          % Check for convergence
        if (iter > 1) && (reslosschange < reslosschangetol)
            flag = 0;
        else
            flag = 1;
        end
         % Check for convergence
        if (flag == 0)
            break;
        end
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
       R=size(L,1);
       X2=X(:);
        for r = 1:R
             if r ==1
                 H = {fac{r}{3},fac{r}{2},fac{r}{1}};
                 J = struct_kron(H);
             else
                H = {fac{r}{3},fac{r}{2},fac{r}{1}};
                J = [J, struct_kron(H)];
            end
        end
        M = geninv(J)*X2;
%         if fac_type == 1
%            Jt=J';
%            [L1,U1] =lu(Jt*J);
%            x = U1\(L1\Jt);
%            M=  x*X2 ;
%         else
%             [Q,R]=qr(J,0);           
%             M =  R\(Q'*X2); 
%         end
        K = size(M,1);
        ll = size(fac{1}{end}(:),1);
        q = [1:ll:K K+1];
        for r = 1:R
          fac{r}{end}= reshape(M(q(r):q(r+1)-1,1),L{r});  
        end
          
 end
 
function z = updateFactor(X,z,d)
         % unfold the tensors in all three modes
          R=size(z,2);
          T = tens2mat(X,d);% T_JKxI
          G = pkron(z,d);
          M = T*pinv(G);
          K = size(M,2);
          ll = size(z{1}{d},2);
          L = [1:ll:K K+1];
          for r = 1:R
              z{r}{d}= M(:,L(r):L(r+1)-1);  
          end
     
 end

