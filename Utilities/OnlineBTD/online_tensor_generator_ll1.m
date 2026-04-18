function [X_stream,U_stream] = online_tensor_generator_ll1(size_slice,num_slice,rank_L,epsilon)

I            = size_slice(1);
K            = size_slice(2);
T            = num_slice;
R            = length(rank_L);  
L_Component  = sum(rank_L);

A = randn(I,L_Component);
C = randn(K,R);

for t = 1 : T
      
    bt = randn(1,L_Component);
    Xt = zeros(I,K);
    Ut = cell(1,R);
    index = 1;

    for r = 1 : R
        
        ll = index + rank_L(r) - 1;
        A_r  = A(:,index:ll);
        bt_r = bt(:,index:ll);
        c_r  = C(:,r);
        Xt   = Xt + (A_r * bt_r') * c_r';
        index = ll+1;  
        
        Ut{1,r}{1,1} = A_r;
        Ut{1,r}{1,2} = bt_r;
        Ut{1,r}{1,3} = c_r;
        Ut{1,r}{1,4} = eye(rank_L(r));
    end
   
    % save
    X_stream{1,t} = Xt;
    U_stream{1,t} = Ut;
    
    %% time-varying
    NA = randn(I,L_Component);
    NC = randn(K,R);
    A = A + epsilon(t) * NA;
    % C = C + epsilon(t) * NC/norm(NC);
end

end



