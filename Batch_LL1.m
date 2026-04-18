function [Ut,er_Xt] = Batch_LL1(X_stream,X_true,rank_L,size_slice,OPTS)
 
%%
if isfield(OPTS,'window') % block of observations at each time
     W = OPTS.window;
else W = 5;   
end
%%     
I   = size_slice(1);
K   = size_slice(2);
R   = length(rank_L);    

%% Tracking ...

% Initalization using batch LL1 method
t_train  = W; %sum(rank_L);
X_train  = []; 
for ii = 1 : t_train
     X_train(:,ii,:) = X_stream{1,ii};  
end
U_train = ll1(X_train,rank_L);
Ut = U_train;
% for r = 1 : R
%     Ut{1,r}{1,2} = randn(W,rank_L(r));
% end

At = []; 
Ct = []; 
Bt = [];

for r = 1 : R
    At = [At U_train{1,r}{1,1}];
    Ct = [Ct U_train{1,r}{1,3}];
    Bt = [Bt U_train{1,r}{1,2}];
end

t = t_train+1;  
j = 1;
while t+W-1 <= length(X_stream)
    %% Data Collecting ...
    j
    idx = t : t+W-1; 
    Xt  = [];  Xt_true = [];
    for ii = t : t+W-1
        Xt(:,ii-t+1,:)      = X_stream{1,ii};
        Xt_true(:,ii-t+1,:) = X_true{1,ii};   % for evalation 
    end     
    %% Processing ...
    Ut = ll1(Xt,Ut,rank_L);
    %% Evaluation Performance    
    Er_Xt    = ll1res(Xt_true,Ut,rank_L);
    er_Xt(j) = norm(Er_Xt(:)) / norm(Xt_true(:));

    %% 
    t = t + W;
    j = j + 1;
end

end