function [Ut,er_Xt] = SBTD(X_stream,X_true,rank_L,size_slice,OPTS)

%%
if isfield(OPTS,'forgetting_factor') % forgetting factor
     beta = OPTS.forgetting_factor;
else beta = 0.9;
end
if isfield(OPTS,'regular_parameter') % forgetting factor
     rho = OPTS.regular_parameter;
else rho = 1e-4;
end
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
t_train  = sum(rank_L);
X_train  = [];
for ii = 1 : t_train
    X_train(:,ii,:) = X_stream{1,ii};
end
U_train = ll1(X_train,rank_L);
Ut = U_train;

At = [];
Ct = [];
Bt = [];

for r = 1 : R
    At = [At U_train{1,r}{1,1}];
    Ct = [Ct U_train{1,r}{1,3}];
    Bt = [Bt U_train{1,r}{1,2}];
end

% Auxiliary matrices for recursive update
S_A = zeros(sum(rank_L));
D_A = zeros(I,sum(rank_L));
S_C = zeros(R);
D_C = zeros(K,R);

t = t_train+1;
j = 1;
while t+W-1 <= length(X_stream)
    %% Data Collecting ...
    idx = t : t+W-1;
    Xt  = [];  Xt_true = [];
    for ii = t : t+W-1
        Xt(:,ii-t+1,:)      = X_stream{1,ii};
        Xt_true(:,ii-t+1,:) = X_true{1,ii};   % for evalation
    end
    
    %% Processing ...
    HBt     = khatri_rao_b(Ct,At,rank_L);
    % HBt_inv = pinv(HBt);
    HBt_inv = (HBt'*HBt + 0.001*eye(size(HBt,2)))^(-1)*HBt';

    %% Estimate Bt
    for  ii = 1 : W
        Xt_ii  = Xt(:,ii,:);
        Xt_ii  = reshape(Xt_ii,[I K]);
        xt_ii  = Xt_ii(:);
        bt_ii  = HBt_inv *  xt_ii;
        Bt     = [Bt; bt_ii'];
    end
    
    %% Estimate At
    WAt  = khatri_rao_b(Ct,Bt(idx,:),rank_L)';
    Xt_A = ten2mat(tensor(Xt),1)';
    D_A  = beta*D_A + Xt_A' * WAt';
    S_A  = beta*S_A + WAt * WAt';
    At   = D_A * inv(S_A + rho*eye(sum(rank_L))) ;
    
    %% Estimate Ct
    WCt   = [];
    index = 1;
    for r = 1:R
        ll    = index + rank_L(r) - 1;
        AB_r  = khatrirao(At(:, index:ll),Bt(idx, index:ll));
        AB_r  = AB_r * ones(rank_L(r),1);
        index = ll + 1;
        WCt   = [WCt AB_r(:)];
    end
    Xt_C  = ten2mat(tensor(Xt),3)';
    D_C   = beta*D_C +  Xt_C' * WCt;
    S_C   = beta*S_C +  WCt' * WCt;
    C_t   = D_C * inv(S_C + rho*eye(R));
    
    %% save
    index = 1;
    for r = 1:R
        ll    = index + rank_L(r) - 1;
        Ut{1,r}{1,1} = At(:,index:ll);
        Ut{1,r}{1,2} = Bt(idx,index:ll);
        Ut{1,r}{1,3} = Ct(:,r);
        index = ll + 1;
    end
    
    %% Evaluation Performance
    Er_Xt    = ll1res(Xt_true,Ut,rank_L);
    er_Xt(j) = norm(Er_Xt(:)) / norm(Xt_true(:));
    
    %%
    t = t + W;
    j = j + 1;
end

end