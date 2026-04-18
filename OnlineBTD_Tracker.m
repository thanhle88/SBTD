%"onlineBTD: Online Block Term Tensor Decomposition"
function [loss_online] = OnlineBTD_Tracker(X_true, X_stream,rank_L,size_slice,K)

W = 5;               % window
I = size_slice(1);
J = size_slice(2);
L = {[rank_L(1) 1 rank_L(1)];[rank_L(2) 1 rank_L(2)];[rank_L(3) 1 rank_L(3)]}  ;

% Initilization of initial tensor and noisy factor
SNR = 5;
initilization_type = 1;       % 1: from noisy factors, 2: random initilzation
t_training = round(K*0.1);
%[~,U_initial] = getnoisyFactorandSubset(fac,L,SNR,t_training,initilization_type);
U_initial = btd_rnd([I J t_training], L)'; 
T_initial = X_stream(:,:,1:t_training);

% BTD ALS for initil tensor
options.printitn  = 2;        % Show convergence progress every 2 iterations.
options.MaxIter   = 10000;    % Set convergence max iterations criteria.
options.tol       = 1e-4;     % Set convergence stoping  criteria
options.fac_type  = 1;

[facs_btd_initial,~,~,~,~] = btd_als_2(T_initial,U_initial,L,options);


%% online BTD process
itr = 0; % track time and loss
mem_online = 0;
time_online = 0;
facs_btd = facs_btd_initial;

for i = t_training+1 : W : K
    itr = itr+1;
    % fprintf('Slice number: %d\n',i);
    if(i+W-1 > K)
        W = (K-i+1); % help in case last slice size is less than W/
    end
    % create tensor from new slices.
    X_new_slice = X_stream(:,:,i:i+W-1);
    [facs_intermediate,~] = onlineBTD(X_new_slice,facs_btd,L,options.fac_type);
     facs_btd = facs_intermediate;

    % Evaluation 
    tmp = X_true(:,:,i:i+W-1);
    [loss_online(itr),~] = lossCal(tmp,facs_intermediate);
       
end

end




