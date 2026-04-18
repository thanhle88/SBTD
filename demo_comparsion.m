%"onlineBTD: Online Block Term Tensor Decomposition"
clear all;clc;close all;
addpath(genpath('.'));
run_path

%% 

size_slice   = [20 30];    % size of data streams
N            = 300;        % number of data blocks of observations\
W            = 5;          % block of observations at each time 
T            = N*W;        % total number of data 
rank_L       = [2 3 4];    % rank of block terms
R            = length(rank_L); 
noise_fac    = 1e-3;       % noise 
epsilon      = zeros(T,1); % time-varying factor
epsilon(1:W:T) = 1e-4;     % time-varying factor
epsilon(2*W)   = 1;        % abrupt change
  
epsilon2        = zeros(T,1); % time-varying factor
epsilon2(1:W:T) = 1e-4;     % time-varying factor
 


I = size_slice(1);
J = size_slice(2);


%% Data Generation
disp('stream data generating ... ')
disp(' + generate true data  ... ')

[X_true,~]      = online_tensor_generator_ll1(size_slice,T,rank_L,epsilon);
[X_true_LL1,~]  = online_tensor_generator_ll1(size_slice,T,rank_L,epsilon2);

for ii = 1 : T
   X_stream{1,ii} = X_true{1,ii} + noise_fac*randn(size_slice);
   X_stream_LL1{1,ii} = X_true_LL1{1,ii} + noise_fac*randn(size_slice);

end

disp(' +  SBTD  ... ')

OPTS = [];
[~,erX_SBTD] = SBTD(X_stream,X_true,rank_L,size_slice,OPTS);

OPTS_LL1.forgetting_factor = 0.99;
[~,erX_LL1]  = SBTD(X_stream_LL1,X_true_LL1,rank_L,size_slice,OPTS_LL1);

disp(' +  OnlineBTD  ... ')
rank_L = [5 5 5]; N_ex = 3;
erX_OnlineBTD_cell = cell(N_ex,1); 

for n_exp  =  1 : N_ex
    
    [X_true_OBTD,~]  = online_tensor_generator_ll1(size_slice,T,rank_L,epsilon);
    X_stream_OBTD    = zeros(I,J,T);
    X_stream_OBTD_true  = zeros(I,J,T);
    
    for ii = 1 : T
        X_stream_OBTD(:,:,ii)      = X_true_OBTD{1,ii}  + noise_fac*randn(size_slice);
        X_stream_OBTD_true(:,:,ii) = X_true_OBTD{1,ii};
    end
    erX_OnlineBTD_nxp = OnlineBTD_Tracker(X_stream_OBTD_true,X_stream_OBTD,rank_L,[I,J],T);
    erX_OnlineBTD_cell{n_exp,1} = erX_OnlineBTD_nxp;
end

erX_OnlineBTD = erX_OnlineBTD_cell{1,1};
for n_exp  =  2 : N_ex
    erX_OnlineBTD = erX_OnlineBTD + erX_OnlineBTD_cell{n_exp,1};
end
erX_OnlineBTD = erX_OnlineBTD/N_ex;
erX_OnlineBTD = [zeros(1,28) erX_OnlineBTD];
erX_OnlineBTD = erX_OnlineBTD(2:end);
%% plot and compare loss

figure; hold on;
plot(erX_LL1,'k-','LineWidth',2,'MarkerSize',10);
plot(erX_OnlineBTD,'b-','LineWidth',2,'MarkerSize',10);
plot(erX_SBTD,'r-','LineWidth',2,'MarkerSize',10);
h = gca;
lgd = legend('\texttt{LL1} (Batch)','\texttt{OnlineBTD}','\texttt{SBTD} (Proposed)');
lgd.FontSize = 18;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);
xlabel('Time Index - $t$','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('RE($\mathcal{Y}_{tr},\mathcal{Y}_{es}$)','interpreter','latex','FontSize',18,'FontName','Times New Roman');


set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');

set(gca,'YScale', 'log','FontSize', 20);

grid on;
box on;





