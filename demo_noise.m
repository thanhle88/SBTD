clear; clc; 

run_path;

size_slice   = [20 30];    % size of data streams
N            = 500;        % number of data blocks of observations\
W            = 5;          % block of observations at each time 
T            = N*W;        % total number of data 
rank_L       = [2 3 4];    % rank of block terms
R            = length(rank_L); 
noise_fac    = 1e-3;       % noise 
epsilon      = zeros(T,1); % set of time-varying factors
epsilon(1:W:T) = 1e-4;     % time-varying factors 
epsilon(2*W)   = 1;        % time-varying factors   
epsilon(300*W) = 1;        % abrupt change
 

%% Data Generation
disp('stream data generating ... ')
disp(' + generate true data  ... ')
[X_true,~]  = online_tensor_generator_ll1(size_slice,T,rank_L,epsilon);

disp(' + add noise  ... ')
noise_fac1 = 1e-3;
noise_fac2 = 1e-2;
noise_fac3 = 1e-1;

for ii = 1 : T
   X_stream1{1,ii} = X_true{1,ii} + noise_fac1*randn(size_slice);
   X_stream2{1,ii} = X_true{1,ii} + noise_fac2*randn(size_slice);
   X_stream3{1,ii} = X_true{1,ii} + noise_fac3*randn(size_slice);
end 


%% Tracking ...
disp('tracking process ...')

disp(' + Noise 1')
OPTS = [];
[~,erX_SBTD_1] = SBTD(X_stream1,X_true,rank_L,size_slice,OPTS);

disp(' + Noise 2')
[~,erX_SBTD_2] = SBTD(X_stream2,X_true,rank_L,size_slice,OPTS);

disp(' + Noise 3')
[~,erX_SBTD_3] = SBTD(X_stream3,X_true,rank_L,size_slice,OPTS);


%% PLOT RESULTS
disp('plotting ...')

makerSize    = 14;
numbMarkers  = 50;
LineWidth    = 2;

k = 1;
M = length(erX_SBTD_1); 
figure; hold on;
d1 = semilogy(1:k:M,erX_SBTD_1(1,1:k:end),...
    'linestyle','-','color','r','LineWidth',LineWidth);
d11 = plot(1:100:M,erX_SBTD_1(1,1:100:end),...
 'marker','h','markersize',makerSize,...
   'linestyle','none','color','r','LineWidth',LineWidth);
d12 = semilogy(1:1,erX_SBTD_1(1,1),...
    'marker','h','markersize',makerSize,...
    'linestyle','-','color','r','LineWidth',LineWidth);

d2 = semilogy(1:k:M,erX_SBTD_2(1,1:k:end),...
    'linestyle','-','color','b','LineWidth',LineWidth);
d21 = plot(1:100:M,erX_SBTD_2(1,1:100:end),...
 'marker','o','markersize',makerSize,...
   'linestyle','none','color','b','LineWidth',LineWidth);
d22 = semilogy(1:1,erX_SBTD_2(1,1),...
    'marker','o','markersize',makerSize,...
    'linestyle','-','color','b','LineWidth',LineWidth);

d3 = semilogy(1:k:M,erX_SBTD_3(1,1:k:end),...
    'linestyle','-','color','k','LineWidth',LineWidth);
d31 = plot(1:100:M,erX_SBTD_3(1,1:100:end),...
 'marker','d','markersize',makerSize,...
   'linestyle','none','color','k','LineWidth',LineWidth);
d32 = semilogy(1:1,erX_SBTD_3(1,1),...
    'marker','d','markersize',makerSize,...
    'linestyle','-','color','k','LineWidth',LineWidth);

axis([0 N 1e-4 1])
lgd = legend([d12 d22 d32],'$\sigma_n = 10^{-3}$','$\sigma_n = 10^{-2}$','$\sigma_n = 10^{-1}$');
lgd.FontSize = 18;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);
xlabel('Time Index - $t$','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('RE($\mathcal{Y}_{tr},\mathcal{Y}_{es}$)','interpreter','latex','FontSize',18,'FontName','Times New Roman');
% 
h=gca;
set(gca,'YScale', 'log');
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:200:N,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 20);
grid on;
box on;

