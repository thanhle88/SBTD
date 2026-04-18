clear; clc; close all;

run_path;

size_slice   = [20 30];    % size of data streams
N            = 500;        % number of data blocks of observations\
W            = 5;          % block of observations at each time 
T            = N*W;        % total number of data 
rank_L       = [2 3 4];    % rank of block terms
R            = length(rank_L); 
noise_fac    = 1e-3;       % noise 
 



%% Data Generation
disp('stream data generating ... ')
disp(' + generate true data  ... ')


epsilon1        = zeros(T,1); % time-varying factor
epsilon1(1:W:T) = 1e-4;      % time-varying factor
epsilon1(2*W)   = 1;         % abrupt change
epsilon1(300*W) = 1;    

epsilon2        = zeros(T,1); % time-varying factor
epsilon2(1:W:T) = 1e-3;     % time-varying factor
epsilon2(2*W)   = 1;        % abrupt change
epsilon2(300*W) = 1;    


epsilon3        = zeros(T,1); % time-varying factor
epsilon3(1:W:T) = 1e-2;     % time-varying factor
epsilon3(2*W)   = 1;        % abrupt change
epsilon3(300*W) = 1;    

epsilon4        = zeros(T,1); % time-varying factor
epsilon4(1:W:T) = 1e-1;     % time-varying factor
epsilon4(2*W)   = 1;        % abrupt change
epsilon4(300*W) = 1;    


[X_true1,~]  = online_tensor_generator_ll1(size_slice,T,rank_L,epsilon1);
[X_true2,~]  = online_tensor_generator_ll1(size_slice,T,rank_L,epsilon2);
[X_true3,~]  = online_tensor_generator_ll1(size_slice,T,rank_L,epsilon3);
[X_true4,~]  = online_tensor_generator_ll1(size_slice,T,rank_L,epsilon4);

disp(' + add noise  ... ')
noise_fac = 1e-3;

for ii = 1 : T
   X_stream1{1,ii} = X_true1{1,ii} + noise_fac*randn(size_slice);
   X_stream2{1,ii} = X_true2{1,ii} + noise_fac*randn(size_slice);
   X_stream3{1,ii} = X_true3{1,ii} + noise_fac*randn(size_slice);
   X_stream4{1,ii} = X_true4{1,ii} + noise_fac*randn(size_slice);

end 


%% Tracking ...
disp('tracking process ...')

disp(' + time-varying 1')
OPTS = [];
[~,erX_SBTD_1] = SBTD(X_stream1,X_true1,rank_L,size_slice,OPTS);

disp(' + time-varying 2')

[~,erX_SBTD_2] = SBTD(X_stream2,X_true2,rank_L,size_slice,OPTS);

disp(' + time-varying 3')

[~,erX_SBTD_3] = SBTD(X_stream3,X_true3,rank_L,size_slice,OPTS);

disp(' + time-varying 4')

[~,erX_SBTD_4] = SBTD(X_stream4,X_true4,rank_L,size_slice,OPTS);

disp('plotting ...')




%% PLOT RESULTS
makerSize    = 14;
numbMarkers  = 50;
LineWidth    = 2;

k = 1;
M = size(erX_SBTD_1,2); 
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


d4 = semilogy(1:k:M,erX_SBTD_4(1,1:k:end),...
    'linestyle','-','color','g','LineWidth',LineWidth);
d41 = plot(1:100:M,erX_SBTD_4(1,1:100:end),...
 'marker','^','markersize',makerSize,...
   'linestyle','none','color','g','LineWidth',LineWidth);
d42 = semilogy(1:1,erX_SBTD_4(1,1),...
    'marker','^','markersize',makerSize,...
    'linestyle','-','color','g','LineWidth',LineWidth);


axis([0 N 0.5e-4 5])
lgd = legend([d12 d22 d32 d42],'$\varepsilon = 10^{-4}$','$\varepsilon = 10^{-3}$', ...
                               '$\varepsilon = 10^{-2}$','$\varepsilon = 10^{-4}$');
lgd.FontSize = 18;
set(lgd, 'Interpreter', 'latex', 'Color', [0.95, 0.95, 0.95]);
xlabel('Time Index - $t$','interpreter','latex','FontSize',13,'FontName','Times New Roman');
ylabel('RE($\mathcal{Y}_{tr},\mathcal{Y}_{es}$)','interpreter','latex','FontSize',18,'FontName','Times New Roman');
% 
h=gca;
set(gca,'YScale', 'log');
set(h,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle','-','MinorGridLineStyle','-','FontName','Times New Roman');
set(h,'Xtick',0:100:N,'FontSize',16,'XGrid','on','YGrid','on','GridLineStyle',':','MinorGridLineStyle','none',...
    'FontName','Times New Roman');
set(h,'FontSize', 20);
grid on;
box on;

