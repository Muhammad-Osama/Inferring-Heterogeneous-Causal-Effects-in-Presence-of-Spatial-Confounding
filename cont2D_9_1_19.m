%This script tests the performance of LS method when estimating effect in
%continuous space-2D when the nusisance function beta(s) is
%high dimensional

%% set seed
clc; 
stream = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(stream);
savedState=stream.State;
%% space grid
s = 0:0.4:10;
[s1, s2] = meshgrid(s);
s_grid = [s1(:) s2(:)];
N = size(s_grid,1);
%% data generating process
%nuisance function 
mn = [-5 -5]; mx = [15 15]; M = 10; sup = [0.2 0.4 0.85]; 
Psi = zeros(N,M^2*length(sup)); lambda = randn(M^2*length(sup),1); lambda(101:end) = 0;
for n = 1:N
    psi = [];
    for s=1:length(sup)
       psi = [psi func_phi_bsplinebasis( s_grid(n,:), mn, mx, M,sup(s))];
    end
    Psi(n,:) = psi;
end
beta = Psi*lambda;
figure;
contourf(s1,s2,reshape(beta,size(s1)));
exp_z_s_true = 0.5.*beta; 

z = exp_z_s_true + randn(N,1);

tau = cos(2*pi.*s1(:)/20+2*pi.*s2(:)/20);

exp_y_s_true = tau.*exp_z_s_true+beta;

y = tau.*z+ beta + 0.2.*randn(N,1);
%% basis for effect
Ns = 10;
%%
[theta_ls,~,Phi] = LS_cont_2D(N,Ns,sup,mn,mx,s_grid,z,y);
tau_est_ls = Phi*theta_ls;

%% ROSCE method

[v_hat,exp_z_s,w_v] = Exp_z_given_x('cont',z,s_grid,Ns,sup,mn,mx);

[r_hat, exp_y_s, w_r] = Exp_y_given_x('cont',y,s_grid,Ns,sup,mn,mx);

c1 = min(v_hat); c2 = max(v_hat);

[tau_est_rosc, w,~] = estimate_gamma('cont',r_hat,v_hat,c1,c2,s_grid,Ns,sup,mn,mx);

%% Bootstrapping
tic;
Nboot = 3000; qaunt94 = floor(Nboot*0.94); qaunt4 = floor(Nboot*0.04);
theta_rosc = zeros((Ns^2)*length(sup)+1,Nboot); theta_ls = zeros(Ns^2*length(sup),Nboot);
for i=1:Nboot
    i
   [~,idx] = datasample(y,N,'Replace',true);
   
   theta_rosc(:,i) = scm_train('cont', y(idx),z(idx),s_grid(idx,:),Ns,sup,sup,mn,mx);
   
   [theta_ls(:,i),~,~] = LS_cont_2D(N,Ns,sup,mn,mx,s_grid(idx,:),z(idx),y(idx));

end
tau_boot_rosc = Phi*theta_rosc(2:end,:);
tau_boot_ls = Phi*theta_ls;
%%
 alpha = 0.1; q1 = floor(Nboot*(1-alpha/2)); q2 = floor(Nboot*alpha/2);

 temp = sort(tau_boot_rosc,2); tau_ci_rosc = [2*tau_est_rosc-temp(:,q1), 2*tau_est_rosc-temp(:,q2)];
 
 temp = sort(tau_boot_ls,2); tau_ci_ls = [2*tau_est_ls-temp(:,q1), 2*tau_est_ls-temp(:,q2)];
%%
tau_est_ls_mat = reshape(tau_est_ls,size(s1));
tau_est_rosc_mat = reshape(tau_est_rosc,size(s1));
tau_true_mat = reshape(tau,size(s1));
tau_ci_ls_low = reshape(tau_ci_ls(:,1),size(s1));
tau_ci_ls_upp = reshape(tau_ci_ls(:,2),size(s1));
tau_ci_rosc_low = reshape(tau_ci_rosc(:,1),size(s1));
tau_ci_rosc_upp = reshape(tau_ci_rosc(:,2),size(s1));
%%
figure;
%ax1 = subplot(1,3,1);
contourf(s1,s2,reshape(tau,size(s1)));
xlabel('$s_1$','interpreter','Latex');ylabel('$s_2$','interpreter','Latex')
xlim([0 10]); ylim([0 10]);
caxis([-1 1]);colorbar;

%ax2 = subplot(1,3,2);
figure;
contourf(s1,s2,reshape(tau_est_ls,size(s1))); xlabel('$s_1$','interpreter','Latex');ylabel('$s_2$','interpreter','Latex')
xlim([0 10]); ylim([0 10]);
caxis([-1 1]);colorbar;

%ax3 = subplot(1,3,3);
figure;
contourf(s1,s2,reshape(tau_est_rosc,size(s1))); xlabel('$s_1$','interpreter','Latex');ylabel('$s_2$','interpreter','Latex')
xlim([0 10]); ylim([0 10]);
caxis([-1 1]);colorbar;

%hLink3 = linkprop( [ax1,ax2,ax3], 'CLim' );
%colorbar;
%%
s = 0:0.4:10;
figure;
h = fill([s,fliplr(s)],[diag(tau_ci_rosc_low)',fliplr(diag(tau_ci_rosc_upp)')],[0.5 0.5 0.5]);
set(h,'FaceAlpha',0.3);hold on;grid on; ylim([-2 2])
%p1 = plot(s,diag(tau_est_rosc_mat),'b','LineWidth',2);
p2 = plot(s,diag(tau_true_mat),'r--','LineWidth',2);
%legend([h p2 p1],{'$95\%~CI~our~model$','$\gamma(s)$','$\hat{\gamma}(s)$'},'Interpreter','latex');
xlabel('$s$','interpreter','Latex'); ylabel('$Effect$','interpreter','Latex');
%%
s = 0:0.4:10;
figure;
h = fill([s,fliplr(s)],[diag(tau_ci_ls_low)',fliplr(diag(tau_ci_ls_upp)')],'y');
set(h,'FaceAlpha',0.3);hold on;grid on;
%p1 = plot(s,diag(tau_est_ls_mat),'b','LineWidth',2);
%p2 = plot(s,diag(tau_true_mat),'r--','LineWidth',2);
%legend([h p2 p1],{'$95\%~CI~our~model$','$\gamma(s)$','$\hat{\gamma}(s)$'},'Interpreter','latex');
%xlabel('$s$','interpreter','Latex'); ylabel('$Effect$','interpreter','Latex');