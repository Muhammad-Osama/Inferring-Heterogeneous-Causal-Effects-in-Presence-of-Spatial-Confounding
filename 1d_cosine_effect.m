%% set the seed
clc; 
stream = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(stream);
savedState=stream.State;

%% Generate data
n = 500; d = 1; 
[y,z,s]= simulated_data(n,d);
%true effect
tau = @(u) cos(2*pi.*u./20);

%% basis parameters
%continuous space
disc_or_cont = 'cont';

Ns = 20; sup = [0.65 1.15]; mn = -5; mx = 20;

%% estimate effect

w_hat = Exp_given_s(disc_or_cont,y,s,Ns,sup,mn,mx);

v_hat = Exp_given_s(disc_or_cont,z,s,Ns,sup,mn,mx);

theta_hat = estimate_theta(disc_or_cont,w_hat,v_hat,s,Ns,sup, mn, mx);

[tau_hat,se] = estimate_effect(disc_or_cont, theta_hat, Ns, sup, mn , mx);

%% confidence interval

tau_CI = compute_CI(disc_or_cont,y,z,s,tau_hat,Ns,sup,mn,mx);

%%
figure;
h = fill([se;flipud(se)], [tau_CI(:,1);flipud(tau_CI(:,2))],[0.5 .5 0.5]);
set(h,'FaceAlpha',0.3);
hold on;
plot(se,tau(se),'r--',se,tau_hat,'b','LineWidth',1.5);
xlabel('$s$','interpreter','Latex'); 
legend({'$95\%~CI$','$\tau(s)$','$\widehat{\tau(s)}$'},'interpreter','Latex');
grid on;
xlim([0, 10])
ylim([-1.5 1.5])
