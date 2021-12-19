function tau_CI = compute_CI(disc_or_cont,y,z,s,tau_hat,Ns,sup,mn,mx)
%Uses 500 bootstrap iterations to compute the 95% confidence interval (CI)
%of the effect estimate

%Input:
%y: n x 1 outcome vector 
%z: n x 1 exposure vector
%s: n x 1 vector of region index in case of discrete space or cartesian 
%coordinates (n x d matrix)if space is R^d
%Ns: number of b-spline components in basis function \phi(s) in case space
%is continuous or number of regions in case space is discrete
%sup: vector of support of b-splines

%Output:
%tau_CI: [] x 2 matrix lower (1st column) and upper bound (2nd column) 
%of the 95% CI 


%%
Nboot = 500; qaunt99 = floor(0.99*Nboot); qaunt4 = floor(0.04*Nboot);

tau_hat_boot = [];

n = length(y);

for i=1:Nboot
    i
    [~,idx] = datasample(y,n,'Replace',true);
     
    w_hat = Exp_given_s(disc_or_cont,y(idx),s(idx,:),Ns,sup,mn,mx);

    v_hat = Exp_given_s(disc_or_cont,z(idx),s(idx,:),Ns,sup,mn,mx);

    theta_hat = estimate_theta(disc_or_cont,w_hat,v_hat,s(idx,:),Ns,sup, mn, mx);

    [temp,~] = estimate_effect(disc_or_cont, theta_hat, Ns, sup, mn , mx);
    
    tau_hat_boot = [tau_hat_boot temp];
end 

dtau_hat = tau_hat_boot - tau_hat; dtau_hat = sort(dtau_hat,2); tau_CI = [tau_hat-dtau_hat(:,qaunt99), tau_hat-dtau_hat(:,qaunt4)];