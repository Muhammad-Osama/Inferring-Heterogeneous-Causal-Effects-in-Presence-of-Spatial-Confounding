function theta_hat = estimate_theta(disc_or_cont,w_hat,v_hat,s,Ns,sup, mn, mx)
%Estimates theta_hat in the min-max optimization problem in paper eq. (12)

%Input:
%dis_or_cont: whether space is continuous (Euclidean) or discrete (regions)
%w_hat: n x 1 residual of outcome 'y'
%v_hat: n x 1 residual of outcome 'z'
%s: n x 1 vector of region index in case of discrete space or cartesian 
%coordinates (n x d matrix)if space is R^d
%Ns: number of b-spline components in basis function \phi(s) in case space
%is continuous or number of regions in case space is discrete
%sup: vector of support of b-splines

%Output:
%theta_hat: d_{theta} x 1 parameter vector 


%%
if strcmp(disc_or_cont,'cont')==1
    [n,d] = size(s);
    d_theta = (Ns^d)*length(sup)+1;
else
    d_theta = Ns+1;
end

theta_hat = zeros(d_theta,1);
    
Gamma_spice = zeros(d_theta);

rho_spice   = zeros(d_theta,1);

kappa_spice = zeros(1,1);
    
U = 1; L=3;

for i=1:n
    if strcmp(disc_or_cont,'cont')==1
        phis=[];
        for j=1:length(sup)
            phis = [phis func_phi_bsplinebasis(s(i,:),mn,mx,Ns,sup(j))];
        end
    else
        phis = zeros(1,Ns);
        phis(s(i)) = 1;
    end
      
    alpha_i = [0,v_hat(i).*phis];
      
    [ theta_hat, Gamma_spice, rho_spice, kappa_spice  ] = func_newsample_covlearn( w_hat(i), alpha_i, Gamma_spice, rho_spice, kappa_spice, theta_hat, i, L, d_theta, U );
end

