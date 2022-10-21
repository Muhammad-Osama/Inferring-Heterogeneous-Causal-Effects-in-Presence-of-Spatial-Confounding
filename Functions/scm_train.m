function [w, w_v, w_r, gam_aprx, exp_y_given_x,r_hat,v_hat] = scm_train(varargin)
%% INPUTS
% discrete_cont. : discrete or continous space
% y : N x 1  response vector. 'N' being the nos. of data points
% z : N x 1 covariate/exposure vector
% x : N x D matrix of spatial coordinates in case of continuous space. 'D'
% being dimension or N x 1 vector of region indicator containing values
% from 1,...,Nx for total Nx number of regions
% Nx : number of bspline basis in case of continous or number of regions in
% case of discrete
% sup : vector equal to number of supports to be used in bspline for
% continous case
% mn : D x 1 vector of minimum coordinate in each spatial dimension
% mx : D x 1 vector of maximum coordinate in each spatial dimension
 y = varargin{2}; z = varargin{3}; s = varargin{4}; discrete_cont = varargin{1}; Nx = varargin{5}; 
 sup = varargin{6}; sup2 = varargin{7}; mn = varargin{8}; mx = varargin{9};
%% OUTPUTS

%%
[N,Nc] = size(z); %nos. of covariates

v_hat = zeros(N,Nc); 

for i=1:Nc
    [v_hat(:,i),~,w_v] = Exp_z_given_x(discrete_cont,z(:,i),s,Nx,sup,mn,mx);
end

[r_hat, exp_y_given_x, w_r] = Exp_y_given_x(discrete_cont,y,s,Nx,sup,mn,mx);

c1 = min(v_hat,[],1); c2 = max(v_hat,[],1);

[gam_aprx, w,~] = estimate_gamma(discrete_cont,r_hat,v_hat,c1,c2,s,Nx,sup2,mn,mx); % on training data
