function [y,z,s] = scm(n,d)
%Generates simulated data a/c to the following structure causal model
%(scm):
% space:     s \sim Uniform[0,100]
% exposure:  z = \alpha*s + N(0,1)
% outcome:   y = \tau(s)*z + \beta*s + N(0,0.2^2)
% where      \tau(s) = cos(2*pi*s/200) in 1D & \tau(s) =
%            cos(2*pi*s1/200)*cos(2*pi*s2/200) in case of 2D

%Input:
%n : nos. of points
%d: dimension of space (1D or 2D)

%Output:
%y: n x 1 outcome vector
%z: n x 1 exposure vector
%s: n x d space coordinates
%%
alpha = 0.5.*ones(d,1); 
beta = ones(d,1); 
tau = {@(a) cos(2*pi.*a/200), @(a) cos(2*pi.*a/200)};
  
s = 100*rand(n,d);
tau_s = ones(n,1);
for i=1:Nc
    for d = 1:d
         tau_s(:,i) = tau_s(:,i).*tau{i}(s(:,d));
    end
end
v = randn(n,1);
z = s*alpha + v;
u = 0.2.*randn(n,1);
f = sum(tau_s.*z,2) + s*beta;
y = f + u ;