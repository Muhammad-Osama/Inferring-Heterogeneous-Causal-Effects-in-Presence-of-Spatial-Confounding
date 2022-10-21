function [theta_ls,lambda_ls,Phi] = LS_cont_2D(N,Ns,sup,mn,mx,x,z,y)

Phi = zeros(N,(Ns^2)*length(sup));
for n=1:N
   phi = [];
   for s = 1:length(sup)
       phi = [phi func_phi_bsplinebasis(x(n,:),mn,mx,Ns,sup(s))];
   end
   Phi(n,:) = phi;
end

A = [z.*Phi Phi]; A = A+0.0005.*eye(size(A));
coeff = (A'*A)\(A'*y);
theta_ls  = coeff(1:Ns^2*length(sup)); lambda_ls = coeff(Ns^2*length(sup)+1:end);
