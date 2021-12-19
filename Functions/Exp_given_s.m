function res_hat = Exp_given_s(dis_or_cont,q,s,Ns,sup,mn,mx)
%computes the residual \widehat{v}=z-\widehat{E[z|s]} for q = z or 
%\widehat{w} = y-E[y|s] q = y by solving the regularized problems in paper

%Input:
%dis_or_cont: whether space is continuous (Euclidean) or discrete (regions)
%q: n x 1 vector of exposure/outcome variable
%s: n x 1 vector of region index in case of discrete space or cartesian 
%coordinates (n x d matrix)if space is R^d
%Ns: number of b-spline components in basis function \phi(s) in case space
%is continuous or number of regions in case space is discrete
%sup: vector of support of b-splines

%Output:
%res_hat: n x 1 vector of residuals
%%
if strcmp(dis_or_cont,'cont')==1
    [n,d] = size(s);
    P = (Ns^d)*length(sup)+1; 
else
    P = Ns+1;
end

wres = zeros(P,1);
    
Gamma_spice = zeros(P);

rho_spice   = zeros(P,1);

kappa_spice = zeros(1,1);
    
U = 1; L=3;

for i=1:n
   if strcmp(dis_or_cont,'cont')==1
      %define basis vector for continuous space
      phis=[];
      for j=1:length(sup)
          phis = [phis func_phi_bsplinebasis(s(i,:),mn,mx,Ns,sup(j))];
      end
   else
      %define basis vector for discrete space
      phis = zeros(1,Ns);
      phis(s(i)) = 1;
   end
        
   alpha_i = [1,phis];
 
   [wres, Gamma_spice, rho_spice, kappa_spice  ] = func_newsample_covlearn( q(i), alpha_i, Gamma_spice, rho_spice, kappa_spice, wres, i, L, P, U );
end

%residual \widehat{v}
res_hat = zeros(n,1);

%\widehat{E[z|s]}
exp_hat = zeros(n,1);

for i=1:n
      if strcmp(dis_or_cont,'cont')==1
         %define basis vector for continuous space
         phis=[];
         for j=1:length(sup)
             phis = [phis func_phi_bsplinebasis(s(i,:),mn,mx,Ns,sup(j))];
         end
      else
         %define basis vector for discrete space
         phis = zeros(1,Ns);
         phis(s(i)) = 1;
      end
      
      alpha_i = [1,phis];
      
      exp_hat(i) = alpha_i*wres;
       
      res_hat(i) = q(i)-exp_hat(i);
end
