function [gam_aprx, w, brk_pnts] = estimate_gamma(disc_or_cont,r_hat,v_hat,c1,c2,x,Nx,sup,mn,mx,varargin)

%Nx: nos. of regions in case of DISCRETE
%v_hat: N x Nc matrix where N is the nos. of points and Nc is nos. of
%covariates
[~,Nc] = size(v_hat);

if strcmp(disc_or_cont,'cont')==1
    [~,D] = size(x);
end

if nargin>10
    q = varargin{:};
    Nv = q;
    brk_pnts = break_points_eta(v_hat,q);
else
    Nv = Nc;
    brk_pnts = [];
end

if strcmp(disc_or_cont,'cont')==1
    P = (Nx^D)*length(sup)*Nv+1;
else
    P = Nx*Nv+1;
end

w = zeros(P,1);
    
Gamma_spice = zeros(P);

rho_spice   = zeros(P,1);

kappa_spice = zeros(1,1);
    
U = 1; L=3;

for n=1:size(x,1)
    
    if strcmp(disc_or_cont,'cont')==1
        phix=[];
        for i=1:length(sup)
            phix = [phix func_phi_bsplinebasis(x(n,:),mn,mx,Nx,sup(i))];
        end
    else
        phix = zeros(1,Nx);
        phix(x(n)) = 1;
    end
      
      if nargin>10
        phi_v = eta_basis(v_hat(n),brk_pnts);
      else
        phi_v = zeros(1,Nc);  
        for j=1:Nc
            phi_v(j)  = v_hat(n,j)*(v_hat(n,j)>=c1(j) && v_hat(n,j)<=c2(j)) + c1(j)*(v_hat(n,j)<c1(j)) + c2(j)*(v_hat(n,j)>c2(j));
        end
      end
      
      alpha_n = [0,kron(phi_v,phix)];
      
      %Regression between space, residuals and outcome
      [ w, Gamma_spice, rho_spice, kappa_spice  ] = func_newsample_covlearn( r_hat(n), alpha_n, Gamma_spice, rho_spice, kappa_spice, w, n, L, P, U );
end

if nargin<=10
    gam_aprx = zeros(size(x,1),Nc);
for c=1:Nc
    for n=1:size(x,1)
   
     if strcmp(disc_or_cont,'cont')==1
        phix=[];
        for i=1:length(sup)
            phix = [phix func_phi_bsplinebasis(x(n,:),mn,mx,Nx,sup(i))];
        end
        gam_aprx(n,c) = phix*w((c-1)*(Nx^D)*length(sup)+2:c*(Nx^D)*length(sup)+1) ;
     else
        phix = zeros(1,Nx);
        phix(x(n)) = 1;
        gam_aprx(n,c) = phix*w((c-1)*Nx+2:c*Nx+1) ;
     end
 
        
    
    end
end
else
    
    gam_aprx =[];
    
end
