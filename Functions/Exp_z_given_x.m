function [residual,exp, w_v] = Exp_z_given_x(dis_or_cont,z,x,Nx,sup,mn,mx,varargin)

%Nx : would be equal to nos. of regions in case of discrete

if strcmp(dis_or_cont,'cont')==1
    [~,D] = size(x);
end
%using least square

% if nargin>2
%     w_v = varargin{:};
%     if length(x)==1
%         x=x.*ones(length(z),1);
%     end
%     B = [ones(length(x),1) x];
%     residual = z-B*w_v;
%     exp = []; 
% else
%     B = [ones(length(x),1) x]; 
%     w_v = (B'*B+0.1.*eye(size(B'*B)))\(B'*z);
%     exp = B*w_v;
%     residual = z-exp;
% end

if nargin>7
    w_v = varargin{:};
    
    if size(x,1) == 1
        if strcmp(dis_or_cont,'cont')==1
            phix=[];
            for i=1:length(sup)
                phix = [phix func_phi_bsplinebasis(x,mn,mx,Nx,sup(i))];
            end
        else
            phix = zeros(1,Nx);
            phix(x) = 1;
        end
        alpha_n = [1,phix];
        residual = z - alpha_n*w_v;
        exp = [];
    else
       residual = zeros(length(z),1); 
       
       for n=1:length(z)
           if strcmp(dis_or_cont,'cont')==1
                phix=[];
                for i=1:length(sup)
                    phix = [phix func_phi_bsplinebasis(x(n,:),mn,mx,Nx,sup(i))];
                end
           else
                phix = zeros(1,Nx);
                phix(x(n)) = 1;
           end
            alpha_n = [1,phix];
            residual(n) = z(n) - alpha_n*w_v;
       end
       exp = [];
    end
else
    
    if strcmp(dis_or_cont,'cont')==1
        P = (Nx^D)*length(sup)+1;
    else
        P = Nx+1;
    end
     
    w_v = zeros(P,1);
    
    Gamma_spice = zeros(P);

    rho_spice   = zeros(P,1);

    kappa_spice = zeros(1,1);
    
    U = 1; L=3;

    for n=1:length(x)
        if strcmp(dis_or_cont,'cont')==1
            phix=[];
            for i=1:length(sup)
                phix = [phix func_phi_bsplinebasis(x(n,:),mn,mx,Nx,sup(i))];
            end
        else
            phix = zeros(1,Nx);
            x(n)
            phix(x(n)) = 1;
        end
        
      alpha_n = [1,phix];
 
      [ w_v, Gamma_spice, rho_spice, kappa_spice  ] = func_newsample_covlearn( z(n), alpha_n, Gamma_spice, rho_spice, kappa_spice, w_v, n, L, P, U );
    end
    
    residual = zeros(length(z),1); 

    exp = zeros(length(z),1);

    for n=1:length(z)
      
      if strcmp(dis_or_cont,'cont')==1
            phix=[];
            for i=1:length(sup)
                phix = [phix func_phi_bsplinebasis(x(n,:),mn,mx,Nx,sup(i))];
            end
        else
            phix = zeros(1,Nx);
            phix(x(n)) = 1;
        end
      
      alpha_n = [1,phix];
      
      exp(n) = alpha_n*w_v;
       
      residual(n) = z(n)-exp(n);
    end
end