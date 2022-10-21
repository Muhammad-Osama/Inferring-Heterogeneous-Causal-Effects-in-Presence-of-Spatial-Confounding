function [residual, exp, w_r] = Exp_y_given_x(disc_or_cont,y,x,Nx,sup,mn,mx,varargin)

if strcmp(disc_or_cont,'cont')==1
    [~,D] = size(x);
end

if nargin>7
    w_r = varargin{:};
    exp = zeros(length(y),1);
    
    if size(x,1) ==1
        x = repmat(x,length(y),1);
    end
    for n=1:length(y)
        
      if strcmp(disc_or_cont,'cont')==1
        phix=[];
        for i=1:length(sup)
            phix = [phix func_phi_bsplinebasis(x(n,:),mn,mx,Nx,sup(i))];
        end
      else
        phix = zeros(1,Nx);
        phix(x(n)) = 1;
      end
      
      alpha_n = [1,phix];
      
      exp(n) = alpha_n*w_r;
    end
    
    residual = [];
else

    if strcmp(disc_or_cont,'cont')==1
        P = (Nx^D)*length(sup)+1;
    else
        P = Nx+1;
    end

    w_r = zeros(P,1);
    
    Gamma_spice = zeros(P);

    rho_spice   = zeros(P,1);

    kappa_spice = zeros(1,1);
    
    U = 1; L=3;

    for n=1:length(x)
 
      if strcmp(disc_or_cont,'cont')==1
            phix=[];
            for i=1:length(sup)
                phix = [phix func_phi_bsplinebasis(x(n,:),mn,mx,Nx,sup(i))];
            end
        else
            phix = zeros(1,Nx);
            phix(x(n)) = 1;
      end
      
      alpha_n = [1,phix];
 
      [ w_r, Gamma_spice, rho_spice, kappa_spice  ] = func_newsample_covlearn( y(n), alpha_n, Gamma_spice, rho_spice, kappa_spice, w_r, n, L, P, U );
    end

    residual = zeros(length(y),1); 

    exp = zeros(length(y),1);

    for n=1:length(y)
      
      if strcmp(disc_or_cont,'cont')==1
            phix=[];
            for i=1:length(sup)
                phix = [phix func_phi_bsplinebasis(x(n,:),mn,mx,Nx,sup(i))];
            end
        else
            phix = zeros(1,Nx);
            phix(x(n)) = 1;
      end
      
      alpha_n = [1,phix];
      
      exp(n) = alpha_n*w_r;
       
      residual(n) = y(n)-exp(n);
    end

end