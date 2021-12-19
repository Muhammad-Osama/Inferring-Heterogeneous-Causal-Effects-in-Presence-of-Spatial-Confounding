function [Phi,cen] = func_phi_bsplinebasis( x, mn, mx, M,support)
%% Phi(x) - using cubic Bspline with rectangular boundaries
% x row vector R^d
% mn row vector R^d
% mx row vector R^d

%% Initialize
D   = size(x,2);
Phi = ones(1,(M^D));
cen = zeros(D,M);

L_vect = mx-mn; 

%% Construct Phi
j_vec = ones(1,D); %index

for c = 1:M^D
    
    %Product
    for k = 1:D
        %stp = L_vect(k)/M;
        cen(k,:) = linspace(mn(k),mx(k),M);
        [~,i] = min(abs(cen(k,:))); cen(k,i) = 0;
        
        %cen(k,end) = mx(k);
        
        %z    =   3/(support*L_vect(k)); %for quadratic bspline
        
        z    =   4/(support*L_vect(k)); %for cubic bspline
        
        %z    =   4/(support); %for cubic bspline
        
        Phi(c) = Phi(c) * cubic_bspline_at_point(x(k),cen(k,:),j_vec(k),z);
    end
    
    %Index update
    j_vec(1) = j_vec(1) + 1;
    if D>1
        for k = 2:D
            if (mod(j_vec(k-1),M+1) == 0)
                j_vec(k-1) = 1;
                j_vec(k)   = j_vec(k) + 1;
            end
        end
    end

end

end

