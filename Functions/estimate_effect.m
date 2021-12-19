function [tau_hat,se] = estimate_effect(disc_or_cont, theta_hat, Ns, sup, mn , mx)
%Compute effect estimate \widehat{tau}(s) from learned parameters theta_hat

%%
if strcmp(disc_or_cont,'disc')==1
    %effect in different regions of space
    tau_hat = theta_hat(2:end);
else
    %dimension of space
    d = length(mn);
    %define fine grid to evaluate effect on
    if d==1
        %1d
        se = linspace(mn,mx,100)';
    else
       %2d 
       [s1,s2] = meshgrid(linspace(mn,mx,25));
       se = [s1(:) s2(:)];
    end
    %nos. of points
    n = size(se,1);  tau_hat = zeros(n,1);
    for i=1:n
        phis=[];
        for j=1:length(sup)
            phis = [phis func_phi_bsplinebasis(se(i,:),mn,mx,Ns,sup(j))];
        end
        %effect over continuous space
        tau_hat(i) = phis*theta_hat(2:end);
    end
        
end