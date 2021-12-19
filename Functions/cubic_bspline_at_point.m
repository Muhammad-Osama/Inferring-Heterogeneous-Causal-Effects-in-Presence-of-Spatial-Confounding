function f=cubic_bspline_at_point(x,cen,j,z)
% x scalar value
% cen a 1xM vector of centers of different bsplines
% j index for different bsplines
% z scaling factor to adjust support of bspline

%% Initialize
b=z*(cen(j)-2/z);

%% Evaluate bspline at x
    
    if ((0+b)/z<=x)&&(x<(1+b)/z)
        f=1/6*(z*x-b)^3;
    elseif ((1+b)/z<=x)&&(x<(2+b)/z)
        f=-1/2*(z*x-b)^3 + 2*(z*x-b)^2 - 2*(z*x-b) + 2/3;
    elseif ((2+b)/z<=x)&&(x<(3+b)/z)
        f=1/2*(z*x-b)^3 - 4*(z*x-b)^2 + 10*(z*x-b) - 22/3;
    elseif ((3+b)/z<=x)&&(x<(4+b)/z)
        f=-1/6*(z*x-b)^3 + 2*(z*x-b)^2 - 8*(z*x-b) + 32/3;
    else
        f=0;
    end
