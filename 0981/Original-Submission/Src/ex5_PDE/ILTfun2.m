function u = ILTfun2(X,Y,t)
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                EXAMPLE 5

        Application of Talbot's method to solve
        the following PDE problem

                u_t  = u_xx  + u_yy,            0 < x,y < 1,   t>0
                u(x,y,0+) = x*(x-1) + y*(y-1)
                u(0,y,t) = u(1,y,t) = 4*t + y*(y-1)
                u(x,0,t) = u(x,1,t) = 4*t + x*(x-1)

        The analytical solution is

                u(x,y,t) = 4*t + x*(x-1) + y*(y-1)

        and Laplace Transform is

                          4    x*(x-1) + y*(y-1)
                U(x,s) = --- + -----------------
                         s^2           s

        where s = 0 is a double pole.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}

%   u(h,k) = u(X(h),Y(h),t(k))
    u = zeros(numel(X),numel(t));
    for h=1:numel(X)
        for k=1:numel(t)
            u(h,k) = 4*t(k) + X(h)*(X(h)-1) + Y(h)*(Y(h)-1);
        end
    end
end
