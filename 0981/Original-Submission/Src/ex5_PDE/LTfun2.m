function U = LTfun2(X,Y,S)
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

%   U(h,k) = U(X(h),Y(h),S(k))
    U = zeros(numel(X),numel(S));
	for k=1:numel(S)
		U(:,k) = 4/S(k)^2 + ( X.*(X-1) + Y.*(Y-1) )/S(k);
	end
end
