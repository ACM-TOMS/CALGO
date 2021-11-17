function funcs = Build_funcs(fx, type, a, b, rho, tau)
% Builds the subfunctions needed to run the different parts of the problem.
% Running Build_funcs instantiates all subfunctions listed below with the
% correct parameters of the problem in question.

% Reference:
% Ratnanather, J. T., Kim, J. H., Zhang, S., Davis, A. M. J., and Lucas, 
% S. K. 2013. Algorithm XXX: IIPBF, a MATLAB toolbox for infinite integral 
% of products of two Bessel functions. To appear in ACM Transactions on 
% Mathematical Software.

% Input:
% fx: function handle
% type: character either 'JJ', 'JY', or 'YY', from IIPBF.m
% a,b: integer order of the two bessel functions
% rho: parameter multiplier of argument of first bessel function
% tau: parameter multiplier of argument of second bessel function.
%
%   Reference: Ratnanather, J. T., Kim, J. H., Zhang, S., Davis, A. M. J.,
%              and Lucas, S. K. 2011. Algorithm: IIPBF, a MATLAB toolbox
%              for infinite integral of products of two Bessel functions.
%              Pending review from ACM Transactions on Mathematical
%              Software
% Output:
% funcs: 5x1 cell of function handles

% Author: Jung Hun Kim
% Date: August 30, 2012

    function y = fxBx_typeJJ(x)
        y = fx(x).*besselj(a,rho*x).*besselj(b,tau*x);
    end

    function y = fxBx_typeJY(x)
        y = fx(x).*besselj(a,rho*x).*bessely(b,tau*x);
    end

    function y = fxBx_typeYY(x)
        y = fx(x).*bessely(a,rho*x).*bessely(b,tau*x);
    end

    function y = fxh1_typeJJ(x)
        y = fx(x).*jj_h1(x);
    end

    function y = fxh2_typeJJ(x)
        y = fx(x).*jj_h2(x);
    end

    function y = fxh1_typeJY(x)
        y = fx(x).*jy_h1(x);
    end

    function y = fxh2_typeJY(x)
        y = fx(x).*jy_h2(x);
    end

    function y = fxh1_typeYY(x)
        y = fx(x).*yy_h1(x);
    end

    function y = fxh2_typeYY(x)
        y = fx(x).*yy_h2(x);
    end

    function y = jj_h1(x)
        y=0.5.*(besselj(a,rho.*x).*besselj(b,tau.*x) -...
            bessely(a,rho.*x).*bessely(b,tau.*x));
    end

    function y = jj_h2(x)
        y=0.5.*(besselj(a,rho.*x).*besselj(b,tau.*x) +...
            bessely(a,rho.*x).*bessely(b,tau.*x));
    end

    function y = jy_h1(x)
        y=0.5.*(besselj(a,rho.*x).*bessely(b,tau.*x) +...
            bessely(a,rho.*x).*besselj(b,tau.*x));
    end

    function y = jy_h2(x)
        y=0.5.*(besselj(a,rho.*x).*bessely(b,tau.*x) -...
            bessely(a,rho.*x).*besselj(b,tau.*x));
    end

    function y = yy_h1(x)
        y=-0.5.*(besselj(a,rho.*x).*besselj(b,tau.*x) -...
            bessely(a,rho.*x).*bessely(b,tau.*x));
    end

    function y = yy_h2(x)
        y=0.5.*(besselj(a,rho.*x).*besselj(b,tau.*x) +...
            bessely(a,rho.*x).*bessely(b,tau.*x));
    end

switch type
    case 'JJ'
        funcs{1} = @fxBx_typeJJ;
        funcs{2} = @fxh1_typeJJ;
        funcs{3} = @fxh2_typeJJ;
        funcs{4} = @jj_h1;
        funcs{5} = @jj_h2;
    case 'JY'
        funcs{1} = @fxBx_typeJY;
        funcs{2} = @fxh1_typeJY;
        funcs{3} = @fxh2_typeJY;
        funcs{4} = @jy_h1;
        funcs{5} = @jy_h2;
    case 'YY'
        funcs{1} = @fxBx_typeYY;
        funcs{2} = @fxh1_typeYY;
        funcs{3} = @fxh2_typeYY;
        funcs{4} = @yy_h1;
        funcs{5} = @yy_h2;
    otherwise
        error('type is ''JJ'', ''JY'', or ''YY'' only')
end

end
