function sinh_bvptwp()
% sinh_bvptwp
%  compares the final mesh obtained from four solvers without continuation.
% The problem under  consideration is
%    y'' = e sinh(e y),
%   with boundary conditions
%    y(0) = 0,    y(1)=1.
%
% For this problem: 
%      Function f is vectorized
%      Analytical Jacobians are not implemented
%      RelTol is the same for both components
%      Printing statistics after solving is enabled
%
%
%       Authors:
%
%       Jeff R. Cash 
%            (Department of Mathematics, Imperial College,  London, England.)
%       Davy  Hollevoet 
%            (Vakgroep Toegepaste Wiskunde en Informatica, Universiteit Gent, Belgium.)
%       Francesca Mazzia  
%            (Dipartimento di Matematica, Universita' di Bari, Italy)
%       Abdelhameed Nagy Abdo
%            (Dipartimento di Matematica, Universit\`a di Bari, Italy)
%            (Dept. of Mathematics, Faculty of Sciences, Benha  University,Egypt)
%            
%
    options = bvptwpset('Vectorized','on',...
		'RelTol',1e-8,...
		'Linear','off',...
		'NMax',1000,...
		'Stats','on');

    % A guess for the initial mesh and the solution
    sol = bvpinit(linspace(0,1,18),[0,0]);

    e = 5.32675;

    % uses deferred correction with MIRK methods
    options = bvptwpset(options,'Solver','twpbvp_m');
    sol1 = bvptwp(@sinhODE,@sinhBC,sol,options);

    % uses deferred correction with MIRK methods and conditioning
    options = bvptwpset(options,'Solver','twpbvpc_m');
    sol2 = bvptwp(@sinhODE,@sinhBC,sol,options);

    % uses deferred correction with Lobatto methods
    options = bvptwpset(options,'Solver','twpbvp_l');
    sol3 = bvptwp(@sinhODE,@sinhBC,sol,options);

    % uses deferred correction with Lobatto methods and conditioning
    options = bvptwpset(options,'Solver','twpbvpc_l');
    sol4 = bvptwp(@sinhODE,@sinhBC,sol,options);



    %%% Nested functions, e is shared with the outer function

    function dydx = sinhODE(x,y)
        % Evaluate the ODE function (vectorized)
        dydx = [ y(2,:)
                 e*sinh(e*y(1,:)) ];
    end

    function res = sinhBC(ya,yb)
        % Evaluate the residual in the boundary conditions
        res = [ ya(1)
                yb(1)-1 ];
    end
end

