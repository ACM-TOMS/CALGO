function [norm, prob] = polyliftedNorm(V,z,tol,logfile,opts)
% 
% [NORM, PROB] = POLYLIFTEDNORM(V,z,TOL)
%
%      Computes the polytope norm of point z in the lifted space: 
%      V{i} is the ith vertex of the polytope norm, 
%      and |z|= inv(sup{l:lz \leq \sum{V{i}\lambda_i}})
%                s.t.   \sum \lambda_i = 1
%                            \lambda \geq 0
%
%      TOL is the tolerance for the SDP solver
%
%      PROB = 1 signals some problem in the resolution of
%               the SDP program, allows to be conservative 
%               with v in jsr_norm_conitope
%      PROB = 0 if a feasible solution has been found
%
% [NORM, PROB] = POLYLIFTEDNORM(V,z,logfile,opts)
%       Is similar but also uses :
%
%      LOGFILE as a file identifier in which it prints the possible
%              messages
%
%      OPTS    is the options structure used with the algorithm, for use
%              with msg in order to see if messages should be printed in
%              the command line
% 
% Requires SeDuMi [by default] ( http://sedumi.ie.lehigh.edu/ )

if (nargin<4)
    logfile = -1;
    opts.verbose=2;
end

% Parameters
m = length(V);
n = size(V{1},1);
opts.SDPsolver.solverOptions.fid = 0;
opts.SDPsolver.solverOptions.bigeps = 1e-3;
tolSolver = tol;


% ---- SDP program construction 
%
% variables : x = [l \lambda_1 ... \lambda_m vec(P)]'
% 
%  Sedumi form :
%     min c'*x
%        s.t. A*x=b
% 
 
A = sparse([],[],[],1+n^2,1+m+n^2);
A(1,2:2+m-1) = ones(1,m); % Constaint \sum \lambda_i = 1

% Constraint P = sum{V(:,:,i)\lambda_i} - l z
A(2:end,1) = -vec(sparse(z));
for i = 1:m
    A(2:end,i+1) = vec(sparse(V{i}));
end
A(2:end,m+2:end) = spdiags(-ones(n^2,1),0,n^2,n^2);

b = spdiags([1],0,1+n^2,1);
c = spdiags([-1],0,1+m+n^2,1);

% Specify cone on which one optimizes the variables :
K.l = m+1;
K.s = n;

% ---- Resolution
[result,dual,feas,stopFlag] = solve_semi_definite_program(A,b,c,K,tolSolver,opts);

    if (feas == 0),
        msg(logfile,opts.verbose>1,' solve_semi_definite_program:  Failure - Infeasible problem');
        prob = 1;
    elseif (stopFlag > 0),
        msg(logfile,opts.verbose>1,' solve_semi_definite_program: Failure - Numerical problems : stopFlag = %d',stopFlag);
        prob = 1;
    else
        msg(logfile,opts.verbose>1,' solve_semi_definite_program: Success - Found feasible solution');
        prob = 0;
    end

    if (~prob)
        norm = 1/result(1);
    else
        norm = 1/result(1);
    end

end