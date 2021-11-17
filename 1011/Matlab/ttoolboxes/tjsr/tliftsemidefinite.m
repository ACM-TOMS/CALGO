function [M, JSR] = tliftsemidefinite(M, k, varargin)
% tliftsemidefinite(M, k, [options])
% Computes semidefinite liftings of the matrices M
% Notes: The function may abort if there is not enough memory
%
% Input:
%   M                   the matrices
%   k                   lifting factor
%
% Options:
%   'verbose',val       verbose-level
%
% Output:
%   Mk      the lifted matrices.
%           JSR(Mk)^(1/2k)=JSR(M)
%   JSR     obtained interval for the JSR from lifting
%
% E.g.: tliftsemidefinite({[1 1;0 1] [1 0; 1 0]}, 2)
%
% Written by: Jungers, JSR-toolbox
% Changed by: tommsch, 2018 
%
% See also: tliftproduct
%
% REFERENCES
%  [1]  V.D.Blondel and Y.Nesterov,
%         "Computationally efficient approximations of the joint spectral radius",
%         SIAM J. of Matrix Analysis, 27(1):256-272, 2005
%  [2]  V.D.Blondel, R. Jungers and V.Protasov,
%         "On the complexity of computing the capacity of codes that avoid forbidden 
%          difference patterns" 
%         IEEE Transactions on Information Theory, 52(11):5122-5127, 2006
%  [3]  R.Jungers, 
%         "The Joint Spectral Radius: Theory and Applications" 
%         Vol. 385 section 2.3.6 in Lecture Notes in Control and Information
%         Sciences, Springer-Verlag. Berlin Heidelberg, June 2009

% Initialization
nM = length(M);
dim = size(M{1}, 1);

verbose=parsem({'verbose','v'},varargin,0);
JSR=[0 inf];

% Iteration
    for i = 1:k
       
        numEntries = length(M{1})^4; % number of entries in Mj below  % Correction

        numBytesNeeded = numEntries * 8; % supposes double (8 Bytes) 

        if numBytesNeeded >  0.85*tavailable_memory % limit to 85% of physical memory
            error('Too much memory would be needed.')
        end
        
        vprintf('Starting lift number %2.0f \n', i,'imp',[1 verbose]);
        
        MM = cell(nM, 1);
        nn = dim*(dim+1)/2;
        S = zeros(nn, nn);
        keep = true(dim^2, 1);
        
        % Build representative matrices
        for j = 1:nM
            Mj = kron(M{j}, M{j});
            for k = 1:dim
                keep(k*dim-dim+k+1:k*dim) = 0;
                Mj(:, k*dim+k:dim:dim*dim-dim+k) = Mj(:, k*dim+k:dim:dim*dim-dim+k) + Mj(:, k*dim-dim+k+1:k*dim);
            end
            MM{j} = Mj(keep, keep);
            S = S + MM{j};
        end
        JSR(2) = max(abs(eig(S)))^(1/2^i);
        JSR(1) = 1/(nM^(1/2^i)) * JSR(2);
        
        vprintf('Iteration %d - current bounds: [%.15g, %.15g] \n', i, JSR(1), JSR(2),'imp',[2 verbose]);
        
        M = MM;
        dim = nn;
    
    end
    
   
    
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.