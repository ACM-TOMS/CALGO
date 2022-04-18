function a = kronm(b, c, varargin)
% a = kronm( A1, ..., An )
% Kronecker product of arbitrary many matrices.
%
% A multi-argument version of kron, defined by 
%
%   kronm(A1, A2, A3) = kron(kron(A1,A2), A3)
%
% and so on.  This is identical to the behaviour of
% Octave's built-in kron function, but not Matlab's.
    
    a = kron(b, c);
    
    for i = 1:(nargin-2)
        a = kron(a, varargin{i});
    end
    
end