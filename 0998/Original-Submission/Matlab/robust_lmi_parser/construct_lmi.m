function LMIs = construct_lmi(varargin)
% LMIs = construct_lmi(Term,ineq)
%Function to construct the LMIs "Term 'ineq' 0"
%Inputs: Term    -> Struct (Term{i,j}) containing the i-th row and j-th 
%                    column of the LMI.
%        ineq    -> Signal of the inequality. Can be '>','<','>=','<='. If
%                    ineq is a different string, such string will be the
%                    label of the resultant matrix.

if (nargin == 2)
    LMIs = construct_lmi_terms(varargin{1},varargin{2});
elseif (nargin == 3)
    LMIs = construct_lmi_terms(varargin{1},varargin{2},varargin{3});
elseif (nargin == 4) %(Term,vertices,vecpoly,'<')
    LMIs = construct_lmi_terms(varargin{1},varargin{2},varargin{4});
else
    error('construct_lmi is an obsolete function');
    LMIs = [];
end


return