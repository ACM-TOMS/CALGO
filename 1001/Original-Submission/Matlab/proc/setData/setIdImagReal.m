%% setIdImagReal
% Define transformation operators to identify real and complex
% vectors/matrices.
%
%% Syntax
%
%   seti = setIdImagReal(seti)
%
%% Description
% |seti = setIdImagReal(seti| adds to struct seti 
% seti.S and seti.T to identify real and complex vectors/matrices
% (S maps C to R x R and T maps R x R to C) as well as
% auxiliary functions seti.R and seti.I to extract the real and imaginary
% part.
%
% * The input struct seti does not need fields, but has to exist to add the
% mentioned new fields.
%
%% Input Arguments
%
% * seti    :   structure array (no fields are needed)
%
%% Output Arguments
%
% * seti.S  :   function to identify complex with real vectors/matrices
%               (input: complex vector/matrix; output: real vector/matrix)
% * seti.R  :   extracts the originally real (R) part of the output of seti.S
% * seti.I  :   extracts the originally imag (I) part of the output of seti.S
% * seti.T  :   function to identify real with complex vectors/matrices
%               (input: real vector/matrix (R x R); output: complex vector/matrix)
%
%% More About:
%
% The transformation operators are, see Section 4.5 in [1]:
%
% * $\texttt{seti.S}
%  = T_{\bf{C} \to \bf{R}^2} : \bf{C} \to \bf{R} \times \bf{R}$, $\quad$
% $T_{\bf{C} \to \bf{R}^2}(x) = (\mathrm{real}(x),\,  \mathrm{imag}(x))$.
%
% * $\texttt{seti.T} 
% = T_{\bf{R} \to \bf{C}} : \bf{R} \times \bf{R} \to \bf{C}$, $\quad$
% $T_{\bf{R} \to \bf{C}} = 
%  y^{\mathrm{real}} + \mathrm{i}\,y^{\mathrm{imag}} 
%  \mathrm{\ where\ } y = (y^{\mathrm{real}},y^{\mathrm{imag}})$.
%
% For $y = (y^{\mathrm{real}},y^{\mathrm{imag}})$ the auxiliary functions
% are defined as:
% 
% * $\texttt{seti.R}(y) = y^{\mathrm{real}}$,
% * $\texttt{seti.I}(y) = y^{\mathrm{imag}}$.
%
%
%% References
% * [1] Florian B&uuml;rgel, Kamil S. Kazimierski, and Armin Lechleiter. A sparsity regularization and total variation based computational framework for the inverse medium problem in scattering. _Journal of Computational Physics_, 339:1-30, 2017.
%
%% See Also
% * <setGeomSim.html>
% * <pda.html>
%
%% Code
function seti = setIdImagReal(seti)
% Identification of imaginary vector with real (and inverse), i.e.
% for real and complex vectors: S: C -> R x R and T: R x R -> C.
%
% Notation:
% hs = h_S = [hr;hi] = [real(h); imag(h)];
% hr = h_R, hi = h_I
% h = hz (complex)

seti.S = @(hz) [real(hz); imag(hz)]; % S: C -> R x R
% More about seti.S:
% real and imag part as vector (V): S: h |-> [hr; hi] (C -> R x R) 
% (vector or matrix possible)

seti.R = @(hs) hs(1:size(hs,1)/2,:,:,:); % real (R) part
% (real part is stored in the first half of real vector hs)

seti.I = @(hs) hs(size(hs,1)/2+1:size(hs,1),:,:,:); % imag (I) part 
% (imaginary part is stored in the second part of real vector hs)

% More about seti.R and seti.I:
% * Note that size(hs,1) is faster than end 
%   (when you have matrices with 3 dimensions).
% * Note that R and I have 4 input arguments 
%   because the matrix gu to store grad(u) has 4 dimensions in case of 3D.

seti.T = @(hs) seti.R(hs) + 1i.*seti.I(hs); % T: R x R -> C
% More about seti.T:
% real(h) + i imag(h): T: hs |-> real + 1i imag (R x R -> C) 
% (vector or matrix possible)

% Connection of S, T, R and I:
% hs = S(hz); hz = T(hz); hr = R(hs); hi = I(hs)

end
