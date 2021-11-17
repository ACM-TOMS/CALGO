%% derivative
% Convenience function to compute the auxiliary matrices JA and JB to
% compute the Jacobian matrix that represents the Frechet derivative for
% finite-dimensional spaces.
%
%
%% Syntax
%
%   [JA,JB,DFFqh,seti] = derivative()
%   [JA,JB] = derivative(seti,qROI)
%   [JA,JB,DFFqh] = derivative(seti,qROI,h)
%
%% Description
% 
% |[JA,JB,DFFqh,seti] = derivative()| uses some default values (given in |seti|) 
% and computes the corresponding auxiliary matrices |JA| and |JB| such that
% the Jacobian matrix of $\mathcal{F}$ at contrast $q$ evaluted on $h$ is 
% |DFFqh| = |JA*diag(h)*JB| = $\mathcal{F}'(q)[h]$.
%
% |[JA,JB] = derivative(seti,qROI)| computes the matrices |JA| and
% |JB| corresponding to |JA*diag(h)*JB| = $\mathcal{F}'(q)[h]$.
%
% |[JA,JB,DFFqh] = derivative(seti,qROI,h) computes the matrices 
% |JA| and |JB| corresponding to |DFFqh| = |JA*diag(h)*JB| = $\mathcal{F}'(q)[h]$.
%
%
%% Example
%
%   [JA,JB,DFFqh] = derivative();
%   figure(101); imagesc(real(seti.G(DFFqh))); axis xy; axis square;
%
%
%% Input Arguments
%
% * |seti|     :    struct as described in <mimo.html>
% * |qROI|     :    contrast in ROI (region of interest)
%                   (complex matrix as vector of size seti.nROI^seti.dim x 1)
% * |h|        :    Update of contrast in ROI (region of interest)
%                   (complex matrix as vector of size seti.nROI^seti.dim x 1)
%                   (new contrast will be q := q + h)
%
%
%% Output Arguments
%
% * |JA|, |JB| :    Auxiliary matrices to compute the Jacobian matrix (derivative)
%                   (JA is a complex matrix of size seti.measNb x seti.nROI^seti.dim,
%                    JB is a complex matrix of size seti.nROI^seti.dim x seti.incNb)
% * |DFFqh|    :    Jacobian matrix of $\mathcal{F}$ at contrast $q$ evaluated on h, 
%                   i.e. |DFFqh| = $\mathcal{F}'(q)[h]$ 
%                   (complex matrix of size seti.measNb x seti.incNb)
% * |seti|     :    struct with settings
%                   (interesting if |derivative()| is called without input arguments)
%
% 
%% Best Practice
%
% This is a convenience function. To avoid unnecessary calls, use
% <mimo.html>, e.g.
%
%   [JA,JB] = mimo(seti,qROI,'jacobian');   % Auxiliary matrices JA and JB
%   DFFq = @(h) JA*diag(h)*JB;              % function handle
%   DFFqh = DFFq(h);                        % evaluate F'(q)[h]
%
%
%% More About
%
% Read <mimo.html>.
%
%% See Also
%
% * <mimo.html>
%
%% Code
%
function varargout = derivative(varargin)
%%
% Preparation for mimo:
if nargin == 0
    init;
    if ~exist('seti','var')
        seti = struct;
    end
    seti = setGeomSim(seti);
    varargout{4} = seti;
    qROI = seti.qROIexact;
    h = rand(size(qROI(:)));
elseif nargin == 2 || nargin == 3
    seti = varargin{1}; % seti... see abvoe
    qROI = varargin{2}; % contrast q on ROI as a vector e.g. 8281 x 1
    if size(qROI,1) ~= seti.nROI^seti.dim || size(qROI,2) ~= 1
        error('Make sure that dimension of input qROI is seti.nROI^seti.dim x 1.')
    end
    if nargin == 3
        h = varargin{3};
        if size(h) ~= size(qROI)
            error('Make sure that dimension of input h is the same as qROI, i.e. seti.nROI^seti.dim x 1.')
        end
    end
else
    error('Function "forward" needs 0, 2 or 3 input arguments.')
end
%%
% Call of mimo:
[JA,JB] = mimo(seti, qROI, 'jacobian'); % Auxiliary matrices JA and JB to compute the Jacobian matrix
DFFq = @(h) JA*diag(h)*JB;

if nargin == 0 || nargin == 3
    DFFqh = DFFq(h);
    varargout{3} = DFFqh;
end

varargout{1} = JA;
varargout{2} = JB;

end
