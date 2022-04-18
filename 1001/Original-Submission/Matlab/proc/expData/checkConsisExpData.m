%% checkConsisExpData
%
% Check consistency of parameters in struct |seti| for experimentally 
% measured data. If not set, default is set.
% 
%
%% Syntax
%
%   seti = checkConsisExpData(seti)
%   seti = checkConsisExpData(seti,dispDepth)
%
%% Description
% |seti = checkConsisExpData(seti)| sets parameters in struct |seti| only
% if the field |expData| of the structure array |seti| exists and is
% |'fresnel'| or |'simonetti'|. Then some fixed variables are overwritten 
% and some are defined if they was not set.
% Otherwise the function does nothing.
%
% |seti = checkConsisExpData(seti,dispDepth)| does the same but allows to
% control the displayed messages by |dispDepth|.
%
% * If the field |expData| is set to |'fresnel'|, then the dimension
% |seti.dim| is set to 2, the incident type |seti.incType| to
% |'pointSource'| and the measurement type |seti.measType| to
% |'nearField'|.
%
% See the comments in Code for detailed information.
%
%% Example
%
%   init;
%   seti.expData = 'fresnel';
%   seti = checkConsisExpData(seti);
%
%% Input Arguments
%
% * seti    : structure array
%
% * seti.expData    : This field is not required. If it is not set, the
% structure |seti| is not changed. Possible inputs: |'fresnel'|. (Note that
% |'simonetti'| is not supported in the public version of this ode.)
%
% The following input arguments can be set by user, otherwise default values are
% set (in case of |seti.expData = 'fresnel'|):
%
% * |seti.rCD|  : Size of computational domain [-rCD,rCD)^dim
%                 (default: 0.2 m), see <setGrid.html>.
% * |seti.nCD|  : Number of discretization points for each dimension of CD
%                 (default: 256), see <setGrid.html>.
% * |seti.fresnelFreq|  : Frequency of the Fresnel data set (default: 5 GHz
% = |5*1E9|)
% * |seti.fresnelFile|  : path to the file with data from Institute Fresnel
%                         (default: 'inexpdata/fresnel_opus_1/twodielTM_8f.exp').
%
% * |seti.nuMax|    : match incident field using Hankel functions of first kind and orders
%                     $\nu = -\texttt{nuMax}, ..., -1, 0, 1, ...,
%                     \texttt{nuMax}$ (default: 7)
% * |seti.ampCalc|  : method to compute the coefficients c (1, 2 or 3),
%                     (default: 1), see <matchIncField.html>.
%
% *Optional Input Argument*
%
% * dispDepth   : Controls the displayed messages by a number (between 0
% and 5; in this file we have only two cases: 0 or higher).
%
%% Output Arguments
%
% * seti    : structure array
%
% In case of |seti.expData = 'fresnel'| the output arguments of |seti| are
%
% * |seti.dim = 2|  : The dimension of the problem is 2.
% * |seti.incType =  'pointSource'| : The type of incident field is set to point sources.
% * |seti.measType = 'nearField'|   : The measurement type is set to near field data.
%
% Additional fields of |seti| may set, if they was not set by the user (see
% Input Arguments).
%
%% See Also
%
% * <checkfield.html>
% * <setGrid.html>
% * <fresnel.html>
% * <matchIncField.html>
% * <readRAWData.html>
% * <loadData.html>
%

%% Code
function seti = checkConsisExpData(seti,varargin)

if nargin == 1
    dispDepth = 0;
else
    dispDepth = varargin{1};
end

if isfield(seti,'expData')
    %% Code: Fresnel data
    if strcmp(seti.expData, 'fresnel')
        if dispDepth >= 1
            disp(' - You decided to use Fresnel data. Several parameters will be set:')
        end

        %%
        % *Fixed Fresnel settings*
        %
        % See <fresnel.html> and <readRAWData.html>.
        if dispDepth >= 1
            disp('-- Fixed Fresnel settings --')
        end
        seti.dim = 2;
        if dispDepth >= 1
            disp('   Set dim = 2.')
        end
        seti.incType = 'pointSource';
        if dispDepth >= 1
            disp('   Set incType = pointSource.')
        end
        seti.measType = 'nearField';
        if dispDepth >= 1
            disp('   Set measType = nearField.')
        end

        %%
        % *Changeable Fresnel settings*
        %
        % See <fresnel.html>.
        if dispDepth >= 1
            disp('-- Changeable Fresnel settings --')
        end
        seti = checkfield(seti,'rCD',0.2,dispDepth); % 0.2 m
        seti = checkfield(seti,'nCD',256,dispDepth);
        seti = checkfield(seti,'fresnelFreq',5*1E9,dispDepth); % 5 GHz
        seti = checkfield(seti,'fresnelFile','inexpdata/fresnel_opus_1/twodielTM_8f.exp',dispDepth);
        
        %%
        % *Changeable Fresnel settings to match the incident field*
        %
        % See <matchIncField.html>.
        seti = checkfield(seti,'nuMax',7,dispDepth);
        seti = checkfield(seti,'ampCalc',1,dispDepth);
        
        %%
        % *Parameters set in loadData*
        %
        % using readRAWData (<readRAWData.html>) 
        % (parameters depend on used fresnel data set)
        %
        % * |seti.incNb|
        % * |seti.measNb|
        % * |seti.radSrc|
        % * |seti.radMeas|
        %

    %% Code: Simonetti data
    % (not available in public version of package)
    elseif strcmp(seti.expData,'simonetti')
        if dispDepth >= 1
            disp(' - You decided to use Simonetti"s data. Several parameters will be set:')
            disp('   *Simonetti data are experimentally implemented...*');
        end
        seti.dim = 2;
        if dispDepth >= 1
            disp('    - Set dim = 2.')
        end
        
    else
       error('seti.expData is unknown ... ?!?'); 
    end
end

end
