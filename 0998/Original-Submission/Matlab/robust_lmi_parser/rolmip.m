function [varargout] = rolmip(varargin)

% This function is the core of the ROLMIP, which stores all the
% variables used. The first argument of the function is the action 
% to be performed, which may be one of the following:
% - 'getvar': Gets a specified variable, or all the variables if none
%   is specified.
%   Examples:
%     vecpoly = rolmip('getvar');
%     poly_A = rolmip('getvar','A');
%
% - 'setvar': Sets a variable.
%   Example:
%     rolmip('setvar',poly);
%
% - 'clearvar': Clear all the variables
%   Example:
%     rolmip('clearvar');
%
% - 'setobj': Sets the objective function to be considered, which is
%   an SDPVAR object
%   Example:
%     rolmip('setobj',obj);
%
% - 'getobj': Gets the objective function 
%   Example:
%     obj = rolmip('getobj');
%
% - 'setopt': Sets the options for the parameter 'sdpsettings'
%   Example:
%     rolmip('setopt','''verbose'',0');
%
% - 'getopt': Gets the options for the parameter 'sdpsettings'
%   Example:
%     opt = rolmip('getopt');
%
% - 'setauxfile': Sets the fid for the auxiliary file that contains
%   the definitions of the auxiliary variables
%   Example:
%     rolmip('setauxfile',fid);
%
% - 'getauxfile': Gets the fid of the auxiliary file that contains
%   the definitions of the auxiliary variables
%   Example:
%     fid = rolmip('getauxfile');


persistent vecpoly obj opt fid

command = varargin{1};

if (strcmp(command,'setvar'))
    poly = varargin{2};
    found = false;
    cont = 1;
    while ((cont <= length(vecpoly)) && (~found))
        if (strcmp(vecpoly{cont}.label,poly.label))
            vecpoly{cont} = poly; %replace the variable by the new one
            found = true;
        else
            cont = cont + 1;
        end
    end
    if (~found)
        vecpoly{length(vecpoly)+1} = poly;
    end
elseif (strcmp(command,'getvar'))
    if (nargin == 1)
        varargout{1} = vecpoly;
    elseif (nargin == 2)
        label = varargin{2};
        found = false;
        cont = 1;
        while ((cont <= length(vecpoly)) && (~found))
            if (strcmp(vecpoly{cont}.label,label))
                found = true;
            else
                cont = cont + 1;
            end
        end
        if (found)
            varargout{1} = vecpoly{cont};
        else
            varargout{1} = [];
        end
    end
elseif (strcmp(command,'clearvar'))
    clear vecpoly;
elseif (strcmp(command,'setobj'))
    obj = varargin{2};
elseif (strcmp(command,'getobj'))
    varargout{1} = obj;
elseif (strcmp(command,'setopt'))
    opt = varargin{2};
elseif (strcmp(command,'getopt'))
    varargout{1} = opt;
elseif (strcmp(command,'setauxfile'))
    fid = varargin{2};
elseif (strcmp(command,'getauxfile'))
    varargout{1} = fid;
else
    error(['Command ', varargin{1},' not recognized']);
end

return
