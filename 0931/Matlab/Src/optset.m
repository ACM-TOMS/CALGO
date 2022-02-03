function options=optset(varargin)
%   OPTIONS = OPTSET('PARAM1',VALUE1,'PARAM2',VALUE2,...) creates an
%     options structure OPTIONS in which the named parameters have
%     the specified values.  Any unspecified parameters are set to 
%     the default value for that parameter when OPTIONS is passed 
%     to the MULTIPLICITY function. It is sufficient to type only the 
%     leading characters that uniquely identify the parameter. 
%     Case is ignored for parameter names.
%     NOTE: For values that are strings, the complete string is required.
%  
%
% OPTSET PARAMETERS for MULTIPLICITY
%  Display - Level of display [ 0 | 1 | 2]
%  Tol - Termination tolerance,Any value below is considered zero 
%  EqnType - The type of equation [Poly | Nonlinear]
%  MaxInt - The max multiplicity
%
% Examples:
%
%       To create an options structure with Tol equal to 1e-10
%         options = optset('TolFun',1e-10);
%       To change the EqnType value of options to 'Nonlinear'
%         options = optset(options,'EqnType','Nonlinear');
 
    options=defaultstructure();
    for i=1:2:length(varargin)
        switch varargin{i}
            case 'EqnType'
                if (i+1<=length(varargin))
                    if (strcmp(varargin{i+1},'Poly') || ...
                            strcmp(varargin{i+1},'Nonlinear'))
                        options.EqnType=varargin{i+1};
                    else
                        disp([varargin{i+1},...
                            'is not a value for EqnType',...
                            ',MulLab set EqnType as "Poly" as defalut!'])
                    end
                else
                    error(['Please input the type of equation computed,',...
                        'Poly or Nonlinear!']);
                end
            case 'Display'
                if (i+1<=length(varargin))
                    if (varargin{i+1}==0 || varargin{i+1}==1 || varargin{i+1}==2)
                        options.Display=varargin{i+1};
                    else
                        disp(['MulLab just accept 0 ,1 or 2 as Display,',...
                            'right now Display=0 as defalut!'])
                    end
                else
                    error('Please input the Display, 0 ,1 or 2!');
                end
            case 'Tol'
                if (i+1<=length(varargin))
                    options.Tol=varargin{i+1};
                else
                    error('Please input the Tol!');
                end
            case 'MaxInt'
                if (i+1<=length(varargin))
                    if (varargin{i+1}>0)
                          options.MaxInt=varargin{i+1};
                     else
                          error('MaxInt should be a postive number!');
                     end
                else
                    error('Please input the max multiplicity!');
                end
        end
    end
end

function options=defaultstructure()
    options.EqnType='Nonlinear';
    options.Display=0;
    options.Tol=1e-8;
    options.MaxInt=1000;
end