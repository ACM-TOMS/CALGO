function options=bvptwpset(varargin)
%bvptwpset  Create/alter BVP OPTIONS structure.
%   OPTIONS = bvptwpset('NAME1',VALUE1,'NAME2',VALUE2,...) creates an integrator
%   options structure OPTIONS in which the named properties have the
%   specified values. Any unspecified properties have default values. It is 
%   sufficient to type only the leading characters that uniquely identify the
%   property. Case is ignored for property names. 
%   
%   OPTIONS = bvptwpset(OLDOPTS,'NAME1',VALUE1,...) alters an existing options
%   structure OLDOPTS. 
%   
%   OPTIONS = bvptwpset(OLDOPTS,NEWOPTS) combines an existing options structure
%   OLDOPTS with a new options structure NEWOPTS. Any new properties overwrite 
%   corresponding old properties. 
%   
%   bvptwpset with no input arguments displays all property names and their
%   possible values. 
%   
%bvptwpset PROPERTIES
% RelTol -
%  positive scalar {1e-3} or vector.
%  Relative tolerance for the error. If RelTol is a vector
%  a different tolerance is used for all components of the solution. 
%  If vector, RelTol(i)=0 means that the
%  corresponding component is not used in the error estimation.
%
%FJacobian - Analytical partial derivatives of ODEFUN 
%          [ function_handle | matrix | cell array ]
%   For example, when solving y' = f(x,y), set this property to @FJAC if
%   DFDY = FJAC(X,Y) evaluates the Jacobian of f with respect to y.
%   In the absence of this option, the Jacobian is
%   calculated numerically with odenumjac from MATLAB.
%
%BCJacobian - Analytical partial derivatives of BCFUN 
%           [ function_handle | cell array ]
%   For example, for boundary conditions bc(ya,yb) = 0, set this property to
%   @BCJAC if [DBCDYA,DBCDYB] = BCJAC(YA,YB) evaluates the partial
%   derivatives of bc with respect to ya and to yb. 
%
%Stats - Display computational cost statistics  [ on | {off} ].
%
%MaxNumberOfMeshes - [positive integer{floor(50000/n)}]
%  Maximum number of mesh points allowed.
%
%Linear -  [on| {off}]
%  Used to indicate whether or not the problem is linear. 
%  If enabled, the problem is solved taking into account the linear behavior.
%
%Solver -'twpbvp m'| {'twpbvp_l'} (for linear problems) |'acdc' |'twpbvpc_m' | {'twpbvpc_l'} (for non linear problems)|'acdcc'
%   Use deferred correction based on MIRK (twpbvp_m) or LOBATTO methods (twpbvp_l, default
%   for linear problems), or a continuation strategy based on LOBATTO methods (acdc). The
%   implementations based on conditioning can be selected with twpbvpc_m, twpbvpc_l (default
%   for nonlinear problems) or acdcc.
%
%
%Vectorized - JVectorized ODE function  [ on | {off} ]
%   Set this property 'on' if the derivative function 
%   ODEFUN([x1 x2 ...],[y1 y2 ...]) returns [ODEFUN(x1,y1) ODEFUN(x2,y2) ...].  
%
%
%JVectorized - Vectorized Jacobian ODE function  [ on | {off} ]
%   Set this property 'on' if the derivative function 
%   FJAC([x1 x2 ...],[y1 y2 ...]) returns [FJAC(x1,y1) FJAC(x2,y2) ...].  
%  
%LambdaStart: positive scalar {0.5} 
%   Starting value for the continuation parameter, used only for acdc(c) solver.
%
%LambdaMin: positive scalar {1e-5} 
%   Final value for the continuation parameter, used only for acdc(c) solver.
%
%MaxNumberOfMeshes: positive integer {100} 
%   This is the maximum number of different meshes that could be used during computation. 
%   It is important to avoid loops.
%
%MaxNumberOfContStep: positive integer {150} 
%   This is the maximum number of continuation steps (only for acdc/acdcc solvers).
%
%Debug: on|{off}
%    Enable or disable debugging information (very verbose!).
%
%   See also bvptwpget, bvptwp, bvpMtest.  
%  
%
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

if (nargin == 0) && (nargout == 0)
  fprintf('      RelTol: [ positive scalar {1e-3} or a vector.]\n');  
  fprintf('      FJacobian: [ function_handle ]\n');
  fprintf('      BCJacobian: [ function_handle ]\n');
  fprintf('      Stats: [ on | {off} ]\n');
  fprintf('      NMax: [ nonnegative integer {floor(50000/n)} ]\n'); 
  fprintf('      MaxNumberOfMeshes: [ nonnegative integer {100} ]\n'); 
  fprintf('      MaxNumberOfContSteps: [ nonnegative integer {150)} ]\n');
  fprintf('      LambdaMin: [final value of the continuation parameter   {1d-5} ]\n'); 
  fprintf('      LambdaStart: [initial value of the continuation parameter   {0.5} ]\n'); 
  fprintf('      Linear: [on | {off})} ]\n'); 
  fprintf('      Solver : [ twpbvp_m | {twpbvpc_m} (for non linear problems) | {twpbvp_l}  (for linear problems) | twpbvpc_l | acdc | acdcc  ]\n');
  fprintf('      JVectorized: [ on | {off} ]\n');
  fprintf('      Vectorized: [ on | {off} ]\n'); 
  fprintf('      Debug : [ on | {off} ]\n'); 
  fprintf('\n');
  return;
end

 Names=[        'RelTol              '
                'FJacobian           '
                'BCJacobian          '
                'Stats               '
                'Nmax                '
                'MaxNumberOfMeshes   '
                'MaxNumberOfContSteps'
                'LambdaMin           '
                'LambdaStart         '
                'Linear              '
                'Solver              '
                'JVectorized         '
                'Vectorized          '
                'Debug               '
                ];
            
m = size(Names,1);
names = lower(Names);

% Combine all leading options structures o1, o2, ... in odeset(o1,o2,...).
options = [];
i = 1;
while i <= nargin
  arg = varargin{i};
  if ischar(arg)                         % arg is an option name
    break;
  end
  if ~isempty(arg)                      % [] is a valid options argument
    if ~isa(arg,'struct')
      error('MATLAB:bvptwpset:NoPropNameOrStruct',...
            ['Expected argument %d to be a string property name '...
                     'or an options structure\ncreated with bvptwpset.'], i);
    end
    if isempty(options)
      options = arg;
    else
      for j = 1:m
        val = arg.(deblank(Names(j,:)));
        if ~isequal(val,[])             % empty strings '' do overwrite
          options.(deblank(Names(j,:))) = val;
        end
      end
    end
  end
  i = i + 1;
end
if isempty(options)
  for j = 1:m
    options.(deblank(Names(j,:))) = [];
  end
end

% A finite state machine to parse name-value pairs.
if rem(nargin-i+1,2) ~= 0
  error('MATLAB:bvptwpset:ArgNameValueMismatch',...
        'Arguments must occur in name-value pairs.');
end
expectval = 0;                      % start expecting a name, not a value
while i <= nargin
  arg = varargin{i};    
  if ~expectval
    if ~ischar(arg)
      error('MATLAB:bvptwpset:InvalidPropName',...
        'Expected argument %d to be a string property name.', i);
    end
    lowArg = lower(arg);
    j = strmatch(lowArg,names);
    if isempty(j)                       % if no matches
      error('MATLAB:bvptwpset:InvalidPropName',...
            'Unrecognized property name ''%s''.', arg);
    elseif length(j) > 1                % if more than one match
      % Check for any exact matches (in case any names are subsets of others)
      k = strmatch(lowArg,names,'exact');
      if length(k) == 1
        j = k;
      else
        msg = sprintf('Ambiguous property name ''%s'' ', arg);
        msg = [msg '(' deblank(Names(j(1),:))];
        for k = j(2:length(j))'
          msg = [msg ', ' deblank(Names(k,:))];
        end
        error('MATLAB:bvptwpset:AmbiguousPropName', '%s).', msg);
      end
    end
    expectval = true;                      % we expect a value next    
  else
    options.(deblank(Names(j,:))) = arg;
    expectval = false;      
  end
  i = i + 1;
end

if expectval
  error('MATLAB:bvptwpset:NoValueForProp',...
        'Expected value for property ''%s''.', arg);
end

