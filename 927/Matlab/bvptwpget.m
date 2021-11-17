function value=bvptwpget(options,key,default)
%bvptwpget  Get BVP OPTIONS parameters.
%   VAL = bvptwpget(OPTIONS,KEY) extracts the value of the named property
%   from integrator options structure OPTIONS, returning an empty matrix if
%   the property value is not specified in OPTIONS. It is sufficient to type
%   only the leading characters that uniquely identify the property. Case is
%   ignored for property names. [] is a valid OPTIONS argument. 
%   
%   VAL = bvptwpget(OPTIONS,KEY,DEFAULT) extracts the named property as
%   above, but returns VAL = DEFAULT if the named property is not specified
%   in OPTIONS. For example 
%   
%       val = bvptwpget(opts,'RelTol',1e-4);
%   
%   returns val = 1e-4 if the RelTol property is not specified in opts.
%   
%   See also bvptwpset, bvptwp,  bvpMtest .
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

if nargin < 2
  error('MATLAB:bvptwpget:NotEnoughInputs', 'Not enough input arguments.');
end
if nargin < 3
  default = [];
end

if ~isempty(options) && ~isa(options,'struct')
  error('MATLAB:bvptwpget:OptsNotStruct',...
        'First argument must be an options structure created with TWPBVPSET.');
end

if isempty(options)
  value = default;
  return;
end
Names = [
    'RelTol              '
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
names = lower(Names);

lowName = lower(key);
j = strmatch(lowName,names);
if isempty(j)               % if no matches
  error('MATLAB:bvptwpget:InvalidPropName',...
        ['Unrecognized property name ''%s''.  ' ...
         'See TWPBVPSET for possibilities.'], key);
elseif length(j) > 1            % if more than one match
  % Check for any exact matches (in case any names are subsets of others)
  k = strmatch(lowName,names,'exact');
  if length(k) == 1
    j = k;
  else
    msg = sprintf('Ambiguous property name ''%s'' ', key);
    msg = [msg '(' deblank(Names(j(1),:))];
    for k = j(2:length(j))'
      msg = [msg ', ' deblank(Names(k,:))];
    end
    error('MATLAB:bvptwpget:AmbiguousPropName', '%s).', msg);
  end
end

if any(strcmp(fieldnames(options),deblank(Names(j,:))))
  value = options.(deblank(Names(j,:)));
  if isempty(value)
    value = default;
  end
else
  value = default;
end


