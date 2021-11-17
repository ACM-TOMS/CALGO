function ascertain(x__,s__)
%ASCERTAIN  Stop if condition does not hold
%  ASCERTAIN(X) stops with an error message if X is 0 (false) but
%  continues quietly if X holds. ASCERTAIN(X,S) uses error message S.
%  ASCERTAIN('EXPRESSION') evaluates EXPRESSION in the calling workspace 
%  and quits, using the the EXPRESSION itself as an error message, if
%  it evaluates to 0 (false).
if ischar(x__)
	ok__ = evalin('caller', x__);
	if ~ok__,
    disp ' '
    error(['Assertion ' x__ ' failed']); 
  end
else
	if nargin == 1, s__ = 'Assertion failed'; end
	if ~x__, 
    disp ' '
    error(s__); 
  end
end