function [strOut, maxLen] = breakLine(strIn, breaknum, delim, spaceNum, ...
    initSpaceNum, adjustSpaceNum)


if nargin<6, adjustSpaceNum=0; end
if nargin<5, initSpaceNum=0; end
if nargin<4, spaceNum=0; end

% Initialize
strOut = ''; 

% Call sub-function, considering the first line adjustment initSpaceNum
[s, strIn] = breakLineSub(strIn, breaknum-initSpaceNum, delim);

if isempty(strIn),  maxLen = min(breaknum, length(s)+5+initSpaceNum);
else maxLen = breaknum;
end

% Main loop
while 1
    % Concatenate the string
    strOut = [strOut, s];
    % If anything left, remove the previous blanks
    indx = find(strIn~=' ', 1);
    % Stopping criteria
    if isempty(strIn) || isempty(indx), break; end
    % Intercept the string that lasts
    strLast = strIn(indx:end);
    % Add adjusting blanks
    strIn = [blanks(spaceNum+adjustSpaceNum), strLast];
    % Call sub-function
    [s, strIn] = breakLineSub(strIn, breaknum, delim);
end
end

% Sub-functino for breakLine
function [strOut1, strOut2] = breakLineSub(strIn, breaknum, delim)
len = length(strIn);
if len<=breaknum
    strOut1 = strIn; strOut2 = [];
else
    indx = find(strIn == delim);
    n = find(indx<=breaknum, 1, 'last');
    strOut1 = sprintf('%s\n', strIn(1:indx(n)));
    strOut2 = strIn((indx(n)+1):end);
end
end