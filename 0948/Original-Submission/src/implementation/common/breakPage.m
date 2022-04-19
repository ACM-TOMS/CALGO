function strOut = breakPage(strIn, numLine)

if nargin<2, numLine = 30; end

% Initialize
i = 1;

% Main loop
while ~isempty(strIn)
    tf = find(isstrprop(strIn, 'cntrl'));
    if length(tf)>numLine
        indx = tf(numLine);
        strOut{i} = strIn(1:indx-1);
        i = i+1;
        strIn = strIn(indx+1:end);
    else
        strOut{i} = strIn;
        break;
    end
end
end