function figureSize = getFigureSize(n,m)

% Decide Figure's Size, determined by size of the problem
scrsz = get(0, 'ScreenSize'); % Get screen Size
if nargin==1, m=n; end

figureSize = [scrsz(3)/2-scrsz(4)/4-scrsz(4)/200*min(n,40) ...
    scrsz(4)/2-scrsz(4)/4-scrsz(4)/200*min(m,40) ...
    scrsz(4)/2+scrsz(4)/100*min(n,40) ...
    scrsz(4)/2-30+scrsz(4)/100*min(m,40)];
end