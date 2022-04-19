function DAESAshowStructIllPosed(sadata, rc)

fcn = getDAEfhandle(sadata);
sigma = getSigma(sadata);

%% Ill posed Positions
[m, n] = size(sigma);
illPosedPosn = zeros(n);
sigmaFinite = isfinite(sigma);
% Equations not evaluated
row = any(sigmaFinite, 2)==0;
illPosedPosn(row, :) = 1;
% Variables not used
col = any(sigmaFinite, 1)==0;
illPosedPosn(:, col) = 1;

%% Show Part of sigma
if nargin==2
    checkRC(n, rc);
    if length(rc)==4
        rows = rc(1):rc(2);
        cols = rc(3):rc(4);
        sigma = sigma(rows, cols);
        illPosedPosn = illPosedPosn(rows, cols);
    elseif length(rc)==2
        rc(3)=rc(1) ; rc(4)=rc(2);
        rows = rc(1):rc(2); cols = rows;
        sigma = sigma(rows, cols);
        illPosedPosn = illPosedPosn(rows, cols);
    else
        error('Wrong length!');
    end
else
    rows = 1:m; cols = 1:n;
end

%% Visualization
[m, n] = size(sigma);
siz = max(m, n);

% Find finite entries in sigma
[i, j] = find(isfinite(sigma));
Aij = sigma(i + (j-1) * m);

figureSize = getFigureSize(n,m);

% Use function name or file name to set figure title
figure(gcf);
if ~isempty(fcn)
    if isa(fcn, 'function_handle'),  fcn = func2str(fcn);  end
    set(gcf, 'Name' ,fcn, 'NumberTitle', 'off', 'Position', figureSize);
else
    set(gcf, 'Position', figureSize);
end
clf; hold on;

% Set up axes.
axis equal
set(gca, 'XAxisLocation','top','YDir','reverse','XLim',[0 n], 'YLim', [0 m]);
% axis('square');

fcn = regexprep(fcn, '(_)', '\\_');
% Structural Analysis information, on the title line above the matrix
if nargin==1
    title({['\bf\fontsize{12}' upper(fcn)], '\rm\fontsize{10} Structurally ill posed' , []});
elseif nargin==2
    title({['\bf\fontsize{12}' upper(fcn) ': \Sigma(' num2str(rc(1)) ':' num2str(rc(2)) ',' ...
        num2str(rc(3)) ':' num2str(rc(4)) ')']...
        '\rm\fontsize{10}Structurally ill posed' , []});
end

% Set different fontSize for problems of different sizes
fontSize = getFontSize(n);

% Form arrays of cell-centre coordinates:
xpos = (1:n)'-0.5; ypos = (1:m)'-0.5;

% Put text value of each matrix entry in centre of cell on graph.
strs = num2str(Aij);
text(xpos(j), ypos(i), strs,'FontSize',fontSize);

load 'showStructDat'; % Load visualization data
if siz < 41
    set(gca,'Position', showStructDat(n,5:8));
else
    set(gca, 'Position', [0.05,0.05,0.9,0.75]);
end

if m < 41
    set(gca,'YTick',(.5:1:n),'YTickLabel',int2str(rows'), ...
        'TickLength',[0 0],'FontSize',fontSize);
    ylabel('Indices of Equations','FontSize',10);
else
    set(gca, 'YTick', []);
end

if n < 41
    set(gca,'XTick',(.5:1:n),'XTickLabel',int2str(cols'), ...
        'TickLength',[0 0],'FontSize',fontSize);
    xlabel('Indices of Variables','FontSize',10);
else
    set(gca, 'XTick', [])
end

[illI, illJ] = find(illPosedPosn);
text(xpos(illJ), ypos(illI), ' ', 'BackgroundColor', 'r', 'FontSize', fontSize);

hold off;
box on;
end
