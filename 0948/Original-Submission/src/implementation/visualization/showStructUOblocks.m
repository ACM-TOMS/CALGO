function showStructUOblocks(sadata)

exitflag = getExitflag(sadata);
fcn = DAESAgetDAEfhandle(sadata);

switch exitflag
    case 0
        error('%s is structurally well-posed (SWP). No figure is displayed.\n', ...
            func2str(fcn));
%     case -2
%         error('The DAE ''%s'' has missing equation or variable. No figure is displayed.\n', ...
%             func2str(fcn));
end


sigma = DAESAgetSigma(sadata);
fcn = DAESAgetDAEfhandle(sadata);

[p,q,~,s] = DAESAgetBTF(sadata);

%% Visualization
[m, n] = size(sigma);
siz = max(m, n);

% Find finite entries in PERMUTED sigma
sigma = sigma(p, q);
[i, j] = find(isfinite(sigma));
Aij = sigma(i + (j-1)*m);

scrsz = get(0,'ScreenSize'); % Get screen Size
figureSize = [scrsz(3)/2-scrsz(4)/4-scrsz(4)/200*min(n,40) ...
    scrsz(4)/2-scrsz(4)/4-scrsz(4)/200*min(m,40) ...
    scrsz(4)/2+scrsz(4)/100*min(n,40) ...
    scrsz(4)/2-30+scrsz(4)/100*min(m,40)];
figure(gcf);

% Use function name or file name to set figure title
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

% Structural Analysis information, on the title line above the matrix

fcn = regexprep(fcn, '(_)', '\\_');
switch exitflag
    case -1
        title({['\bf\fontsize{12}' upper(fcn) ': Diagnostic BTF']...
            ,['\rm\fontsize{10} Structurally ill posed'], ..., size: ', num2str(n)], ...
            ['Shaded: under- and over-determined']}, ...
            'HorizontalAlign', 'center');
        
    case -2
        title({['\bf\fontsize{12}' upper(fcn) ': Diagnostic BTF']...
            ,['\rm\fontsize{10} Structurally ill posed'], ..., size: ', num2str(n)], ...
            ['Missing equation/variable']}, ...
            'HorizontalAlign', 'center');
end

% Set different fontSize for problems of different sizes
if siz <= 20, fontSize = 10;
elseif siz <= 30, fontSize = 8;
elseif siz <= 40, fontSize = 7;
elseif siz <= 55, fontSize = 6;
elseif siz <= 70, fontSize = 5;
else fontSize = 4;
end

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
    set(gca,'YTick',(.5:1:n),'YTickLabel',int2str(p'), ...
        'TickLength',[0 0],'FontSize',fontSize);
    ylabel('Indices of Equations','FontSize',10);
else
    set(gca, 'YTick', []);
end

if n < 41
    set(gca,'XTick',(.5:1:n),'XTickLabel',int2str(q'), ...
        'TickLength',[0 0],'FontSize',fontSize);
    xlabel('Indices of Variables','FontSize',10);
else
    set(gca, 'XTick', [])
end

%% Plot lines to show UNDER-, WELL-, and OVER-determined blocks.
siz = [1 5];
r = s(1,:); s(1,:) = [];
% Horizontal lines lines from (0, r(k)-1) to (n, r(k)-1) for each k.
x = [r-1; r-1; NaN(siz)]; x = x(:);
y = [zeros(siz); repmat(n,siz); NaN(siz)]; y = y(:);
plot(y, x, ':b');

% Vertical lines from (s(k)-1, 0) to (s(k)-1, n-1) for each k.
x = [s-1; s-1; NaN(siz)]; x = x(:);
y = [zeros(siz); repmat(m,siz); NaN(siz)]; y = y(:);
plot(x, y, ':b');

%% Color
% Define color
darkPink = [255 64 64]/256;
lightPink = [255 193 193]/256;

% Under-determined part
partSigma = sigma(1:r(2)-1,1:s(3)-1);   
[i, j] = find(isfinite(partSigma));
Aij = sigma(i + (j-1)*m);
strs = num2str(Aij(:));
text(xpos(j), ypos(i), strs, 'BackgroundColor', darkPink, 'FontSize', fontSize);

line([0 0],[0 r(2)-1],'LineWidth',3,'Color',darkPink)
line([0 s(3)-1],[0 0],'LineWidth',3,'Color',darkPink)
line([0 s(3)-1 s(3)-1],[r(2)-1 r(2)-1 0],'LineWidth',2,'Color',darkPink)

% Over-determined part
partSigma = sigma(r(3):n,s(4):n);       
[i, j] = find(isfinite(partSigma));
i = i+r(3)-1;
j = j+s(4)-1;
Aij = sigma(i + (j-1)*m);
strs = num2str(Aij(:));
text(xpos(j), ypos(i), strs, 'BackgroundColor', darkPink, 'FontSize', fontSize);

line([n n],[r(3)-1 n],'LineWidth',3,'Color',darkPink)
line([s(4)-1 n],[n n],'LineWidth',3,'Color',darkPink)
line([s(4)-1 s(4)-1 n],[n r(3)-1 r(3)-1],'LineWidth',2,'Color',darkPink)

[meqn, mvar] = DAESAgetMissing(sadata);
nmeqn = length(meqn);
nmvar = length(mvar);
text(xpos(repmat((1:nmvar)',n,1)), ypos(reshape(repmat(1:n,nmvar,1),n*nmvar,1)), ...
    repmat(' ', nmvar*n, 1), 'BackgroundColor', lightPink)

text(xpos(repmat((1:n)',nmeqn,1)), ypos(n+1-reshape(repmat(1:nmeqn,n,1),n*nmeqn,1)), ...
    repmat(' ', nmeqn*n, 1), 'BackgroundColor', lightPink)

hold off;
box on;
end
