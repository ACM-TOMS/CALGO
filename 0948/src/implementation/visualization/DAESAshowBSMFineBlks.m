function DAESAshowBSMFineBlks(sadata, rc)

sigma = getSigma(sadata);
fcn = getDAEfhandle(sadata);
JNZ = getJNZ(sadata);
[p, q, ~, rF] = DAESAgetBTF(sadata);
index = getIndex(sadata);
DOF = getDOF(sadata);

numBlks = length(rF)-1;
if any(rc<1) || any(rc>numBlks)
    error('Wrong block number');
end

if length(rc)==4
    % checkRC(n, rc);
    rows = rF(rc(1)):(rF(rc(2)+1)-1);
    cols = rF(rc(3)):(rF(rc(4)+1)-1);
    r = rF(rc(1):rc(2)) - rF(rc(1)) + 1;
    s = rF(rc(3):rc(4)) - rF(rc(3)) + 1;
elseif length(rc)==2
    % checkRC(n, rc);
    rows = rF(rc(1)):(rF(rc(2)+1)-1); cols = rows;
    rc(3) = rc(1); rc(4) = rc(2);
    r = rF(rc(1):rc(2)) - rF(rc(1)) + 1;
    s = r;
else
    error('Wrong Size!')
end

sigmapm = sigma(p, q);
sigblk = sigmapm(rows, cols);
m = length(rows); n = length(cols);
siz = max(m, n);

figureSize = getFigureSize(n,m);

fcn = func2str(fcn);
titleStr = [fcn ' - fine (' num2str(rc(1)) ':' num2str(rc(2)) ',' ...
    num2str(rc(3)) ':' num2str(rc(4)) ')'];

figure(gcf);
set(gcf, 'Name', titleStr , 'Position', figureSize);
clf, hold on;

set(gca, 'XAxisLocation','top','YDir','reverse','XLim',[0 n], 'YLim', [0 m]);
axis('equal');
axis([0 n 0 m]);

fcn = regexprep(fcn, '(_)', '\\_');
title({['\bf\fontsize{12}' upper(fcn) ': Fine Blocks \Sigma(' num2str(rc(1)) ':' num2str(rc(2)) ',' ...
    num2str(rc(3)) ':' num2str(rc(4)) ')']...
    ,['\rm\fontsize{10} Size ', num2str(getSize(sadata)), ...
    ', structural index ', num2str(index) ...
    ', DOF ' num2str(DOF)] ...
    ,'Shaded: structural nonzeros in system J' ...
    ,'Boxed: positions that contribute to det(J)' ...
    , []
    });

fontSize = getFontSize(n);

xpos = (1:n)'-0.5; ypos = (1:m)'-0.5;

% Find finite entries
[i, j] = find(isfinite(sigblk));
Aij = sigblk(i + (j-1)*m);
i=i(:); j=j(:); Aij=Aij(:);
strs = num2str(Aij);
text(xpos(j), ypos(i), strs, 'FontSize', fontSize);

% Find structual non-zeros
JNZ = JNZ(p, q);
JNZblk = JNZ(rows, cols);
[i, j] = find(JNZblk==1);
Aij = sigblk(i + (j-1)*m);
i=i(:); j=j(:); Aij=Aij(:);
strs = num2str(Aij);
text(xpos(j), ypos(i), strs, 'BackgroundColor', 'y', 'FontSize', fontSize, ...
    'EdgeColor', 'black', 'Margin',1);

% Find HVT in sigma
[i, j] = find(JNZblk==-1);
Aij = sigblk(i + (j-1)*m);
strs = num2str(Aij);
text(xpos(j), ypos(i), strs, 'BackgroundColor', 'g', 'FontSize', fontSize, ...
    'Margin', 1);

load 'showStructDat'; % Load visualization data
if siz < 41
    set(gca,'Position', showStructDat(n,5:8));
else
    set(gca, 'Position', [0.05,0.05,0.9,0.75]);
end

if m < 41
    set(gca,'YTick',(.5:1:m),'YTickLabel',int2str((p(rows))'), ...
        'TickLength',[0 0],'FontSize',fontSize);
    ylabel('Indices of Equations','FontSize',10);
else
    set(gca, 'YTick', []);
end

if n < 41
    set(gca,'XTick',(.5:1:n),'XTickLabel',int2str((q(cols))'), ...
        'TickLength',[0 0],'FontSize',fontSize);
    xlabel('Indices of Variables','FontSize',10);
else
    set(gca, 'XTick', [])
end

%% Plot lines to show blocks.
% Note vectors x, y are re-used.
% HorizontalAlignmentontal lines from (0, r(k)-1) to (n, r(k)-1) for each k.
siz = size(r);
x = [r-1; r-1; NaN(siz)]; x = x(:);
y = [zeros(siz); repmat(n,siz); NaN(siz)]; y = y(:);
plot(y, x, ':b');
% Vertical lines from (s(k)-1, 0) to (s(k)-1, n-1) for each k.
siz = size(s);
x = [s-1; s-1; NaN(siz)]; x = x(:);
%siz = siz-1;
y = [zeros(siz); repmat(m,siz); NaN(siz)]; y = y(:);
plot(x, y, ':b');

hold off;
box on;
end