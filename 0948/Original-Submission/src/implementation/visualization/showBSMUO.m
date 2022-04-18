function showBSMUO(sadata, rc)

siz = getSize(sadata);
sigma = getSigma(sadata);
fcn = getDAEfhandle(sadata);
[p,q,~,cc] = DAESAgetBTF(sadata);
rr = cc(1,:); cc(1,:) = [];

if any(rc<1) || any(rc>4)
    error('Wrong block number');
end

if length(rc)==4
    % checkRC(n, rc);
    rows = rr(rc(1)):(rr(rc(2)+1)-1);
    cols = cc(rc(3)):(cc(rc(4)+1)-1);
    r = rr(rc(1):rc(2)) - rr(rc(1)) + 1;
    s = cc(rc(3):rc(4)) - cc(rc(3)) + 1;
elseif length(rc)==2
    rows = rr(rc(1)):(rr(rc(2)+1)-1); 
    cols = cc(rc(1)):(cc(rc(2)+1)-1);
    rc(3) = rc(1); rc(4) = rc(2);
    r = rr(rc(1):rc(2)) - rr(rc(1)) + 1;
    s = cc(rc(1):rc(2)) - cc(rc(1)) + 1;
else
    error('Wrong Size!')
end

m = length(rows); n = length(cols);
if m*n==0
    fprintf('\nThe block is %i-by-%i. No figure is displayed.\n\n', m,n);
    return
end

% Mark under- and over-determined parts
% [B11 B12] and [B34; B44]
JNZ = sparse(siz,siz);
JNZ(1:rr(2)-1, 1:cc(3)-1) = 1;
JNZ(rr(3):n, cc(4):n) = 1;
JNZblk = JNZ(rows, cols);

sigmapm = sigma(p, q);
sigblk = sigmapm(rows, cols);
siz = max(m, n);

% if siz==0
%     fprintf(['\nThe DAE ''%s'' does not have a well-determined part. ' ...
%         '\nNo figure is displayed.\n\n'],func2str(DAESAgetDAEfhandle(sadata)))
%    return 
% end

scrsz = get(0, 'ScreenSize'); % Get screen Size
figureSize = [scrsz(3)/2-scrsz(4)/4-scrsz(4)/200*min(n,40) ...
    scrsz(4)/2-scrsz(4)/4-scrsz(4)/200*min(m,40) ...
    scrsz(4)/2+scrsz(4)/100*min(n,40) ...
    scrsz(4)/2-30+scrsz(4)/100*min(m,40)];

fcn = func2str(fcn);
titleStr = [fcn ' - coarse (' num2str(rc(1)) ':' num2str(rc(2)) ',' ...
    num2str(rc(3)) ':' num2str(rc(4)) ')'];

figure(gcf);
set(gcf, 'Name', titleStr , 'Position', figureSize);
clf, hold on;

set(gca, 'XAxisLocation','top','YDir','reverse','XLim',[0 n], 'YLim', [0 m]);
axis('equal');
axis([0 n 0 m]);

fcn = regexprep(fcn, '(_)', '\\_');
title({['\bf\fontsize{12}' upper(fcn) ': Diagnostic blocks \Sigma(' ...
    num2str(rc(1)) ':' num2str(rc(2)) ',' ...
    num2str(rc(3)) ':' num2str(rc(4)) ')']...
    ,['\rm\fontsize{10} Structurally ill posed'] ...
    , []
    });

if siz <= 20, fontSize = 10;
elseif siz <= 30, fontSize = 8;
elseif siz <= 40, fontSize = 7;
elseif siz <= 55, fontSize = 6;
elseif siz <= 70, fontSize = 5;
else fontSize = 4;
end

xpos = (1:n)'-0.5; ypos = (1:m)'-0.5;

% Find finite entries
finsigblk = isfinite(sigblk);
[i, j] = find(finsigblk);
Aij = sigblk(i + (j-1)*m);
i=i(:); j=j(:); Aij=Aij(:);
strs = num2str(Aij(:));
text(xpos(j), ypos(i), strs, 'FontSize', fontSize);

% Find structual non-zeros
[i, j] = find(JNZblk.*finsigblk==1);
Aij = sigblk(i + (j-1)*m);
i=i(:); j=j(:); Aij=Aij(:);
strs = num2str(Aij(:));
text(xpos(j), ypos(i), strs, 'BackgroundColor', 'r', 'FontSize', fontSize);

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
% Vertical lines from (s(k)-1, 0) to (s(k)-1, m-1) for each k.
siz = size(s);
x = [s-1; s-1; NaN(siz)]; x = x(:);
%siz = siz-1;
y = [zeros(siz); repmat(m,siz); NaN(siz)]; y = y(:);
plot(x, y, ':b');

hold off;
box on;
end