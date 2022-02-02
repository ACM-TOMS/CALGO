function DAESAshowSM(sadata, rc, fcn, HVT, JNZ, n)

if isa(sadata, 'SAdata')
    
    exitflag = getExitflag(sadata);
    sigma = getSigma(sadata);
    fcn = getDAEfhandle(sadata);
    
    if exitflag==-3 % No visualization
        fprintf('The DAE ''%s'' is not SWP. No figure is displayed.\n', func2str(fcn));
        return;
    elseif exitflag==-2 || exitflag==-1 % Not SWP, sigma still displayed
        if nargin==1
            DAESAshowStructIllPosed(sadata);
        else DAESAshowStructIllPosed(sadata, rc);
        end
        return;
    end
    
    fcn = func2str(fcn);
    
    HVT = getHVT(sadata);
    JNZ = getJNZ(sadata);
    sigsiz = getSize(sadata);
    n = getSize(sadata);
    index = getIndex(sadata);
    DOF = getDOF(sadata);
    
elseif isa(sadata, 'numeric')
    sigma = sadata;
    sigsiz = n;
end

if length(rc)==4
    checkRC(n, rc);
    rows = rc(1):rc(2);
    cols = rc(3):rc(4);
elseif length(rc)==2
    checkRC(n, rc);
    rows = rc(1):rc(2); cols = rows;
    rc(3) = rc(1); rc(4) = rc(2);
else
    error('Wrong Size!')
end

sigblk = sigma(rows, cols);
m = length(rows); n = length(cols);
siz = max(m, n);

figureSize = getFigureSize(n,m);

titleStr = [fcn ' (' num2str(rc(1)) ':' num2str(rc(2)) ',' ...
    num2str(rc(3)) ':' num2str(rc(4)) ')'];

figure(gcf);
set(gcf, 'Name', titleStr , 'Position', figureSize);
clf, hold on;

set(gca, 'XAxisLocation','top','YDir','reverse','XLim',[0 n], 'YLim', [0 m]);
axis('equal');
axis([0 n 0 m]);

fcn = regexprep(fcn, '(_)', '\\_');
title({['\bf\fontsize{12}' upper(fcn) ': \Sigma(' ...
    num2str(rc(1)) ':' num2str(rc(2)) ',' ...
    num2str(rc(3)) ':' num2str(rc(4)) ')']...
    ,['\rm\fontsize{10} Size ', num2str(getSize(sadata)), ...
    ', structural index ', num2str(index) ...
    ', DOF ' num2str(DOF)] ...
    ,['Shaded: structural nonzeros in system J'] ...
    ,['Boxed: positions that contribute to det(J)'] ...
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
JNZblk = JNZ(rows, cols);
[i, j] = find(JNZblk);
Aij = sigblk(i + (j-1)*m);
i=i(:); j=j(:); Aij=Aij(:);
strs = num2str(Aij);
text(xpos(j), ypos(i), strs, 'BackgroundColor', 'g', 'FontSize', fontSize, ...
    'EdgeColor', 'white', 'Margin', 1);

% Find HVT in sigma
HVTmat = sparse(HVT, 1:sigsiz, ones(sigsiz,1), sigsiz, sigsiz);
HVTblk = HVTmat(rows, cols);
[i, j] = find(HVTblk);
Aij = sigblk(i + (j-1)*m);
strs = num2str(Aij);
text(xpos(j), ypos(i), strs, 'BackgroundColor', 'y', 'FontSize', fontSize, ...
    'EdgeColor', 'black',  'Margin', 1);

load 'showStructDat'; % Load visualization data
if siz < 41
    set(gca,'Position', showStructDat(n,5:8));
else
    set(gca, 'Position', [0.05,0.05,0.9,0.75]);
end

if m < 41
    set(gca,'YTick',(.5:1:m),'YTickLabel',int2str(rows'), ...
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

hold off;
box on;
end