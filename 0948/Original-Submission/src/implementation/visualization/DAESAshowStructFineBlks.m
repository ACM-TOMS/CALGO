function DAESAshowStructFineBlks(sadata, varargin)
% function showStructFineBlks(sadata)
% function showStructFineBlks(sadata, rc)
%
% Input:
%        sadata: an SAdata class structure
% rc (optional): a row vector specifying a sub-block of the signature matrix.
%
%   If rc is a vector of length 4, submatrix [rc(1):rc(2), rc(3):rc(4)] is
%   displayed.
%   If rc is a vector of length 2, the diagonal submatrix
%   [rc(1):rc(2), rc(1):rc(2)] is displayed.
%
%   If rc is an integer (1>=rc<=b, where b is number of blocks), diagonal
%   block with number rc is displayed. The numbering is along the main
%   diagonal blocks, starting from the top.
%
%  Output:
%       A figure that visualizes the structure of signature matrix permuted in
%       block upper triangular form, displaying problem name, structural index,
%       DOF, global, local offsets, structurally zeros in Jacobian

%% Extract information in sadata
if isa(sadata, 'SAdata')
    
    exitflag = getExitflag(sadata);
    fcn = getDAEfhandle(sadata);
    
    if exitflag<0 % No visualization
        fprintf('%s is structurally ill posed. Fine blocks cannot be displayed.\n'...
            , func2str(fcn));
        return;
    end
    
    if nargin == 2
        error('Impossible to have two arguments only!');
    elseif nargin == 3
        if varargin{1}==1
            DAESAshowSMFineBlks(sadata, varargin{2}); return;
        elseif varargin{1}==2
            DAESAshowBSMFineBlks(sadata, varargin{2}); return;
        end
    end
    
    n = getSize(sadata);
    sigma = getSigma(sadata);
    DOF = getDOF(sadata);
    [c, d, clocal, dlocal] = getOffsets(sadata);
    [p, q, ~, r] = DAESAgetBTF(sadata);
    JNZ = getJNZ(sadata);
else
    error('Input should be a sadata object!');
end

% Compute structural index by the definition in J.D. Pryce's BIT paper:
index = max(c)+(~isempty(find(d==0,1)));

%% Produce a plot showing the DMPERM blocks:

% Decide Figure's Size, determined by size of the problem
figureSize = getFigureSize(n);
fcn = func2str(fcn);

figure(gcf);
set(gcf, 'Name', [fcn ' (fine)'], 'Position', figureSize);
clf, hold on;

% Set up axes.
set(gca, 'XAxisLocation','top','YDir','reverse','XLim',[0 n], 'YLim', [0 n]);
axis('square');

fcn = regexprep(fcn, '(_)', '\\_');
title({['\bf\fontsize{12}' upper(fcn) ': Fine BTF'], ...
   ... ['\rm\fontsize{10}Permuted to block-triangular, irreducible diagonal blocks'], ...
    ['\rm\fontsize{10} Size ', num2str(n), ...
    ', structural index ', num2str(index) ...
    ', DOF ' num2str(DOF)], ...
    'Shaded: structural nonzeros in system Jacobian J' ...
    'Boxed: positions that contribute to det(J)' ...
    },'HorizontalAlign', 'center');

% Set different fontSize for problems of different sizes
fontSize = getFontSize(n);

m = n;
% JNZ stores structurally non-zeros
[iNZ, jNZ] = find(JNZ==1);
[iZ, jZ] = find(JNZ==-1);

% Find entries for finite entries
[i, j] = find(isfinite(sigma));
Aij = sigma(i+(j-1)*n);
% Find inverse permutations:
p_inv = invperm(p);
q_inv = invperm(q);

% Permute A, c and d into the block triangular form.
i = p_inv(i); j = q_inv(j); % Aij stays same

% Entries that have equality and contribute to Jacobian
iNZ = p_inv(iNZ); jNZ = q_inv(jNZ);
strsNZ = num2str(sigma(JNZ==1));

% Entries that have equality but don't contribute to Jacobian
iZ = p_inv(iZ); jZ = q_inv(jZ);
strsZ = num2str(sigma(JNZ==-1));

% Permute offsets
c = c(p); d = d(q);
clocal = clocal(p); dlocal = dlocal(q);

% Form arrays of cell-centre coordinates:
xpos = (1:n)'-0.5; ypos = (1:m)'-0.5;

% Put text value of each matrix entry in centre of cell on graph.
strs = num2str(Aij);
text(xpos(j), ypos(i), strs, 'FontSize', fontSize);

% Highlight in green the set where structually non-zero:
% The union of all HVTs in the permuted matrix are the positions where
% equality hold.
text(xpos(jNZ), ypos(iNZ), strsNZ, 'BackgroundColor', 'y', 'FontSize', fontSize, ...
    'EdgeColor', 'black', 'Margin',1);

% Over-write the plot again, hightlight in yellow the set where structually zero:
text(xpos(jZ), ypos(iZ), strsZ, 'BackgroundColor', 'g', 'FontSize', fontSize, ...
    'Margin',1);

%% Display Offsets
if n < 41
    % Load visualization data
    load 'showStructDat';
    % Set up plot size
    set(gca, 'Position', showStructDat(n, 5:8));
    % Set up tick labels and use XLABEL, YLABEL to describe them
    set(gca, 'XTick', (.5:1:n), 'XTickLabel', int2str(q(:)), ...
        'TickLength', [0 0], 'FontSize', fontSize);
    set(gca, 'YTick', (.5:1:n), 'YTickLabel', int2str(p(:)), ...
        'TickLength', [0 0], 'FontSize', fontSize);
    
    % Visualize according to problem size, offsets for offsets
    X = showStructDat(n, 1); Y = showStructDat(n, 2);
    Z = showStructDat(n, 3); W = showStructDat(n, 4);
    
    % Write & label the offsets c,d on the picture, beyond border of array
    text(repmat(n + X, size(xpos)), xpos, num2str(clocal(:)), 'FontSize', fontSize);
    text('Interpreter','latex','String','$$\rm\widehat{c}_i$$', 'Position',...
        [n+X,0], 'VerticalAlignment', 'Bottom', 'FontSize', 13);
    
    text(repmat(n + Y, size(xpos)), xpos, num2str(c(:)), 'FontSize', fontSize);    
    text('Interpreter','latex','String','$$\rm{c_i}$$', 'Position',...
        [n+Y,0], 'VerticalAlignment', 'Bottom', 'FontSize', 13);
    
    text(ypos,repmat(m + Z, size(ypos)),num2str(dlocal(:)), ...
        'FontSize', fontSize, 'HorizontalAlignment','center');
    text('Interpreter','latex','String','$$\rm\widehat{d}_j$$', 'Position',...
        [0,m+Z], 'HorizontalAlignment','right', 'FontSize', 13);
    
    text(ypos, repmat(m + W, size(ypos)), num2str(d(:)), ...
        'FontSize', fontSize, 'HorizontalAlignment', 'center');
    text('Interpreter','latex','String','$$\rm{d_j}$$', 'Position',...
        [0,m+W], 'HorizontalAlignment','right', 'FontSize', 13);
    
    xlabel('Indices of Variables', 'FontSize', 10);
    ylabel('Indices of Equations', 'FontSize', 10);
    
else
    set(gca, 'Position', [0.05,0.05,0.9,0.75]);
    set(gca, 'YTick', [], 'XTick', []) % Hide axis
end

%% Plot lines to show blocks.
% Note vectors x, y are re-used.
% HorizontalAlignmentontal lines from (0, r(k)-1) to (n, r(k)-1) for each k.
siz = size(r);
x = [r-1; r-1; NaN(siz)]; x = x(:);
y = [zeros(siz); repmat(n,siz); NaN(siz)]; y = y(:);
plot(y, x, ':b');
% Vertical lines from (s(k)-1, 0) to (s(k)-1, n-1) for each k.
siz = size(r);
x = [r-1; r-1; NaN(siz)]; x = x(:);
%siz = siz-1;
y = [zeros(siz); repmat(n,siz); NaN(siz)]; y = y(:);
plot(x, y, ':b');

hold off;
box on;
end