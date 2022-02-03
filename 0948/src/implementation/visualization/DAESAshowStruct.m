function DAESAshowStruct(sadata, rc)
% function showStruct(sadata)
% function showStruct(sadata, rc)
% function showStruct(sigma)
% function showStruct(filename)
%
% Input:
%        sadata: an SAdata class structure
%    (or) sigma: signature matrix
% (or) filename: file that contains a signature matrix
% rc (optional): a row vector specifying a sub-block of the signature matrix.
%
% If rc is a vector of length 4, submatrix [rc(1):rc(2), rc(3):rc(4)] is
% displayed.
% If rc is a vector of length 2, the diagonal submatrix
% [rc(1):rc(2), rc(1):rc(2)] is displayed.
%
% Output:
%       A figure that visualizes the structure of signature matrix in
%       original order, displaying problem name, structural index, DOF,
%       offsets, HVT and structurally non-zeros in Jacobian.

if isa(sadata, 'SAdata')
    
    exitflag = getExitflag(sadata);
    fcn = getDAEfhandle(sadata);
    sigma = getSigma(sadata);
    
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
    
    n = getSize(sadata);
    DOF = getDOF(sadata);
    HVT = getHVT(sadata);
    [c, d] = getOffsets(sadata);
    JNZ = getJNZ(sadata);
    
    if nargin == 2, showSM(sadata, rc); return; end
    
    %% Disabled for now
else
    if isa(sadata,'numeric')
        fcn = inputname(1);
        sigma = sadata;
        
    elseif isa(sadata,'char')
        [pathstr, fcn, ~] = fileparts(sadata);
        oldFolder = cd (pathstr);
        sigma = eval(fcn);
        cd(oldFolder);
    else
        error('Input should either be a sadata object, a matrix or a file containing the Sigma matrix.');
    end
    [c, d, ~, ~, HVT, DOF, ~, ~, ~, ~, ~, exitflag] = lapdm(sigma);
    if exitflag==-2
        if nargin==1
            showStructIllPosed(fcn, sigma);
        else showStructIllPosed(fcn, sigma, rc);
        end
        return;
    end
    n = size(sigma, 1);
    c = c';
    [i, j] = find(repmat(d,[n 1]) - repmat(c', [1 n]) == sigma);
    JNZ = sparse(i, j, ones(length(i), 1), n, n);
    
    if nargin == 2
        showBlock(sigma, rc, fcn, HVT, JNZ, n)
        return;
    end
end

% Find finite entries in sigma
[i, j] = find(isfinite(sigma));
Aij = sigma(i + (j-1) * n);
m = n; % Used in the plotting. In future we may allow m~=n.

% Compute structural index by the definition in John D. Pryce's BIT paper:
index = max(c) + (~isempty(find(d == 0,1)));

%% Produce a plot showing sigma

% Decide Figure's Size, determined by size of the problem
figureSize = getFigureSize(n);

% Use function name or file name to set figure title
figure(gcf);

if ~isempty(fcn)
    if isa(fcn, 'function_handle'),  fcn = func2str(fcn);  end
    set(gcf, 'Name', fcn , 'Position', figureSize);
else
    set(gcf, 'Position', figureSize);
end
clf, hold on;

% Set up axes.
axis equal
set(gca, 'XAxisLocation','top','YDir','reverse','XLim',[0 n], 'YLim', [0 m]);

fcn = regexprep(fcn, '(_)', '\\_');
title({['\bf\fontsize{12} ' upper(fcn)]...
    ,['\rm\fontsize{10} Size ', num2str(n), ...
    ', structural Index ', num2str(index) ...
    ', DOF ' num2str(DOF)] ...
    ,'Shaded: structural nonzeros in system Jacobian J' ...
    ,'Boxed: HVT' ...
    , []
    }, 'HorizontalAlign', 'center');

% Set different fontSize for problems of different sizes
fontSize = getFontSize(n);

% Form arrays of cell-centre coordinates:
xpos = (1:n)'-0.5; ypos = (1:m)'-0.5;

% Put text value of each matrix entry in centre of cell on graph.
strs = num2str(Aij);
text(xpos(j), ypos(i), strs,'FontSize',fontSize);

[iNZ, jNZ] = find(JNZ);
Aij = iNZ + (jNZ-1)*n;
strsNZ = num2str(sigma(Aij));
text(xpos(jNZ), ypos(iNZ), strsNZ, 'BackgroundColor', 'g', 'FontSize', fontSize, ...
    'EdgeColor', 'white', 'Margin', 1);

% Obtain positions and entries in HVT
strs = num2str(diag(sigma(HVT, :)));

% Over-write the plot, to highlight the set where d(j)-c(i)=A(i,j):
text(xpos(1:n), ypos(HVT), strs, 'BackgroundColor', 'y', 'FontSize', fontSize, ...
    'EdgeColor', 'black',  'Margin', 1);

if n < 41
    load 'showStructDat'; % Load visualization data
    % Set up plot size
    set(gca,'Position', showStructDat(n,5:8));
    
    set(gca,'XTick',(.5:1:n),'XTickLabel',int2str((1:n)'), ...
        'TickLength',[0 0],'FontSize',fontSize);
    set(gca,'YTick',(.5:1:n),'YTickLabel',int2str((1:n)'), ...
        'TickLength',[0 0],'FontSize',fontSize);
    
    xlabel('Indices of Variables','FontSize',10);
    ylabel('Indices of Equations','FontSize',10);
    
    % Offsets for offsets
    X = showStructDat(n, 1);
    Z = showStructDat(n, 3);
    
    % Write & label the offsets c,d on the picture, beyond border of array
    text(repmat(n + X, size(xpos)), xpos, num2str(c(:)), 'FontSize', fontSize);
    text('Interpreter','latex','String','$$\rm{c_i}$$', 'Position',...
        [n+X,0], 'VerticalAlignment', 'Bottom', 'FontSize', 13);
    
    text(ypos,repmat(m + Z, size(ypos)), num2str(d(:)), 'FontSize', fontSize, 'HorizontalAlignment', 'center');
    text('Interpreter','latex','String','$$\rm{d_j}$$', 'Position',...
        [0,m+Z], 'HorizontalAlignment','right', 'FontSize', 13);
    
else
    set(gca, 'Position', [0.05,0.05,0.9,0.75]);
    set(gca, 'YTick', [], 'XTick', []) % Hide axis for very large problems
end
box on;
hold off;
end