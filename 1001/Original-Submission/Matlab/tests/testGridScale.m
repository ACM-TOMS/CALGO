%% testGridScale
% Test of functions gridUp and gridDown to scale a grid up and down.
%
%% Syntax
%
%   testGridScale
%
%% Description
% |testGridScale| tests the grid scaling (up and down) by comparing
% scaled matrices in 2D and 3D.
%
% This test should run in environment runtests, see <runtests.html> via
%
%   runtests('g')
%
%% Input Arguments
%
% dim   :   dimension (2 or 3) of the matrix
%
%% Output Arguments
%
% In case of 2D some figures, see <runtests.html>.
%
%
%% See More
% * <runtests.html>
% * <gridUp.html>
% * <gridDown.html>
%
%% Code
%
% *Input*
%
% dim = 2 or 3

disp('test gridScale: gridUp and gridDown')

% input
if ~exist('dim','var')
    dim = 2; % 2 or 3
end
seti.fileSuffix = sprintf('_testGridScale%dD',dim);

%%
% *Dimensions of test matrices*
%
% * nROI: length of full matrix (upscaled)
% * nInvv = length of down scaled matrix

if dim == 3
    nROI = 32; % 3D
    nInvv = 15; % 3D (use nInvv to test different nInv (as vector))
else
    nROI = 512;
    %nInvv = 3:1:nROI; % test different nInv (vector as input is possible)
    nInvv = 127;
end

ndu = 2; % number of down- and upscalings in 2)
% ndu > 1 should have the same result as ndu = 1

% test imresize
%gridUp = @(B,nROI,nInv) imresize(B,nROI/nInv);
%gridDown = @(A,nROI,nInv) imresize(A,nInv/nROI);

% test up- and downscaling

disp(' ')
%%
disp(' 1) test up- and downscaling')

D = zeros(1,length(nInvv));

for ii = 1:length(nInvv)

    nInv = nInvv(ii);
    %nROI = 5;
    %nInv = 3;

    if dim == 2
      A = zeros(nROI,nROI);
      B = zeros(nInv);
    elseif dim == 3
      A = zeros(nROI,nROI,nROI);
      B = zeros(nInv,nInv,nInv);
    end

    % A: test matrix
    if 0
    for i = 1:nCD
        for j = 1:nCD
            A(i,j) = i*10+j;
        end
    end
    end

    % B: test matrix (2D and 3D)
    if 1
    for i = 1:nInv
        for j = 1:nInv
            if dim == 3
                for k = 1:nInv % 3D
                    B(i,j,k) = i*100+j*10+k;
                end
            end
        B(i,j) = i*10+j; % 2D
        end
    end
    end

    %B

    % test up- and downscaling
    disp('grid upscaling')
    tic
    A = gridUp(B,nROI,nInv); % grid up-scaling
    toc
    disp('grid downscaling')
    tic
    B2 = gridDown(A,nROI,nInv); % grid down-scaling
    toc

    fprintf('dim (length(size(A)): %g | dim (length(size(B2))): %g\n',length(size(A)),length(size(B)))
    
    if 0
        B
        A
        B2
        unique(B-B2); % is 0 if everything is correct...
        % Id(.) = gridDown(gridUp( . )), so B = gridDown(gridUp(B))
    end

    D(ii) = max(max(max(abs(B-B2)))); % should be 0 (then correct)
    fprintf('nInv = %g, D = %g | ',nInv,D(ii))
    if floor(ii/5) == ii/5
        fprintf('\n')
    end

end
fprintf('\n');

%%
% *Test down- and upscaling (several times, only in 2D)*

if dim == 2

    disp(' ')
    disp(' 2) test down- and upscaling (several times...) (2D)')
    fprintf('number of down- and upscalings: ndu = %g\n',ndu)

    % Smooth picture...
    % Probability density function of two-variate (or two dimensional) normal distribution (cf. [1, p. 845])
    % [1] Bronstein; Semendjajew; Musiol; M&uuml;hlig (2008): _Taschenbuch der Mathematik._
    %     7., vollst&auml;ndig &uuml;berarbeitete und erg&auml;nzte Auflage. Frankfurt: Harri Deutsch.

    pdf_twovarnormdist = @(mu1,mu2,sig1,sig2,rho,x1,x2) 1/(2*pi*sig1*sig2*sqrt(1-rho*rho)).*exp(-1/(2*(1-rho*rho))*((x1-mu1).^2./(sig1.^2)+(x2-mu2).^2./(sig2.^2)-2*rho*(x1-mu1).*(x2-mu2)/(sig1*sig2)));

    gv = linspace(-5,5,nROI);
    [x,y] = meshgrid(gv);
    mu1 = 1; mu2 = -1.5; sig1 = 1; sig2 = 1; rho = 0;
    z1 = pdf_twovarnormdist(mu1,mu2,sig1,sig2,rho,x,y);
    mu1 = -1; mu2 = 0.5; sig1 = 1.5; sig2 = 1.5; rho = 0;
    z2 = pdf_twovarnormdist(mu1,mu2,sig1,sig2,rho,x,y);
    z = z1+1/2*z2;
    %surf(x,y,z);

    A = z;
    normA = norm(A);

    % mexican hat 2D:
    % source: http://de.mathworks.com/help/wavelet/ref/mexihat.html#examples
    % lb = -5;
    % ub = 5;
    % N = 100;
    % [psi,xval] = mexihat(lb,ub,N);
    % plot(xval,psi)
    % title('Mexican Hat Wavelet');

    err = zeros(length(nInvv),1); % store error
    for ii = 1:length(nInvv)
        A2 = A;
        nInv = nInvv(ii);
        for idu = 1:ndu % iteration: several down- and up-scalings...
            disp('gridDown...')
            tic
            B = gridDown(A2,nROI,nInv); % grid down-scaling
            toc
            disp('gridUp...')
            tic
            A2 = gridUp(B,nROI,nInv); % grid up-scaling
            toc
        end
        err(ii) = norm(A-A2)/normA; % relative error
        %D(ii) = max(max(max(abs(A-A2))));

        fprintf('nInv = %g, err = %g | ',nInv,err(ii))
        if floor(ii/3) == ii/3
            fprintf('\n')
        end
    end
    fprintf('\n');

    figure(31); imagesc(A); axis xy; title('A original (2D)'); colorbar;
    figure(32); imagesc(B); axis xy; title('B = gridDown(A) (2D)'); colorbar;
    figure(33); imagesc(A2); axis xy; title('A2 = gridUp(gridDown(A)) (2D)'); colorbar;

    savePngFig(31,0,seti);
    savePngFig(32,0,seti);
    savePngFig(33,0,seti);

end
