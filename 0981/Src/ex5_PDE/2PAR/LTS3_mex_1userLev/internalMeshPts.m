function XY = internalMeshPts( NXYval, XYmin, XYmax )
%% INTERNAL MESH POINTS: P(i,j),  i,j = 2,...,N-1
%  Return a (col-wise) matrix XY of size (n^2, 2)
%  where n = NXYval-2.
%  The two columns of XY contain the abscissas and ordinates
%  of the internal mesh points
%
    n = NXYval-2;     % number of mesh points along each spatial dimension
    h = 1/(NXYval-1); % mesh step
    [Y,X]=meshgrid(linspace(XYmin+h,XYmax-h,n));
    X=X(:); Y=Y(:); % they are arrays of length n^2
    XY = [X Y]; % col-wise matrix
end

