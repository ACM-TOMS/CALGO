function plotm( varargin)
% plotm( data, [linespec], [options])
% Wrapper function which calls stem/plot/plot3 depending on the dimension of varargin{1}
% Also handles interval data
%
% Input:
% =============
%   data        dim x N array
%   [linespec]  Matlab Linespec string (e.g. 'r--' or 'xk')
%               The color specification might not work in many cases.
%   [options]   nearly everything which can be parsed by plot/plot2/plot3/stem/stem3/scatter3/mesh3/contour3/...
%
% Helper options:
% ==================
%   'equal',str             (experimental), chars, default='', Unifies axes
%                           If given, all properties (abbreviated) in str of the plots will be set to the same value
%                           str is a concatenation of the charaters: 'x', 'y', 'z', 'c' (for x, y, z, c-axis), 
%   'link',str              (experimental), chars, default='', Links xyz-axes and camera position
%                           If given, all properties (abbreviated) in str of the plots will be linked
%                           str is a concatenation of the charaters: 'x', 'y', 'z', 'c' (for x, y, z, c-axis), 'v' (for camera view properties)
%   'imp',[imp val]         1x2 vector, default=[0 0], Plots only if imp>=val
%   'verbose',val           Verbose level
%   'label',val             (experimental), default=enabled, adds labels to the axis
%   
%
% Options:
% =============
%   'box',val               int, Plots the hypercube of dimension val with volume 1
%   'spherical',[type off]  [int/string double], Plots points given in hyperspherical coordinates of type type (default: t=0), with offset off from the unit-sphere
%   'unit',str              string, default='', 
%                           Transforms last row prior of data prior further processing, depending on str. Possible values for str are:
%                           'abs', 'lin2log', 'lin2abslog', 'log2lin', 'log2abslin', 'lin2db', 'lin2absdb', 'db2lin', 'db2abslin', 
%   'xkcd'                  xkcd-ify all plots
%                       
% Options for 1-dim:
% ========================
%   'height',val            scalar or 1x2 vector, default=[0 1],
%                           determines the length of the lines to be plotted
%                               scalar: Line goes from 0 to val
%                               vector: Line goes from val(1) to val(2) 
%
% Options for 2-dim:
% ========================
%       'hull'              Plots also the convex hull
%       'boundary',val      Plots also the boundary 
%                               -1 <= val <= 0: computed with boundary (matlab function)
%                               val < -1:       points which have few neighbours are computed
%                               0 < val:        another strange method using Delauney triangulation
%
% Options for 3-dim:
% ========================
%       'contour'           Plots contour lines
%       'resolution',val    scalar, default depends on the input
%                           determines the resolution
%                               If val==0, the points are plotted
%                               otherwise val determines the resolution of the interpolated grid
%       'hull'              Plots also the convex hull
%       'boundary',val      Plots also the boundary 
%
% Options for 4-dim (experimental):
% ==================================
%
%
% Options for subdivision schemes
% ===================================
%
%
% Info:
% ===========
%   1-dim data:
%       Each point is plotted as a dot at height 1. 
%       The height can be changed by Option <'height',val>
%       To change the marker, the option <'Marker','something'> should be used
%
%   2-dim/3-dim data 
%       A point cloud is plotted
%
%
% E.g.: clf; hold on; clear all; load trimesh3d; plotm( [x y z]','resolution',500); plotm( [x y z]','resolution',0); view(3); plotm('rotate',0);
%       plotm(randn(2,40),'.-')
%       plotm(randn(1,40),'height',4);
%
%
% See also: surfm, plot, plot3, stem
%
%
% Written by: tommsch, 2016

% Changelog: tommsch,   2019-04-01      Added option 'hull', 'boundary' for plotm3
%                                       Added function plotm4
%            tommsch,   2019-11-25      Refactored functions to remove code duplication
%            tommsch,   2019-12-06      Added experimental options 'spherical','unit','imp'
%            tommsch,   2020-01-08      Added parsing of Matlabs LineSpec
%                                       Improved plotm4
%            tommsch,   2020-01-24      Added experimental options 'link' and 'equal'
%            tommsch,   2020-03-06      Added experimental option 'cylindrical' (expects coordinates in order
%            tommsch,   2020-04-24      Added experimental option 'delaunay'

% XX make plotm plot sequences
% XX write tests for unitflag, spherical, cylindrical
% XX add option 'filename',filename

    parg = varargin;
    [imp,varargin] = parsem( {'imp'}, varargin,[0 0]); %whether to plot or not
    if(imp(1)<imp(2) );
        return; end;

    %Parsing
    [verbose,varargin] = parsem( {'verbose','v'}, varargin, 0 ); %Verbose level
    [xkcdflag,varargin] = parsem( {'xkcd'}, varargin ); %Verbose level
    [linkaxisflag,varargin] = parsem( {'link','linkaxis'}, varargin, '' ); %pointcloud
    [sameaxisflag,varargin] = parsem( {'same','sameaxis','equal','equalize','equalizeaxis','equalaxis'}, varargin, '' ); %pointcloud
    [height,varargin] = parsem( {'height','h'}, varargin, [] ); %height
    [pointflag,varargin] = parsem( {'point','points','pt','pts'}, varargin ); %pointcloud
    [contourflag,varargin] = parsem( 'contour', varargin ); %contour
    [surfaceflag,varargin] = parsem( 'surface', varargin ); %surface
    [hullflag,varargin] = parsem( 'hull', varargin ); %hull
    [delaunayflag,varargin] = parsem( 'delaunay', varargin ); %delaunay triangulation
    [boundaryflag,varargin] = parsem( 'boundary', varargin, [] ); %boundary
    [res,varargin] = parsem( {'resolution','res'}, varargin, -inf );
    [LABEL,varargin] = parsem( {'label'}, varargin, 0 );
    [MARKERSIZE,varargin] = parsem( {'MarkerSize','markersize','ms'}, varargin, 5 );
    [boxflag,varargin] = parsem( 'box', varargin,0 );     
    [unitflag,varargin] = parsem( {'unit'}, varargin, [] );    
    if( ~isempty(unitflag) );
        switch unitflag;
            case {''};
                %do nothing
            case {'abs'};
                varargin{1}(end,:) = abs( varargin{1}(end,:) );
            case {'lin2log'};
                varargin{1}(end,:) = log( varargin{1}(end,:) );
            case {'lin2abslog'};
                varargin{1}(end,:) = abs( log(varargin{1}(end,:)) );
            case {'log2lin'};
                varargin{1}(end,:) = exp( varargin{1}(end,:) );
            case {'log2abslin'};
                varargin{1}(end,:) = abs( exp(varargin{1}(end,:)) );
            case {'lin2db','lin2dB'};
                varargin{1}(end,:) = 10*log10( varargin{1}(end,:) );
            case {'lin2absdb','lin2absdB'};
                varargin{1}(end,:) = abs( 10*log10(varargin{1}(end,:)) );
            case {'db2lin','dB2lin'};
                varargin{1}(end,:) = 10.^(varargin{1}(end,:)./10 );
            case {'db2abslin','dB2abslin'};
                varargin{1}(end,:) = abs( 10.^(varargin{1}(end,:)./10) );
            otherwise
                warning( 'plotm:unit', 'Wrong value for ''unit''. Option ignored.' ); end;
    end;
    [cylindricalflag,varargin] = parsem( {'cylindrical','cylinder','hypercylindrical','hypercylinder'}, varargin, [] );
    if( ~isempty(cylindricalflag) )
        varargin{1} = varargin{1}([2:end 1],:); 
        cylindricalflag = [cylindricalflag -inf]; end;
    [sphericalflag,varargin] = parsem( {'spherical','sphere','hyperspherical','hypersphere','hs'}, varargin, cylindricalflag );
    if( ~isempty(sphericalflag) )
        if( numel(sphericalflag)==1 );
            sphericalflag(2) = 0; end;
        if( isnan(sphericalflag(2)) );
            sphericalflag(2) = min( varargin{1}(end,:) ); end;
        if( sphericalflag(2)==-inf ); 
            varargin{1} = sph2cartm( varargin{1}, sphericalflag(1) );
        elseif( sphericalflag(1)<inf )
            varargin{1}(end,:) = varargin{1}(end,:)-sphericalflag(2);
            varargin{1} = sph2cartm( varargin{1}, sphericalflag(1) );
        elseif( sphericalflag(1)==inf ); 
            varargin{1}(end,:) = sph2cartm2( varargin{1}(1:end-1,:), sphericalflag(1) );
        else
            warning( 'plotm:spherical', 'Wrong value for ''spherical''. Option ignored.' ); end;
    end
    [LINESPEC,COLORSPEC] = deal('');
    if( numel(varargin)>=2 ); 
        [LINESPEC,LINESTYLE,~,COLORSPEC] = parselinespec( varargin{2} );
        if( ~isempty(LINESPEC) ); 
            varargin(2) = []; end; end;
    if(isempty(COLORSPEC)); 
        val = get( gca, 'colororder' ); val = val(1,:);
    else
        val = COLORSPEC; end;
    [COLOR,varargin] = parsem( {'Color','color','c','Colour','colour'}, varargin, val );    

    %Preprocessing
    switch numel( boxflag )
        case 1;
            switch boxflag
                case 0;
                    %do nothing
                case 1;
                    plotm( [0 1], LINESPEC );
                case 2;
                    if( isempty(LINESPEC) );
                        LINESPEC = '-.'; end;
                    plotm( [0 0;1 0;1 1;0 1; 0 0]', LINESPEC, 'point' );
                case 3;
                    %plotm( [0 0 0;1 0 0;1 1 0;0 1 0; 0 0 0;0 0 1;1 0 1;1 1 1;0 1 1; 0 0 1;0 0 0; 1 0 0;1 0 1;0 0 1; 0 0 0;0 1 0;1 1 0;1 1 1;0 1 1; 0 1 0;0 0 0;0 1 0;0 1 1;0 0 1;0  0 0;1 0 0;1 1 0;1 1 1;1 0 1;1  0 0;0 0 0]','hull');
                    plotm( [0 0 0; 0 0 1;0 1 0;0 1 1;1 0 0;1 0 1;1 1 0;1 1 1].','hull');
                otherwise;
                    if( isempty(LINESPEC) );
                        LINESPEC = 'x'; end;
                    plotm( mixvector([0 1],boxflag), LINESPEC ); end;
        case 2;
            plotm( boxflag );
        case 4;
            if( isempty(LINESPEC) );
                        LINESPEC = '-'; end;
            x1 = boxflag(1); 
            x2 = boxflag(2); 
            y1 = boxflag(3); 
            y2 = boxflag(4);
            plotm( [x1 y1; x1 y2; x2 y2; x2 y1; x1 y1].',LINESPEC,'c',COLOR);
        otherwise
            warning( 'plotm:box', 'Wrong value for ''box''. Option ignored.' ); end;

    %Plot functions
    if( ~isempty(varargin) ); 
        if( isS(varargin{1}) )
            blf( varargin{:} );
        elseif( iscell(varargin{1}) );
            dim = size( varargin{1}{1}, 1 );
            if( dim==1 ); 
                plotinterval1( varargin{:} );
            else; 
                warning( 'plotm:dim', 'no intervalplot for dim>1.' ); end;
                
        elseif( isa(varargin{1},'sequence') )
            dim = ndims( varargin{1} );
            if( dim==1 ); 
                plotsequence1( parg{:} );
            elseif( dim==2 ); 
                plotsequence2( parg{:} );
            else; 
                warning( 'plotm:dim', 'no sequenceplot for dimension %i', dim ); end;
        else
            dim = size( varargin{1}, 1 );
            if( dim==1 ); 
                plotm1( varargin{:} );
            elseif( dim==2 ); 
                plotm2( varargin{:} ); 
            elseif( dim==3 ); 
                plotm3( varargin{:} );
            elseif( dim==4 ); 
                plotm4( varargin{:} );
            else; 
                warning( 'plotm:dim', 'no ''plotm'' for dimension %i', dim ); end; end; end;
    
    if( iscell(LABEL) );
        xlabel( LABEL{1} );
        ylabel( LABEL{2} );
        zlabel( LABEL{3} );
    elseif( ~isequal(LABEL,0) );
        xlabel( 'x' );
        ylabel( 'y' );
        zlabel( 'z' ); end;
    
    if( ~isempty(sameaxisflag) ); 
        sameaxis( sameaxisflag ); end;    
    
    if( ~isempty(linkaxisflag) );
        linkaxis( linkaxisflag ); end;
    
    if( xkcdflag );
        xkcd; end;

    function plotinterval1( Q, varargin )
        
        if( isempty(height) );
            height = 0; end;

        n = numel(Q);
        if( numel(height)==1 ); 
            height = [height height+1]; end;
        X = vertcat( Q{:} )'; 
        X = [X; flipud(X)];
        Y = [height(1)*ones( 2, n ); height(2)*ones( 2, n )];
        try;
            patch( 'XData',X, 'YData',Y, varargin{:} );
        catch;
            warning( 'plotm:option', 'Wrong plot options. All options are ignored.' );
            patch( 'XData',X, 'YData',Y ); end;

        if( verbose>=1 );
            val1 = sum( diff(vertcat(Q{:}),1,2) );
            val2 = sum( abs(diff(vertcat(Q{:}),1,2)) );
            if( val1~=val2 );
                if( verbose>=2 )
                    fprintf( 'Area (oriented/unoriented): %i / %i\n', val1, val2 ); end;
            else;
                if( verbose>= 2);
                    fprintf('Area: %i\n',val1); end; end; end;

    end

    function plotsequence1( c, varargin);
        plotm( [((c.idx):(c.idx+numel(c)-1) );(c.c).'], varargin{:} );
    end

    function plotsequence2( c, varargin);
        Z = c.c;
        Z = Z(:).';
        Z(Z==0) = [];
        if( isempty(LINESPEC) );
            plotm( [supp(c); Z], varargin{:} );
        else
            plotm( [supp(c); Z], LINESPEC, varargin{:} ); end;
    end

    function plotm1( Q, varargin )
        if( isempty(height) );
            height = 1; end;
       
        if( isreal(Q) )
            try
                stem( Q, height*ones(1,size(Q,2)), LINESPEC, varargin{:}, 'Color',COLOR, 'MarkerSize', MARKERSIZE); 
                %quiver kann option 'rx', etc. nicht parsen, deshalb kann es nicht verwendet werden.
            catch
                warning( 'plotm:option', 'Wrong plot options. All options are ignored.' );
                stem( Q, ones(1,size(Q,2)), LINESPEC, 'Color',COLOR, 'MarkerSize',MARKERSIZE); end;
        else
            plotm( [real(Q); imag(Q)], LINESPEC, 'Color',COLOR, 'MarkerSize',MARKERSIZE, varargin{:} );
            xlabel( 'Re' );
            ylabel( 'Im' ); end;
    end

    function plotm2( Q , varargin)
        try 
            if( hullflag )
                idx = convhull( Q(1,:), Q(2,:) );
                try; 
                    plot(Q(1,idx),Q(2,idx), LINESPEC, 'Color',COLOR, 'MarkerSize', MARKERSIZE, varargin{:} ); 
                catch; 
                    patch( Q(1,idx), Q(2,idx), COLOR, varargin{:} );
                end;
            elseif( delaunayflag )
                DT = delaunay( Q(1,:), Q(2,:) );
                triplot( DT, Q(1,:), Q(2,:), LINESPEC, varargin{:} );
            elseif( ~isempty(boundaryflag) )
                if( -1<=boundaryflag && boundaryflag<=0 );
                    boundaryflag=-boundaryflag;
                    idx = boundary(Q(1,:).',Q(2,:).',boundaryflag); %compute boundary using Matlabs function
                    plot(Q(1,idx),Q(2,idx), LINESPEC, 'Color',COLOR, 'MarkerSize', MARKERSIZE, varargin{:} );                
                elseif( boundaryflag<-1 );
                    boundaryflag = -boundaryflag-1;
                    D = pdist2( Q.', Q.', 'cityblock', 'Smallest',10); %compute distance to other points
                    DD = sum( D, 1 ); %compute something like an average distance
                    DDD = cumsum( sort(DD-mean(DD)) );
                    [~,idx] = min( DDD );
                    bar = abs( DDD(idx+1)-DDD(idx) )+mean( DD )*(boundaryflag-1);
                    idx = DD>bar; %choose vertices which have high distance to other points. Thus, they are likely to be at the boundary
                    plot( Q(1,idx), Q(2,idx), LINESPEC, 'Color',COLOR, 'MarkerSize',MARKERSIZE, varargin{:} );
                else            
                    DT = delaunayTriangulation( Q.' ); %compute delauny triangulation
                    AA = DTarea( DT ).';               %compute area of triangles
                    AAA = cumsum( sort(AA-mean(AA)) ); %compute something like an average value of the areas
                    [~,idx] = min( AAA );
                    bar = AAA(idx(1)+1)-AAA(idx(1))+mean( AA )*(boundaryflag);
                    idx = AA<=bar; %select triangles with area less than bar
                    EDGE = [DT.ConnectivityList(idx,[1 2]); DT.ConnectivityList(idx,[2 3]); DT.ConnectivityList(idx,[3 1])]; %edges of chosen triangles
                    EDGE = sort( EDGE.' ).'; %remove all edges which occur multiple times
                    [~,~,idx] = unique(EDGE,'rows');
                    counts = histc(idx,unique(idx) );
                    counts = counts(idx);
                    idx = counts==1;
                    EDGE = EDGE(idx,:);
                    plot( [Q(1,EDGE(:,1) ); Q(1,EDGE(:,2))], [Q(2,EDGE(:,1) ); Q(2,EDGE(:,2))], LINESPEC, 'Color',COLOR, 'MarkerSize',MARKERSIZE, varargin{:} ); end;
            else
                if( numel(COLOR)==size(Q,2) )
                    holdflag = ishold;
                    hold on;
                    if( ~isempty(LINESTYLE) )
                        try; %#ok<TRYNC>
                            plot( Q(1,:), Q(2,:), [LINESTYLE, COLORSPEC], 'MarkerSize', MARKERSIZE, varargin{:} ); end; end;
                    scatter( Q(1,:), Q(2,:), MARKERSIZE, COLOR, varargin{:} );
                    if( ~holdflag );
                        hold off; end; 
                else
                    plot( Q(1,:), Q(2,:), LINESPEC, 'Color',COLOR, 'MarkerSize', MARKERSIZE, varargin{:} ); end; end;
        catch
            warning( 'plotm:option', 'Wrong plot options. All options are ignored.' );
            plot( Q(1,:), Q(2,:) ); end;
        xlabel( 'x' )
        ylabel( 'y' )
        zlabel( 'z' )
    end

    function plotm3( Q, varargin )
        if( numel(res)==1 ); 
            res = [res res]; end;
        x = Q(1,:).';
        y = Q(2,:).';
        z = Q(3,:).';
        Q = []; %#ok<NASGU> %free space

        minx = min( x ); maxx = max( x );
        miny = min( y ); maxy = max( y );
        %test if we plot surface or points.
        %we plot only points be default, if there are more z values for the same pair (x,y) in Q
        if( res(1)==-inf )
            res = [max( 33, ceil(sqrt(length(x))) ), max( 33, ceil(sqrt(length(y))) )]; end;
        resval = [round( x/((maxx-minx)/(3*res(1))) ), round( y/((maxy-miny)/(3*res(2))) )];
        resval = unique( resval, 'rows' );
        if( isequal(contourflag,0) && isequal(surfaceflag,0) && isequal(pointflag,0) && isequal(hullflag,0) && isempty(boundaryflag) )
            if( numel(x)-10>size(resval,1) ); 
                pointflag = 1;
            else
                surfaceflag = 1; end; end;


        if( ~isequal(surfaceflag,0) || ~isequal(contourflag,0) )
            xlin = linspace( minx, maxx, res(1) );
            ylin = linspace( miny, maxy, res(2) );
            [X,Y] = meshgrid( xlin, ylin );
            f = scatteredInterpolant( x, y, z, 'linear', 'none' ); %linear interpolation, none extrapolation
            Z = f( X, Y );
            if( isempty(Z) || ~isempty(LINESPEC) ); 
                surfaceflag = 0; contourflag = 0; pointflag = 1; end; end;
        

        
        if( pointflag )
            try;
                if( parsem('filled',varargin) );
                    scatter3( x, y, z, MARKERSIZE, COLOR, varargin{:} );
                else
                    scatter3( x, y, z, MARKERSIZE, COLOR, LINESPEC, varargin{:} ); end;
            catch; 
                try; 
                    plot3( x, y, z, LINESPEC, 'MarkerSize',MARKERSIZE, 'Color',COLOR, varargin{:} );
                catch
                    warning( 'plotm:option', 'Wrong plot options. All options are ignored.' );
                    scatter3( x, y, z ); end; end; end;
        
        if( ~isequal(surfaceflag,0) )
            try; 
                mesh( X, Y, Z, 'EdgeColor',COLOR, varargin{:} ); %interpolated
            catch; 
                warning( 'plotm:option', 'Wrong plot options. All options are ignored.' );
                mesh( X, Y, Z ); end; end; %interpolated
        
        if( ~isequal(contourflag,0) )
            try; 
                if( isequal(contourflag,1) )
                    contour3( X, Y, Z, LINESPEC, 'LineColor',COLOR, varargin{:} ); %interpolated
                else
                    contour3( X, Y, Z, contourflag, LINESPEC, 'LineColor',COLOR, varargin{:} ); end; %interpolated
            catch; 
                warning( 'plotm:option', 'Wrong plot options. All options are ignored.' );
                contour3( X, Y, Z ); end; end; %interpolated
        
        if( ~isequal(hullflag,0) )
            K = convhulln( [x y z] );
            try
                if( numel(COLOR)==numel(x) );
                    trisurf( K, x, y, z, COLOR, varargin{:} ); 
                else
                    trisurf( K, x, y, z, varargin{:} ); end;
            catch
                warning( 'plotm:option', 'Wrong plot options. All options are ignored.' );
                trisurf( K, x, y, z ); end; end;
        
        if( ~isempty(boundaryflag) )
            k = boundary( [x y z], boundaryflag );
            try
                if( numel(COLOR)==numel(x) );
                    trisurf( k, x, y, z, COLOR, varargin{:} );
                else
                    trisurf( k, x, y, z, varargin{:} ); end;
            catch
                warning( 'plotm:option', 'Wrong plot options. All options are ignored.' );
                trisurf( k, x, y, z); end; end;
        axis tight;
        
    end


    function plotm4( Q , varargin );
    if( numel(res)==1 ); 
        res = [res(1) res(1) res(1)]; 
    elseif( numel(res)==2 );
        res = [res(1) res(1) res(1) res(2)]; 
    end;
        
        x = Q(1,:).';
        y = Q(2,:).';
        z = Q(3,:).';
        v = Q(4,:).';
        clear Q;
        
        minx = min( x ); maxx = max( x );
        miny = min( y ); maxy = max( y );
        minz = min( z ); maxz = max( z );
        if( res(1)==-inf )     
            res=[max( 10, ceil(sqrt(length(x)))) 
                 max( 10, ceil(sqrt(length(y))))
                 max( 10, ceil(sqrt(length(z))))]; end;    
        xlin = linspace( minx, maxx, res(1) );
        ylin = linspace( miny, maxy, res(2) );
        zlin = linspace( minz, maxz, res(3) );
        [X,Y,Z] = meshgrid( xlin, ylin, zlin );
        f = scatteredInterpolant( x, y, z, v, 'linear', 'none' ); %linear interpolation, none extrapolation
        V = f( X, Y, Z );
        
        if( ~ishold ); 
            clf; end;
            
        MINV = minm( V ); MAXV = maxm( V );
        val = MINV*0.99+MAXV*0.01;
        MAXV = MINV*0.01+MAXV*0.99;
        MINV = val;
        if( numel(res)==4 );
            levelres = res(4);
        else
            levelres = max( 5, floor(max(res)^(1/3)) ); end;
        if( levelres<=1 );
            levelres = 2; end;
        if( ~ischar(COLOR) && ~isequal(size(COLOR,1),levelres) );
            colormap( hot(levelres) );
            COLOR = colormap;
        elseif( ischar(COLOR) )
            colormap( COLOR );
            COLOR = colormap;
        else
            %do nothing
        end;        
        p = cell( 1, levelres );
        if( ~isempty(LINESPEC) ); 
            scatter3(x,y,z,3,v,LINESPEC); end;
        for i = 0:levelres-1
            level = MAXV*(1-i/(levelres-1))+MINV*i/(levelres-1);
            alpha = (1*(1-i/(levelres-1))+0.2*i/(levelres-1))^2;
            p{i+1} = patch( isosurface(X,Y,Z,V,level) );
            %isonormals(X,Y,Z,V, p{i+1} );
            set( p{i+1}, 'FaceColor',COLOR(i+1,:), 'EdgeColor','none', 'FaceAlpha',alpha ); end;
        view( 3 ); 
        camlight( 'left' ); 
        lighting gouraud;
        colorbar;
    end

end

function A = DTarea( DT )
    N = size( DT.ConnectivityList, 1 );
    A = zeros( N, 1 );
    for n = 1:N
        idx = DT.ConnectivityList(n,:);
        p1 = DT.Points(idx(1),:);
        p2 = DT.Points(idx(2),:);
        p3 = DT.Points(idx(3),:);
        A(n) = 1/2*abs( (p2(1)-p1(1))*(p3(2)-p1(1))-(p3(1)-p1(1))*(p2(2)-p1(2)) ); end;
end

function linkaxis( linkaxisflag )
    persistent LinkX LinkY LinkZ LinkC LinkCam;
            
    if( isequal(linkaxisflag,1) || isequal(linkaxisflag,2) );
        linkaxisflag = 'xyzc'; end;
    if( isequal(linkaxisflag,3) );
        linkaxisflag = 'xyzcv'; end;    
    
    ax = findobj( 0, '-property', 'xlim', '-not', 'tag', 'colorbar' );
    
    [linkaxisflag,val] = containsstring( linkaxisflag, {'v','Cam','cam'} ); 
    if( ~isempty(val) ); 
        LinkCam = linkprop( ax, {'CameraUpVector', 'CameraPosition', 'CameraTarget'} ); 
        %setappdata( gcf, 'LinkCam', LinkCam );
    elseif( ~isnumeric(LinkCam) );
        removeprop( LinkCam, 'CameraUpVector' ); 
        removeprop( LinkCam, 'CameraPosition' ); 
        removeprop( LinkCam, 'CameraTarget' ); end;
    
    [linkaxisflag,val] = containsstring( linkaxisflag,{'x','XLim','X'} ); 
    if( ~isempty(val) ); 
        LinkX = linkprop( ax, {'XLim'} ); 
        %setappdata(gcf, 'LinkX', LinkX); 
    elseif( ~isnumeric(LinkX) );
        removeprop( LinkX, 'XLim' ); end;
    
    [linkaxisflag,val] = containsstring( linkaxisflag, {'y','YLim','Y'} ); 
    if( ~isempty(val) ); 
        LinkY = linkprop( ax, {'YLim'} ); 
        %setappdata(gcf, 'LinkY', LinkY);
     elseif( ~isnumeric(LinkY) );
        removeprop( LinkY, 'YLim' ); end;
    
    [linkaxisflag,val] = containsstring(linkaxisflag, {'z','ZLim','Z'} ); 
    if( ~isempty(val) ); 
        LinkZ = linkprop( ax, {'ZLim'} ); 
        %setappdata(gcf, 'LinkZ', LinkZ);
    elseif( ~isnumeric(LinkZ) );
        removeprop( LinkZ, 'ZLim' ); end;      
    
    [linkaxisflag,val] = containsstring( linkaxisflag, {'c','CLim','C'} ); 
    if( ~isempty(val) ); 
        LinkC = linkprop( ax, {'CLim'} ); 
        %setappdata( gcf, 'LinkC', LinkC);
    elseif( ~isnumeric(LinkC) );
        removeprop( LinkC, 'CLim' ); end;       
    
    if( sum(linkaxisflag)~=0 );
        warning( 'plotm:link', 'Wrong options for ''link''.' ); end;
end



function sameaxis( XYZC )
% SAMEAXES unifies/synchronizes axis limits on different axes and subplots.
%
% See also PBASPECT, DASPECT, LINKAXES, XLIM, YLIM, ZLIM.
% Created Jun/13 by Johannes Keyser (jkeyser@uni-osnabrueck.de)
% Revised Aug/13 by jkeyser: +arg HNDS to restrict axes search from parents
% Revised Okt/13 by jkeyser: +arg XYZC to generalize to any axis, +comments
% Revised Jan/14 by jkeyser: +polished for Matlab File Exchange publication
% Revised May/14 by jkeyser: +fixed check for handles, +exclusion example
% Changed by tommsch for plotm
    %XYZC = 'xyzc';
    if( isequal(XYZC,1) || isequal(XYZC,2) );
        XYZC = 'xyc'; end;    
    if( isequal(XYZC,3) );
        XYZC = 'xyzc'; end;        
    for xyzc = XYZC(:)' % iterate over x, y, z, c
        lim = [xyzc 'lim'];
        % find axes (== objects with color-limits - except colorbars)
        axs = findobj( 0, '-property',lim, '-not','tag','colorbar' );
        if isempty(axs)
            warning( 'plotm:noSuchLim', 'No children with "%s".', lim )
            continue; end;
        % get() and then set() the pooled min() & max() to all axes objects
        lims = get( axs, lim );
        if( iscell(lims) );
            lims = [lims{:}]; end % unpack if necessary
        set( axs, lim, [min(lims) max(lims)] ); end;
end

function [ LINESPEC, LINESTYLE, MARKER, COLOR ] = parselinespec( spec );
    if( ~ischar(spec) );
        LINESPEC = '';
        return; end;
    [spec,LINESTYLE] = containsstring( spec, {'--','-.',':','-'} );
    [spec,MARKER] = containsstring (spec, {'square','diamond','pentagram','hexagram','+','o','*','.','x','s','d','^','v','>','<','p','h',} );
    [spec,COLOR] = containsstring( spec, {'r','g','b','c','m','y','k','w'} );
    if( sum(spec)==0 );
        LINESPEC = strcat( MARKER, LINESTYLE, COLOR );
    else
        COLOR='';
        MARKER='';
        LINESTYLE='';
        LINESPEC = ''; end;
end

function [ strout, whatout ] = containsstring( strin, whatin );
    whatout = '';
    strout = strin;
    for i = whatin;
        val = strfind( strin, i );
        if( ~isempty(val) );
            whatout = i{1}; 
            strout(val:val+numel(i{1})-1)=0; 
            return; end;
    end;
end

function xkcd
    try;
        if( isequal(exist('xkcdify','file'),2) );
            xkcdify(gca);
        else
            warning( 'plotm:xkcdify', 'xkcdify seems not be installed.', lim ); end;
    catch;
        warning( 'plotm:xkcdify', 'Unkown error in xkcdify.', lim ); end;
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.