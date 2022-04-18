function [ M ] = tjsr_domain(varargin)
% [ M ] = tjsr_domain( M, bound, [mhandle], [options] )
%
% Input:
%   M           6xN matrix, list of coordinates of vertices of triangles. Each column is one triangle and consists of [x1 y1 x2 y2 x3 y3]'.
%            or 2xN matrix, list of points of convex hull. Programm computes Delaunay triangulation from it.
%            or output of a run from this program
%
%   bound       double, bound of JSR which shall be achieved in each triangle
%   mhandle     function handle, default: @(p1,p2) costanzamatrix(2,p1,p2), function returning a set of transition matrices
%            or cell array of symbolic matrices in symvars "mu_ and la_" or "m and l".
%
%
% Options:
%   'verbose',val       Verbose level, default=1
%   'minsize',val       double, default=0
%   'minmaxco',val      double, default=[-inf inf -inf inf], minimum/maximum x/y values of all vertices of triangle needed to be considered
%   'matrix'            function handle, default: @(p1,p2) costanzamatrix(2,p1,p2), function returning a set of transition matrices
%   'plot'              Plots output
%   'onlyplot',val      default=0, Plots only, does not compute.
%   'save',val          string, default=0 (means no saving), Path where to save the output. If 'save'=1, output is saved in the folder of thise file.
%                       The filename is automatically generated.
%   'savename',val      string, default='', add 'savename' to the output file
%   'autosavetime',val  integer, default=60*60, saves output every 'autosavetime' seconds
%   'maxtriangle',val   integer, default=2^16, Maxmimum number of triangles to generate
%   'select',val        integer, default=0, which selection strategy to use for subdividing triangles
%                           0   subdivides triangles which have triangles with the same vertices and jsr less than bound
%                           1                        with large size/diameter/circumference
%                           2                        with smallest lower bound of JSR
%                           3                        which are near to convex hull of JSR<lb area
%                          -1   uses a random method
%   'tjsropt',cell      cell array, opts passed to tjsr
%
% Output:
%   M       List of coordinates of vertices of triangles/JSR/orderings/idx (see below)
%
%

% XX Make interface of jsr_pathcomplete

    %pre-initialization
    %%%%%%%%%%%%%%%%%%%%%%
    rng( str2double(datestr(clock,'HHMMSSFFF')) );
    randn(1000);
    
    %parse input
    %%%%%%%%%%%%%%%%%%%
    [verbose,varargin] = parsem( {'verbose','v'}, varargin, 1 );
    [tjsropt,varargin] = parsem( {'tjsropt'}, varargin,{});
    [minsize,varargin] = parsem( {'minsize','minsze'}, varargin, 0 );
    [minmaxco,varargin] = parsem( {'minmaxco'}, varargin, [-inf inf -inf inf]);
    [maxtriangle,varargin] = parsem( {'maxtriangle','maxnumtriangle','maxt','maxnumt'}, varargin, 2^16 );
    [plotflag,varargin] = parsem( {'plot','p'}, varargin );
    [onlyplotflag,varargin] = parsem( {'onlyplot'}, varargin );
    if( onlyplotflag ); 
        plotflag = 1; 
        verbose = 2; end;
    [saveflag,varargin] = parsem( {'save','s'}, varargin, 0 );
    [savename,varargin] = parsem( {'savename','sname','snme'}, varargin, datestr(clock,'HHMMSSFFF') );
    [autosave,varargin] = parsem( {'autosave','saveauto','sauto'}, varargin, 30*60 );    
    M = varargin{1};
    bound = varargin{2};
     varargin(1:2) = [];
    if( isa(varargin{1},'function_handle') )
        matrixhandle = varargin{1};
        varargin(1) = [];
    elseif( iscell(varargin{1}) && isa(varargin{1}{1},'sym') );
        matrixhandle = @( mu, la ) cellfun( @(x) subs(x,{'mu_' 'la_'},{mu, la}), varargin{1}, 'UniformOutput',false );
        varargin(1) = [];
    elseif( isequal(varargin{1},'B') )
        matrixhandle = @( mu, la ) costanzamatrix( -2, mu, la ); 
        varargin(1) = [];
    else
        if( isequal(varargin{1},'A') );
            varargin(1) = []; end;
        matrixhandle = @( mu, la ) costanzamatrix( 2, mu, la ); end;

    %initialize
    [select,varargin] = parsem( {'select','sel'}, varargin, -1 );
    
    if( ischar(M) );
        load(M,'M');
    elseif( size(M,1)==6 ) %
        %add missing data
        M = [M; zeros( 1, size(M,2) ); inf( 1, size(M,2) ); -ones( 1, size(M,2) ); trianglevalue( M, M, select, bound, verbose )];
    elseif( size(M,1)==SIZE );
        %do nothing
    elseif( size(M,1)==2 )
        %compute Delauney
        DT = delaunay( M(1,:), M(2,:) );
        M  = reshape( M(:,DT'), 6, [] );
        M  = [M; zeros( 1, size(M,2) ); inf( 1, size(M,2) ); -ones( 1, size(M,2) )];
        M  = [M ; trianglevalue( M, M, select, bound, verbose )];
    elseif( size(M,1)==11 ) %some old format
        M(9,:) = [];
        M(SIZE,:) = trianglevalue( M, M, select, bound, verbose );
    else
        if( size(M,2)==2 ); 
            vprintf( 'Did you transpose the input ''M''?\n', 'cpr','err' ); end;
        error( 'Wrong input format.' ); end;

    parsem( varargin, 'test' );
    
    val = saveoutput( M, bound, savename, 0, plotflag, minmaxco, verbose );
    if( onlyplotflag );
        return; end;
    vprintf( ['Job ' val ' started.\n'], 'imp',[1 verbose] );
    

    %Main Loop
    %%%%%%%%%%%%%%
    
    while( true )
        starttime = clock;
        M(SIZE,:) = trianglevalue( M, M, select, bound, verbose );
       
        indices = M(SIZE,:)>minsize & ...
                  M(IDX,:)<0 & (M(JSRUB,:)>bound) & ...
                  ~anym( M(COX,:)<minmaxco(1), 1 ) & ...
                  ~anym( M(COX,:)>minmaxco(2), 1 ) & ...
                  ~anym( M(COY,:)<minmaxco(3), 1 ) & ...
                  ~anym( M(COY,:)>minmaxco(4), 1 ) ;
        indices = find( fliplr(chooseval(fliplr(M(SIZE,:)), [10 50], fliplr(indices))) );
        if( size(M,2)>maxtriangle ); 
            break; end;
        if( any(indices) )
            val = saveoutput( M, bound, savename, 0, 0, minmaxco, verbose );
            vprintf( [repmat('.',[1 numel(indices)]) ' bound: %f, Current triangle-measure: %i, ' val '\n'], bound, max(M(SIZE,indices)), 'imp',[1 verbose] );
            for i = indices
                if( plotflag && verbose>=2 );
                    hold on;
                    fill( M(COX,i), M(COY,i), 'blue', 'EdgeColor','blue' )
                    drawnow; end;
                
                %asserts
                if( M(JSRUB,i)<bound ); 
                    vprintf( '\nError: Triangle with JSR<bound chosen.\n', 'cpr','err' ); end;
                                
                if( ~isequal(autosave,0) && etime(clock,starttime)>=autosave );
                    starttime = clock;
                    saveoutput( M, bound, savename, saveflag, plotflag, minmaxco, verbose ); end;
                
                ret = computedomainjsr( M(:,i), bound, matrixhandle, verbose, tjsropt );
                ret = computedomainjsr_post( ret, bound, verbose );
                
                if( plotflag && verbose>=2 );
                    hold on;
                    if( ret(2)<bound );
                        fill( M(COX,i),M(COY,i), 'green' )    
                    else
                        fill( M(COX,i),M(COY,i), 'red' ); end;
                    drawnow; end;                
                M(JSRLB,i) = ret(1);
                M(JSRUB,i) = ret(2);
                if( ret(2)<bound )
                    %flag triangle as finished
                    M(IDX,i) = 0;
                else 
                    M(IDX,i) = size(M,2)+1;
                    NT = makenewtriangle( M(CO,i) );
                    M = [M [NT; repmat( M(JSR,i), [1 size(NT,2)] ); -ones( 1, size(NT,2) ); ones(1,size(NT,2))]]; end; end; %#ok<AGROW>                
            vprintf( '\n', 'imp',[1 verbose] );
        else
            break; end;
        saveoutput( M, bound, savename, saveflag, plotflag, minmaxco, verbose );
    end
    
    %post-processing
    %%%%%%%%%%%%%%%%%%%
    
    val = saveoutput( M, bound, savename, saveflag, 0, minmaxco, verbose );
    vprintf( ['Job ' val ' finished.\n'], 'imp',[1 verbose] );
    
end

function [ filename ] = saveoutput( M, bound, savename, saveflag, plotflag, minmaxco, verbose )
    filename = ['tjsr_domain' '_'  num2str(savename) '_bd('  num2str( bound, 17 ) ')_time(' datestr( now, 'yyyy_mm_dd-hh_MM_ss' ) ')'];
    if( isequal(saveflag,1) )
        fullname = fullfile( fileparts(mfilename('fullpath')), filename );
    else
        fullname = fullfile( saveflag, filename ); end;

    if( plotflag )
        clf; hold on;
        %idx=M(IDX,:)~=-1;
        idx = M(IDX,:)==-1 | M(IDX,:)==0;
        idxlow = M(JSRUB,:)<bound & idx;
        idxup = M(JSRUB,:)>=bound & idx;
        fill( M(COX,idxup), M(COY,idxup), M(JSRLB,idxup), 'EdgeColor','Red' );
        fill( M(COX,idxlow), M(COY,idxlow), M(JSRUB,idxlow) );
        colorbar;
        set( gca, 'clim', [bound*.95 bound*1.005] );
        val = axis;
        val = [max( minmaxco(1), val(1) ) min( minmaxco(2), val(2) ) max( minmaxco(3), val(3) ) min( minmaxco(4), val(4) )];
        try; %#ok<TRYNC>
            axis(val); end;
        title( ['Bound: ' num2str(bound,17)] );
        drawnow; end;

    if( saveflag )    
        vprintf( '\nSave to disk: %s\n', fullname, 'imp',[2 verbose] );
        save( [fullname '.mat'], 'M' ); end ;

    if( saveflag && plotflag )
        FigList = findobj( allchild(0), 'flat', 'Type', 'figure' );
        for iFig = 1:length( FigList )
            FigHandle = FigList( iFig );
            print( FigHandle, [fullname '.png'], '-dpng' ); end; end;
end

function [ ret ] = trianglevalue( T, M, select, bound, verbose )
    %computes a value of a triangle
    if( size(T,2)==1 );
        ret = 1;
        return; end;
    
    vprintf( 'Selection based on: ', 'imp',[2 verbose] );
    
    val = deal( zeros(1,size(T,2)) );
    if( isequal(select,-1) );
        select = randi( 4,1,3 )-1; end;
    
    value = zeros( 0, size(T,2) );
    
    if( any(select==0) ); %distance to triangle with JSR<bound
        vprintf( 'distance ', 'imp',[2 verbose] );
        T = T(1:6,:);
        parfor i = 1:size( T, 2 );
            idx_1 = M(COX,:)==T(1,i) & M(COY,:)==T(2,i); %#ok<PFBNS>
            idx_2 = M(COX,:)==T(3,i) & M(COY,:)==T(4,i);
            idx_3 = M(COX,:)==T(5,i) & M(COY,:)==T(6,i);
            idx_JSR = M(JSRUB,:)<=bound ;
            val(i) = nnz( (idx_1 | idx_2 | idx_3) & idx_JSR )+1; end;   
        val = round(val/max(val)*5)/5*1.5;
        value = [value; val]; end;
    
    if( any(select==1) ); %size of triangle
        vprintf( 'size ', 'imp',[2 verbose] );
        X2 = T(3,:)-T(1,:);
        Y2 = T(4,:)-T(2,:);
        X3 = T(5,:)-T(1,:);
        Y3 = T(6,:)-T(2,:);
        val = sqrt( 1/2 *abs( X2.*Y3-X3.*Y2 ) ) + ... %area
               (sqrt(X2.^2+Y2.^2)+sqrt(X3.^2+Y3.^2)+sqrt((X3-X2).^2+(Y3-Y2).^2))/3;  %circumference 
        val = round(val/max(val)*3)/3;
        value = [value; val]; end;
    
    if( any(select==2) ); %least lower bound of JSR
        vprintf( 'JSR ', 'imp',[2 verbose] );
        if( size(T,1)>JSRLB );
            val = 1./T(JSRLB,:);
            val = round(val/max(val)*5)/5;
            val(isnan( val )) = 0;
        else
            val = ones( 1, size(T,2) ); end;
        value = [value; val]; 
    end;
    
    if( any(select==3) ); %in convex hull of jsr<bound domain
        vprintf( 'convex hull ', 'imp',[2 verbose] );
        idx = M(JSRUB,:)<bound;
        if( ~any(idx) );
            value = [value; ones( 1, size(T,2) )];
        else
            VV = reshape( M(CO,idx), 2, [] );
            k = convhulln( VV.' );
            VV = VV(:,unique( k ));
            pts = [mean( T(COX,:) ); mean( T(COY,:) )];
            val =  1./computepolytopenorm( pts, VV, 1 );
            val = round(val/max(val)*5)/5;
            value = [value; val]; end; end;
    
    ret = sum( value, 1 ) + 1;
     vprintf( '\n', 'imp',[2 verbose] );
end


function [ JSR ] = computedomainjsr( T, bound, matrixhandle, verbose, tjsropt );
    MAXSMPDEPTH = 8;

    A = [matrixhandle( T(1), T(2) ) matrixhandle( T(3), T(4) ) matrixhandle( T(5), T(6) )];
    [~,~,info] = findsmp( A, 'maxsmpdepth',MAXSMPDEPTH, 'v',verbose-2, 'd', 'bound', [0 bound] );
    JSR = info.jsrbound;
    if( JSR(1)>=bound );
        return; end;
    if( JSR(2)<bound );
        return; end;
       
    opts = jsrsettings( 'verbose',verbose-2 );
    [val,~,info] = jsr_pathcomplete( A, opts );
    if( info.status~=0 || info.stopFlag==2 || info.stopFlag==5 );
        val = [0 inf]; end;
    JSR(2) = min( JSR(2), max(val) );
    if( JSR(2)<bound );
        return; end; 
    vprintf( '!', 'imp',[1 verbose] );
    
    opts = jsrsettings( 'verbose',verbose-2 );
    [val,~,info] = jsr_opti_sos( A, opts );
    if( info.status~= 0);
        val = [0 inf]; end;
    JSR(2) = min( JSR(2), max(val) );
    if( JSR(2)<bound );
        return; end;
    vprintf( '!', 'imp',[1 verbose] ); 
    
    
    [val,info] = tjsr( A, 'bound',bound, 'v',verbose-2, 'delta',-1, 'maxtime',10*60, 'noclassify',1, tjsropt{:} ); %#ok<ASGLU>
    JSR(2) = min( JSR(2), max(val) );
    if( JSR(2)<bound );
        return; end;
end

function [ JSR ] = computedomainjsr_post( JSR, bound, verbose )
    if( numel(JSR)==1 );
        JSR = [JSR JSR]; end;
    if( JSR(2)<bound );
        vprintf( '_', 'imp',[1 verbose] ); 
        return; end;
    if( JSR(1)>=bound );
        vprintf( 'O', 'imp',[1 verbose] ); 
        return; end;
    
    vprintf( 'o', 'imp',[1 verbose] ); 
end



function T = makenewtriangle( T )
    %subdivides triangles into smaller triangles

    % 1 to 4 subdivision
    x1 = T(1); y1 = T(2);
    x2 = T(3); y2 = T(4);
    x3 = T(5); y3 = T(6);
    x4 = (x1+x2)/2; y4 = (y1+y2)/2;
    x5 = (x1+x3)/2; y5 = (y1+y3)/2;
    x6 = (x2+x3)/2; y6 = (y2+y3)/2;

    T=[ x1 y1 x4 y4 x5 y5;
        x2 y2 x6 y6 x4 y4; 
        x3 y3 x5 y5 x6 y6; 
        x4 y4 x6 y6 x5 y5]';
end

function val = CO;      val = 1:6;      end %coordinates of triangles, in the order [x1 y1 x2 y2 x3 y3]
function val = COX;     val = [1 3 5];  end %x-coordinates of triangles
function val = COY;     val = [2 4 6];  end %y-coordinates of triangles 
function val = JSR;     val = [7 8];    end %bounds of JSR of triangle
function val = JSRLB;   val = 7;        end %lb of JSR of triangle
function val = JSRUB;   val = 8;        end %upper bound of jsr of triangle
function val = IDX;     val = 9;        end %index of child. is set to -1 if it has no children yet. is set to 0 if jsr is less than bound
function val = SIZE;    val = 10;       end %value of triangle is always last row


function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 
