function tjsr_plotoutput( type, varargin );
% tjsr_plotoutput( type, [options] );
%   Plots some output
%
% Input:
%   type        the type struct from tjsr containing type.opt.plot
%
% Options:
%   anything else may forwarded to functions which are called
%
% Available strings for type.opt.plot:
%   ['type_' 'zz1_zz2_ ... _zzn'], where 'zzi' is anything from type.cyclictree.zz      
%               plots type.cyclictree.zz1 to type.cyclictree.zzn
%               these data-things should have compatible size
%   'norm'      plots type.cyclictree.norm (in a nice way)
%   'polytope'  plots the polytope (in a nice way)
%   'L'         plots the number remaining vertices
%   'tree'      plots the vertices in a graph
%
% Output:
%   nice pictures
%
% Written by tommsch, 2018

% Changelog: 2019-03-30     Made double logarithmic plot for option 'norm'
%                           Added identifiers 'n', 'graph', etc...

if( isequal(type.opt.plot,'none') || isequal(type.opt.plot,0) ); 
    return; end; %fast return

try
    if( 1==strfind(type.opt.plot,'info_') );
        underscore = strfind(type.opt.plot,'_');
        underscore(end+1) = length( type.opt.plot )+1;
        name = cell( 1, length(underscore)-1 );
        what = cell( 1, length(underscore)-1 );
        for i = 1:numel( underscore )-1;
            name{i} = type.opt.plot(underscore(i)+1:underscore(i+1)-1);
            if( isempty(strfind(name{i},'.')) );
                what{i} = eval( ['type.cyclictree.' name{i}] );
            else
                what{i} = eval( ['type.' name{i}] ); end;
            if(iscell(what{i})); 
                what{i} = [what{i}{:}]; end;
            if( any(strfind(name{i},'norm')) || any(strfind(name{i},'dist')) )
                %do nothing
                %what{i}=log10(what{i});
                %name{i}=['log10(' name{i} ')'];
            elseif( any(strfind(name{i},'time')) || any(strfind(name{i},'dist')) );
                what{i} = log10( cumsum(what{i}) );
                name{i} = ['log10(' name{i} ')']; end; end;
        type.opt.plot = 'info'; end;

    switch type.opt.plot
        case 'info';
            clf; hold on; 
            if(any(diff(cellfun(@length, what)))); %different length
                for i=1:length(what);
                    plot(what{i},'.'); end;
                legend(name{:});
            else
                switch length( what )
                    case 1; 
                        plot( what{1}, '.' ); 
                        xlabel( 'Vertex-number' ); 
                        ylabel( name{1} ); 
                        title( name{1} );
                    case 2; 
                        plot( what{1}, what{2}, '.' ); 
                        xlabel( name{1} ); 
                        ylabel( name{2} ); 
                        title( [name{1} ' vs ' name{2}] );
                    case 3; 
                        plot3( what{1}, what{2}, what{3}, '.' ); 
                        xlabel( name{1} );
                        ylabel( name{2} );
                        zlabel( name{3} );
                        title( [name{1} ' vs ' name{2} ' vs ' name{3}] ); 
                        view( 3 );
                    otherwise; error( 'Too many axes to plot.' ); end; end;

        case {'o'};
            for i = 1:type.counter.numordering
                figure( i ); 
                clf;
                image( type.cyclictree.o{i}.', 'CDataMapping', 'scaled' ); end;

        case {'norm','N'}; 
            idx_alive = tjsr_getalivevertex( type, 1 );
            WARNING = warning( 'query','MATLAB:Axes:NegativeDataInLogAxis' );
            warning( 'off','MATLAB:Axes:NegativeDataInLogAxis' )
            norm_reshaped_live = cellfun( @(x,y,z) mat2cell(x.*z,1,y)', type.cyclictree.norm, type.cyclictree.L, idx_alive, 'UniformOutput',0 );
            norm_reshaped_live = [norm_reshaped_live{:}];
            norm_reshaped_live = arrayfun( @(x) [norm_reshaped_live{x,:}], 1:size(norm_reshaped_live,1), 'UniformOutput',0 );

            norm_reshaped_dead = cellfun( @(x,y,z) mat2cell(x.*(~z),1,y)', type.cyclictree.norm, type.cyclictree.L, idx_alive, 'UniformOutput',0 );
            norm_reshaped_dead = [norm_reshaped_dead{:}];
            norm_reshaped_dead = arrayfun(@(x) [norm_reshaped_dead{x,:}],1:size(norm_reshaped_dead,1),'UniformOutput',0);

            norm_lvl = cell2mat( arrayfun(@(x,y) repmat(x,1,y), type.cyclictree.normlvl, cellfun(@length,norm_reshaped_live), 'UniformOutput',0) );
            norm_live = [norm_reshaped_live{:}];
            norm_dead = [norm_reshaped_dead{:}];
            valy = [norm_lvl; norm_live; norm_dead];
            valy(valy<1) = 0;
            semilogy( log10(valy.'), '.' ); 
            yt = arrayfun( @(x) ['1+' num2str(10^(x)-1)], yticks, 'UniformOutput',0 );
            yticklabels(yt);
            title( 'Norm' );
            ylabel( 'norm' );
            xlabel( 'index of vertex' );
            warning( WARNING.state, 'MATLAB:Axes:NegativeDataInLogAxis' );

        case {'polytope','P','V'}
            VVraw = tjsr_getpolytope( type );
            dim = type.info.dim;
            if( type.info.dim>3 ); 
                vprintf( 'No polytope plot for dimension > 3. Plot only the first 3 dimensions.', 'imp',[3 type.opt.verbose] ); 
                chosendim = sort( randperm(dim,3) );
                VVraw = VVraw(chosendim,:);
                dim = 3;
            elseif( type.info.dim==3 );
                chosendim = [1 2 3];
            else
                chosendim = [1 2]; end;

            switch type.info.algorithm
                case TJSR_CONEFUNCT;
                    VV = VVraw;
                    for i = 1:dim;
                        val = VVraw;
                        val([1:i-1 i+1:end],:) = 0;
                        VV = [VV val]; %#ok<AGROW>
                        val = VVraw;
                        val(i,:) = 0;
                        VV = [VV val]; end; %#ok<AGROW>

                    VV = [VV zeros( dim, 1 )];
                case {TJSR_MINKFUNCT,TJSR_COMPLEXFUNCT}
                        VV = [VVraw -VVraw];
                otherwise
                    error( 'wrong algorithm' ); end;
            clear V;

            clf; hold on;
            if( dim>=3 ); 
                try
                    idx = convhulln( VV.' );
                    trisurf( idx, VV(1,:), VV(2,:), VV(3,:) );
                    %plotm(VV,'x'); 
                    alpha 0.8;
                    if( ~isequal(chosendim,[1 2 3]) );
                        xlabel( ['Axis ' num2str(chosendim(1))] );
                        ylabel( ['Axis ' num2str(chosendim(2))] );
                        zlabel( ['Axis ' num2str(chosendim(3))] );
                    else
                        xlabel( 'x' ); ylabel( 'y' ); zlabel( 'z' ); end;
                catch
                    %plotm(VV,'x');
                end


            elseif( dim==2 );
                clf; hold on;
                try
                    idx=convhull(VV.');
                    plot(VV(1,idx),VV(2,idx),'r-');
                    %plotm(VV,'bx','MarkerSize',5);
                catch
                    %plotm(VV,'bx','MarkerSize',1);
                end 

                %plot eigenplanes.
                 LENGTH = .5;
                 for i = 1:type.counter.numordering
                     if( type.cyclictree.smpflag(i)==2 ); 
                         continue; end;
                     p = type.cyclictree.V{i}(1:2,1);   %coordinates of corresponding point
                     n = type.cyclictree.Vs{i}(1:2,1);  %normal-vector of eigenplane
                     pa = p+LENGTH/2*[n(2); -n(1)];
                     pe = p-LENGTH/2*[n(2); -n(1)];
                     plotm( [pa pe], 'k-' )
                     plotm( p, 'ko' ); end; end;
            title( 'polytope' );
%             for i=1:type.counter.numordering
%                 lvlidx=unique(type.cyclictree.level{i});
%                 for j=lvlidx
%                     idx=type.cyclictree.level{i}==j;
%                     %plotm(type.cyclictree.V{i}(chosendim,idx),'Color',num2color(i+j),'Marker','.','LineStyle','none','MarkerSize',2*j+2);
%                     %plotm(type.cyclictree.V{i}(chosendim,idx),'Color',num2color(i+j),'Marker','.','LineStyle','none');
%                 end
%             end
%            legend;
            if( dim>=3 );
                view( 120, 30 ); end;

        case {'tree','graph'}
            for i = 1:type.counter.numordering
                if( type.cyclictree.L{i}~=0 );
                    figure( i ); 
                    clf;
                    ENTRY = cell(0,2);
                    %ENTRY = [ENTRY; {'rho'} type.cyclictree.rho{i}]; %#ok<AGROW> %add rho to label
                    ENTRY = [ENTRY; {'norm'} type.cyclictree.norm{i}]; %#ok<AGROW> %add norm to label
                    %ENTRY = [ENTRY; {'o'} {cellfun( @(x) mat2str(x'), cellfun(@(x) removezero(x,'compact'),num2cell(type.cyclictree.o{i},1),'UniformOutput',0), 'UniformOutput',0 )}]; %#ok<AGROW>
                    makeorderinggraph( type.cyclictree.o{i}, 'value',ENTRY, 'verbose',type.opt.verbose-2, 'labeldescription',0, 'labelnumber',1, varargin{:} );
                    valx = type.cyclictree.ordering{i}; 
                    if( iscolumn(valx) ); 
                        valx = valx.'; end;
                    title( ['tree (' num2str(i) '), ordering: ' num2str(valx)] ); end; end;

        case {'num','L'}
            val = type.cyclictree.livingvertex;
            Lval = cell2mat( type.cyclictree.L' );
            Lval(end+1,:) = sum( Lval, 1 );
            semilogy( [Lval; val].' ); 
            title( '(L) Number of added vertices, number of remaining vertices' );
            
        case 'none'
            %do nothing
            
        otherwise
            fprintf( 'Wrong argument for ''plot''.\n' ); end;
    drawnow;
catch
    fprintf( '''tjsr_plotoutput'': Error (Possibly out or memory).\n' ); end;

end

function idx_alive = tjsr_getalivevertex( type, cellflag )
% This function belongs to tjsr!
% returns linear indices of vertices which are not inside the polytope and whose children are not all computed yet
% 
% Input:
%   type        the type struct from tjsr
%   cellflag    returns cell array of indices for each tree
%
% Written by tommsch, 2018

idx_alive = cell( 1, type.counter.numordering );
for i = 1:type.counter.numordering
    idx_norm = type.cyclictree.norm{i}>1-type.opt.epslinprog; %vertices which are outside or at the border
    idx_nan = isnan( type.cyclictree.norm{i} ); %vertices without norm
    idx_status = type.cyclictree.status{i}==0; %vertices which are children
    val = type.cyclictree.parent{i}(idx_nan);
    idx_partlyfreshparent = false( size(idx_nan) );
    idx_partlyfreshparent(val) = true; %vertices whose children norms are not fully computed
    idx_alive{i} = idx_norm & (idx_partlyfreshparent | idx_status); end;

if(~cellflag)
    idx_alive = [idx_alive{:}]; end;

end


function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   
