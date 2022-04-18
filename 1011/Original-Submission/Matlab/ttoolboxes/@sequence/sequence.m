% Class which represents the space \ell_0(\ZZ^s), s\in\NN
%
% Properties:
%   c       the sequence
%   idx     index of first entry in c
% 
% Info:
%   Due to stupid Matlab, the values of c can't be accesed using sequence(idx).
%   Instead sequence.c(idx) must be used
%
% Important functions:
%   norm
%   sum
%   supp
%   symbol
%   upsample       
%   conv
%   diffsequence
%
% Notes:
%   Indexing is zero based
%   Class requires m-package, but the functions of this package have no appended m
%   Due to strange Matlab-behaviour, Matlab cannot display cell arrays with sequences as entries
%

classdef sequence
    properties
        c       %the sequence
        idx     %index of first entry in c
    end

%Constructors, Get, Set
    methods
        function obj = sequence( c_, idx_ ) %Constructor 
            if( nargin==1 );
                obj.c = c_;
                obj.idx = zeros( ndimsm(c_), 1 );
            else
                if( isrow(idx_) );
                    idx_ = idx_.'; end;
                obj.c = c_;
                obj.idx = idx_;
                if( ndimsm(c_)>numel(idx_) ); 
                    warning( 'sequence:dim', 'Sequence is higher dimensional than index.' ); end; end;
        end
        
%Indexing
        function c = ref( obj, idx_ )
            %referencing elements in c. Each column in idx_ represents the indices for one element to be returned
            num=size(idx_,2);
            c=zeros( 1, num );
            for i=1:num
                idxcol_ = idx_(:,i);
                if( any(obj.idx>idxcol_) || any(obj.idxmax<idxcol_) )
                    c(i) = 0;
                else
                    idxcol_ = idxcol_-obj.idx+1;
                    idxcol_ = num2cell( idxcol_ );
                    c(i) = obj.c(idxcol_{:} );end; end;
        end
              
%Arithmetic Operators
        function obj = uplus( o1 ); 
            obj = o1; 
        end        
        function obj = plus( o1, o2 ); 
            [A,amin] = addsequence( o1.c, o2.c, o1.idx, o2.idx ); 
            obj = sequence( A, amin ); 
        end
        function obj = uminus( o1 ); 
            obj = sequence( -o1.c, o1.idx ); 
        end
        function obj = minus( o1, o2 ); 
            [A,amin]=addsequence( o1.c, -o2.c, o1.idx, o2.idx ); 
            obj = sequence( A, amin ); 
        end
        function obj = times( o1, o2 ); 
            [C,cmin] = equalizeminidx( {o1.c o1.idx;o2.c o2.idx} ); 
            obj = sequence( C{1}.*C{2}, cmin ); 
        end
        function obj = rdivide( o1, o2 ); 
            [C,cmin] = equalizeminidx( {o1.c o1.idx;o2.c o2.idx} ); 
            obj = sequence( C{1}./C{2}, cmin ); 
        end
        function obj = ldivide( o1, o2 ); 
            [C,cmin] = equalizeminidx( {o1.c o1.idx;o2.c o2.idx} );
            obj = sequence( C{1}.\C{2}, cmin ); 
        end
        function obj = power( o1, o2 ); 
            [C,cmin] = equalizeminidx( {o1.c o1.idx;o2.c o2.idx} ); 
            obj = sequence( C{1}.^C{2}, cmin ); 
        end

        function obj = mtimes( o1, o2 ); 
            if( isa(o1,'sequence') && isscalar(o2) && isnumeric(o2) ) %scalar multiplication
                obj = sequence(o1.c*o2,o1.idx);
            elseif( isscalar(o1) && isa(o2,'sequence')  && isnumeric(o1) ) %scalar multiplication
                obj = sequence( o1*o2.c, o2.idx );
            else
                warning( 'sequence:algebra', 'Multiplication of sequences is probably the wrong thing.' );
                [C,cmin] = equalizeminidx( {o1.c o1.idx;o2.c o2.idx} ); 
                obj = sequence( C{1}*C{2}, cmin ); end;
        end
        
        function obj = mrdivide( o1, o2 ); 
            if( isa(o1,'sequence') && isscalar(o2) && isnumeric(o2) ) %scalar division
                obj = sequence( o1.c/o2, o1.idx );
            elseif( isscalar(o1) && isa(o2,'sequence')  && isnumeric(o1) ) %scalar division
                warning( 'sequence:algebra', 'Dividing sequences is probably the wrong thing.' );
                obj = sequence( o1/o2.c, o2.idx );
            else     
                warning( 'sequence:algebra', 'Dividing sequences is probably the wrong thing.' );
                [C,cmin] = equalizeminidx( {o1.c o1.idx;o2.c o2.idx} ); 
                obj = sequence( C{1}/C{2}, cmin ); end;
        end
        
        function obj = mldivide( o1, o2 );
            if( isa(o1,'sequence') && isscalar(o2) && isnumeric(o2)) %scalar division
                obj = sequence( o1.c\o2, o1.idx );
            elseif( isscalar( o1 ) && isa(o2,'sequence')  && isnumeric( o1 )) %scalar division
                obj = sequence( o1\o2.c, o2.idx );
            else     
                warning( 'sequence:algebra', 'Dividing sequences is probably the wrong thing.' );
                [C,cmin] = equalizeminidx( {o1.c o1.idx;o2.c o2.idx} );
                obj = sequence( C{1}\C{2}, cmin ); end
        end
        
        function obj = mpower( o1, o2 );
            warning( 'sequence:algebra', 'Powers of sequences are probably the wrong thing.' );
            [C,cmin] = equalizeminidx( {o1.c o1.idx;o2.c o2.idx} );
            obj = sequence( C{1}^C{2}, cmin );
        end
        
        %function obj = transpose( o1 ); obj = sequence(o1.c.', o1.idx); end
        %function obj = ctranspose( o1 ); obj = sequence(o1.c', o1.idx); end

%Relational Operators
        function obj = eq( o1, o2 ); 
            [C,cmin] = equalizeminidx( {o1.c o1.idx;o2.c o2.idx} );
            obj = sequence( C{1}==C{2}, cmin ); 
        end
        function obj = ne( o1, o2 ); 
            [C,cmin] = equalizeminidx( {o1.c o1.idx;o2.c o2.idx} );
            obj = sequence( C{1}~=C{2}, cmin ); 
        end
        function obj = gt( o1, o2 ); 
            [C,cmin] = equalizeminidx( {o1.c o1.idx;o2.c o2.idx} );
            obj = sequence( C{1}>C{2}, cmin ); 
        end
        function obj = ge( o1, o2 ); 
            [C,cmin] = equalizeminidx( {o1.c o1.idx;o2.c o2.idx} );
            obj = sequence( C{1}>=C{2}, cmin ); 
        end
        function obj = lt( o1, o2 ); 
            [C,cmin] = equalizeminidx( {o1.c o1.idx;o2.c o2.idx} );
            obj = sequence( C{1}<C{2}, cmin ); 
        end
        function obj = le( o1, o2 ); 
            [C,cmin] = equalizeminidx( {o1.c o1.idx;o2.c o2.idx} );
            obj = sequence( C{1}<=C{2}, cmin ); 
        end        
        function r = isequal( o1, o2 ); 
            C = equalizeminidx({o1.c o1.idx;o2.c o2.idx} );
            r = isequal( C{1},C{2} );
        end  
        
% Functions for arrays
        function r = idxmin( obj );   
            r = obj.idx; 
        end
        function r = idxmax( obj );   
            r = obj.idx+size( obj ).'-1; 
        end
        function r = nnz( o1 ); 
            r = nnz( o1.c ); 
        end  
        function r = numel( o1 ); 
            r = numel( o1.c ); 
        end
        function r = ndims( obj ); 
            r = size(obj.idx,1); 
        end
        
        function r = size( obj, dim, N )
             if( nargin==1 );      
                 r = sizem( obj.c, [], obj.ndims ); %matlab cant display 1d-sequences with this size-function
             elseif( nargin==2 );  
                 r = sizem( obj.c, dim ); 
             elseif( nargin==3 );  
                 r = sizem( obj.c, dim, N ); end;
        end
        
        function varargout = equalize( varargin ) %sets the same idx for all sequences
            C = cell( nargin, 2 );
            for i = 1:nargin; 
                C{i,1} = varargin{i}.c; 
                C{i,2} = varargin{i}.idx; end;
            [C_,idx_] = equalizeminidx( C );
            varargout = cell( size(varargin) );
            for i = 1:nargin; 
                varargout{i} = sequence( C_{i}, idx_ ); end;
            varargout{nargin+1} = idx_;
        end
        
        function obj = shrink( obj ) 
            %removes zeros at the border, does not change logically the sequence
            [c_,origin_] = removezero( obj.c, 'border', 'dontkeepdim' );
            obj.c = c_;
            obj.idx = obj.idx+origin_-1;
        end
               
%usual mathematical functions for \ell(\ZZ^s)
        function r = norm(o1,p); 
            if( nargin==1 ); 
                p = 2; end;
            if( p==inf ); 
                r = max( abs(o1.c(:)) ); 
            elseif( p==-inf ); 
                r = 0; %finitely supported sequences always have zeros
            elseif( p==0 );
                r = nnz( o1 );
            else; 
                r = summ( abs(o1.c).^p, [] )^(1/p); end;
        end      
        
%display
        function disp( o1 );
              fprintf( 'class: sequence\nc=\n[\n' ); 
              disp( o1.c ); 
              fprintf('\b]\n' );
              fprintf( 'idx= ( ' ); 
              fprintf( '%i ', o1.idx ); 
              fprintf( ') \n' );
        end
          
        function plot( o, varargin );
            plotm( o, varargin{:} );
        end
          
%%%%%%%%%%%%%%%%%%%%
%The important stuff

        function r = sum( obj ); 
            r = sum( obj.c(:) ); 
        end
        function r = supp( obj ); 
            r = supp( obj.c, obj.ndims, obj.idx ); 
        end
        function r = symbol( obj ); 
            r = mask2symbol( obj.c, 'amin', obj.idx, 'dim', obj.ndims ); 
        end
        function oout = upsample( oin, M ); 
            [c_,idx_] = upsamplem( oin.c, oin.idx, M ); 
            oout = sequence( c_, idx_ ); 
        end
        
        function oout = conv( varargin )
            %c_=convm(varargin{:});
            c_ = convm( varargin{1}.c, varargin{2}.c );
            idx_ = varargin{1}.idx+varargin{2}.idx;            
            for i = 3:size( varargin, 2 )
                c_ = convm( c_, varargin{i}.c );
                idx_ = idx_+varargin{i}.idx; end;
            oout = sequence( c_, idx_ );
        end 
        
        function [obj] = characteristic( obj )
            obj.c = (obj.c~=0);
        end
        
        function [oout, mu] = diffsequence( oin, DIFF )
            [oout, mu] = diffsequence( oin.c, DIFF, 'cell' );
            for i = 1:length( oout ); 
                oout{i}= sequence( oout{i}, oin.idx ); end;
            if( numel(oout)==1 ); 
                oout = oout{1}; end;
        end
    end
end