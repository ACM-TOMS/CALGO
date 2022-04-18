function S = getS(varargin)
% [ S ] = getS(dim || cellarray || 'name' || list, ['help'], [options])
% Returns subdivision operators in this packages format.
% Format: { mask, dilation-matrix,  digit-set, name}
%
% Input:
%       dim                     return all unnamed subdiv operators with the given dimension
%       cellarray               returns the given subdivision operators
%       list                    List of names plus corresponding values. Eg.: (< 'a',1/2*[1 2 1]','M',2 >) 
%       'name'                  names of subdivision operators
%
% Options:  
%       'help'                  returns the names of all named subdivision operators (and some more strings).
%       'supp'                  returns not the digit sets D, but the support of the masks as the digit sets.        
%       'OmegaRR'               returns not the digit sets D, but the support of the masks minus the digit sets as the digit sets.
%       'characteristic'        returns not the masks a, but the characteristic function of the digit sets as masks.
%       'nocheck'               does not make basic checks whether the subdivision operators are correct or not.
%       'bigcheck'              (experimental) checks some more stuff, but this is slower
%       'verbose',val           verbose level
%
% Info:
%       More elements (than: a,M,D, name) may be added to a subdivision scheme later. The name is always the last element in each row
%
% Output
%       S                       the subdivision operaotrs
%
% E.g.: getS(2)
%       getS('a',1/2*[1 2 1]','M',2)
%       getS('2_butterfly')
%       getS('1_all')
%
% See also: blf, tile, transitionmatrix, constructdigit, characteristic, supp, normalizeS, 
%
% Written by: tommsch, 2016
% For more information write to: <a href="tommsch@gmx.at">tommsch@gmx.at</a>

% XX bigcheck-sum-rules not tested

%#ok<*AGROW>

[bigcheck,varargin] = parsem( 'bigcheck', varargin );
[nocheck,varargin] = parsem( 'nocheck', varargin );
[verbose,varargin] = parsem( {'verbose','v'}, varargin, 1, 'expect',{'clop',[-1 inf]} );
[suppflag,varargin] = parsem( 'supp', varargin );
[OmegaRRflag,varargin] = parsem( {'OmegaRR','Om'}, varargin );
[characteristicflag,varargin] = parsem( 'characteristic', varargin );

    
    NUMELEMENTS = 4 ;%a subdivision scheme consists of 4 elements
    S = cell( 0, NUMELEMENTS );
    if( ~isequal(parsem( {'a','M','D','mask','dilation','digit'},varargin),0) );
        [S,varargin] = getS_namevalue( S, varargin{:} );  
    elseif( isscalar(varargin{1}) ); 
        %UNNAMED SUBDIVISION OPERATORS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dim = varargin{1}; 
        varargin(1) = [];
        if( dim==1 );
            %S{end+1,2} =4; S{end,3} = [0 1 8 9];
            S{end+1,2} = 3; S{end,1} = 1/8*[1 2 3 4 4 4 3 2 1]';
            %S{end+1,2} = 2; S{end,1} = 1/2*[1 2 1]'; S{end,end} = '1_unconnected';
        elseif( dim==2 );
            
            S{end+1,2} = [2 1; 0 2];
            S{end,3} = [0 0; 3 0; 0 1; 3 1]';
                    
        elseif( dim==3 );            
            S{end+1,2} = diag([3 2 2]); 
            a1 = 1/9*[1 3 6 7 6 3 1]'; 
            a2 = 1/8*[1 4 6 4 1]'; 
            a3 = 1/2*[1 2 1]'; 
            S{end,1} = convm( a1, a2, a3, 'outer' );
            S{end,end} = 'tensor3'; end; 
        
    elseif( iscell(varargin{1}) );
        S = varargin{1}; 
        varargin(1) = [];
    else;                           
        [S,varargin] = getS_named( S, varargin{:} ); end;
    
    parsem( varargin, 'test' );
    
%Post-processing
%if(numel(varargin)>1); vprintf('There are some unkown options among these: \n %v\n',varargin,'cpr','err','imp',[0 verbose]); end;    

if( ~isempty(S) );
    dim = size( S{1,2}, 2 );
else
    warning( 'getS:nodata', 'getS: No data returned.' ); end;
sizeS = size( S, 1 );
    
    for i = 1:sizeS;
        if( size(S,2)<NUMELEMENTS ); 
            S{1,NUMELEMENTS} = []; end;
        if( isempty(S{i,2}) ); 
            warning( 'getS:M', 'Missing dilation matrix for entry= %i, name: %s', i, S{i,end} ); end;
        if( isempty(S{i,3}) ); 
            S{i,3} = constructdigit( S{i,2} ); end;
        if( isempty(S{i,1}) || characteristicflag ); 
            S{i,1} = {}; [val1,val2] = characteristic( S{i,3} ); 
            S{i,1} = sequence( val1, val2 ); end;
        if( ~isempty(S{i,1}) && ~isa(S{i,1},'sequence')); 
            S{i,1} = sequence(S{i,1});end;  % make mask to svm if it is not a svm        
        if( suppflag ); 
            S{i,3}=S{i,1}.supp; end;
        if( OmegaRRflag ); 
            S{i,3} = setplus( S{i,1}.supp, -S{i,3} ); end; 
        if( isempty(S{i,end}) ); 
            S{i,end} = 'unnamed'; end; end;
    
    %Do Checks
    if( ~nocheck && sizeS>0 );  
        for i = 1:sizeS;
            if( min(abs(eig(S{i,2})))<=1+eps ); 
                warning( 'getS:notexpanding', 'Dilation matrix not expanding for entry: %i, name: %s', i, S{i,end} ); end;
            val = summ(S{i,1}.c,[]);
            detval = abs(det(S{i,2}) );
            if( abs(val-detval) > 100 * eps ) %%XX Aendern zu Berechnung fuer eps_\mathcal A
                warning( 'getS:sumrule0', 'Mask sums up to %f, instead to determinant=%i for entry: %i, name: %s', val, detval, i, S{i,end} ); end;
            if( isequal(sizem(S{i,2},1),1) && size(S{i,1}.c,2)>1 && size(S{i,1}.c,1)==1 ); 
                warning( 'getS:transpose', 'Seems like you forgot to transpose the mask sequence for entry: %i, name: %s', i, S{i,end} ); end;
            if( ~isequal(size(S{i,2}),[dim dim]) ); 
                warning( 'getS:dimM', 'Dilation Matrix has the wrong size for the entry: %i, name: %s', i, S{i,end} ); end;
            
            try
                if( bigcheck )                
                    %test digit set
                    Cl = constructdigit(S{i,2},'ZZ',S{i,3},'classify','sym'); %Classify entries in digit set
                    Cl = unique(Cl);
                    if( size(Cl,2)~=abs(det(S{i,2}))); 
                        warning( 'getS:digit', 'Digits are no digit set for entry: %i, name: %s\n', i, S{i,end} ); end;
                    %test sum rules
                    SUPP = supp(S{i,1});
                    Cl = constructdigit(S{i,2},'ZZ',SUPP,'classify','sym'); %Classify entries in mask
                    Dtemp = unique(Cl);
                    val = zeros( 1, numel(Dtemp) );
                    for j = 1:numel( Dtemp )
                        val(j) = summ( S{i,1}.c(SUPP(:,Cl==Dtemp(j))+S{i,1}.idx+1) ); end;
                    if( any(diff(val)) );
                        warning( 'getS:sumrule', 'Sum rules not fulfiled.' ); end; end;
            catch
                warning( 'getS:bigcheck', 'Bigcheck failed.' ); end; end; end;
    
end

%=============================================================================================================================
%=============================================================================================================================
%=============================================================================================================================


function [S,varargin] = getS_named(S,varargin)
        [val,varargin] = parsem( '1_all', varargin ); if( val ); [~,varargin] = parsem( 'parse_one',varargin,1,'set'); end;
        [val,varargin] = parsem( '1_rand', varargin ); if( val );                  S(end+1,:)=randS(1); end;
        [val,varargin] = parsem( '1_143', varargin ); if( val );                   S{end+1,2} = 2; S{end,1} = 1/4*[1 4 3]'; S{end,3} = [0 1]; S{end,end} = '1_143_scaled'; end; %just for testing the construction of matrices
        [val,varargin] = parsem( '1_1133', varargin ); if( val );                  S{end+1,2} = 2; S{end,1} = 1/4*[1 1 3 3]'; S{end,end} = '1_1133'; end; %dim U \neq dim V_0
        [val,varargin] = parsem( '1_ternary_C1', varargin ); if( val );            S{end+1,2} = 3; S{end,1} = 1/9*[1 3 6 7 6 3 1]'; S{end,end} = '1_ternary_C1'; end; %ternary, C1
        [val,varargin] = parsem( '1_unconnected_support', varargin ); if( val );   S{end+1,2} = 2; S{end,1} = 1/2*[1  2 1  0 -1 -2 -1 0 1 2 1]'; S{end,end} = '1_unconnected_support'; end;%closure(support) is not connected
        [val,varargin] = parsem( '1_singular_scheme', varargin ); if( val );       S{end+1,2} = 2; S{end,1} = [0.5 0.5 0.5 0.5]';  S{end,end} = '1_singular_scheme'; end;%Singular Scheme to starting seq [-1 1 -1 1 -1 1 -1 1]
        [val,varargin] = parsem( '1_4point', varargin ); if( val );                S{end+1,2} = 2; S{end,1} = 1/16*[-1 0 9 16 9 0 -1]'; S{end,end} = '1_4point'; end;%interpolatory, n>1, n==4: 4point scheme   
        [val,varargin] = parsem( '1_DD', varargin ); if( val );                    S{end+1,2} = 2; S{end,1} = 1/256*[3 0 -25 0 150 256 150 0 -25 0 3]'; S{end,end} = '1_DD'; end;% n>3, n=6: Dubuc Deslauriers
        [val,varargin] = parsem( '1_strange_interpolatory', varargin ); if( val ); %#ok<ALIGN>
            z=-1/32; %strange interpolatory, -1<z<0
            S{end+1,2} =4; 
            S{end,1} = 1/128*[8*z, 72*z, -72*z-15, -8*z-7, 9-24*z, 33-216*z, 216*z+110, 24*z+126, 24*z+126, 216*z+110, 33-216*z, 9-24*z, -8*z-7, -72*z-15, 72*z, 8*z]'; %strange interpolatory
            S{end,end} = '1_strange_interpolatory'; end;
        [val,varargin] = parsem( {'1_HD_3point','1_Hassan_Dodgson_3point'}, varargin ); if( val );              S{end+1,2} = 2; S{end,1} = 1/16*[1 5 10 10 5 1]'; S{end,end} = '1_Hassan_Dodgson_3point'; end; 
        
        [val,varargin] = parsem( '1_devil_stairs', varargin ); if( val );               S{end+1,2} = 3; S{end,1} = 1/2*[1 1 2 1 1]'; S{end,end} = '1_devil_stairs'; end;
        [val,varargin] = parsem( {'1_balanced_ternary','m0p'}, varargin ); if( val );   S{end+1,2} = 3; S{end,3} = [-1 0 1]; S{end,end} = '1_balanced_ternary'; end;
        
        
        [val,varargin] = parsem( {'1_HD_ternarybalanced_ternary','1_Hassan_Dodgson_balanced_ternary'}, varargin ); if( val );     S{end+1,2} = 3; S{end,1} = 1/9*[-1 0 2 8 9 8 2 0 -1]'; S{end,end} = '1_HD_ternarybalanced_ternary';  end;         
       
        [val,varargin] = parsem( {'1_11','1_spline_binary_0','1_haar'}, varargin ); if( val );   S{end+1,2} = 2; S{end,1} = [1 1]'; S{end,end} = '1_spline_binary_0'; end;
        [val,varargin] = parsem( {'1_121','1_spline_binary_1','1_hat'}, varargin ); if( val );   S{end+1,2} = 2; S{end,1} = [1/2 1 1/2]'; S{end,end} = '1_spline_binary_1'; end;
        [val,varargin] = parsem( {'1_1331','1_spline_binary_2'}, varargin ); if( val );          S{end+1,2} = 2; S{end,1} = 1/4*[1 3 3 1]'; S{end,end} = '1_spline_binary_2'; end;
        [val,varargin] = parsem( {'1_14641','1_spline_binary_3'}, varargin ); if( val );         S{end+1,2} = 2; S{end,1} = 1/8*[1 4 6 4 1]'; S{end,end} = '1_spline_binary_3'; end;            

        [val,varargin] = parsem( {'1_lizhang','1_lizhang_combined_ternary_5point'}, varargin ); 
        if( val ); %#ok<ALIGN>
            [param,varargin] = parsem( '1_lizhang', varargin, [] );
            if( isempty(param)); a = -1/243; b = -7/243; c = -2/81; d = -1/81; 
            elseif( size(param,2)==1 ); d = param(1); b = -7*d/3; c = 2*d; a = d/3; 
            elseif( size(param,2)==2 ); b = param(1); d = param(2); c = 2*d; a = (-b+c+d)/2;
            elseif( size(param,2)==3 ); b = param(1); c = param(2); d = param(3); a = (-b+c+d)/2;
            elseif( size(param,2)==4 ); a = param(1); b = param(2); c = param(3); d = param(4);            
            else; error('getS: ''1_lizhang'' needs a vector of zero to 4 entries as second argument'); end;
            S{end+1,2} = 3; 
            S{end,1} = 1/9*convm([1 1 1],[1 1 1],[-9*a   1/9-9*d+18*a   1/3-9*a-9*c+18*d   2/3-9*b+18*c+9*d   7/9+18*b-18*c-36*d   2/3-9*b+18*c+9*d   1/3-9*a-9*c+18*d   1/9-9*d+18*a   -9*a])';
            S{end,end} = '1_lizhang_combined_ternary_5point'; end;
        
        
        [val,varargin] = parsem( '1_paper1', varargin ); if( val );             S{end+1,2} = 2; S{end,1} = [3/8 1 3/4 0 -1/8]'; S{end,3} = [0 1]; S{end,end} = '1_paper1'; end; %fuer example in paper
        [val,varargin] = parsem( '1_three_disjoint_Om', varargin ); if( val );  S{end+1,2} = -2; S{end,1} = [1 0 0 1]'; S{end,3} = [0 3]; S{end,end} = '1_three_disjoint_Om'; end % Example where three smallest disjoint invariant sets Om exist: Om1= [0], Om2= [-1 2]; Om3= [ -2 1];
        [val,varargin] = parsem( '1_two_disjoint_Om', varargin ); if( val );    S{end+1,2} = 2; S{end,1} = sequence([1/2 0 0   1  0 0   1/2]',0); S{end,3} = [ 0 3]; S{end,end} = '1_two_disjoint_Om'; end; % Example where two smallest disjoint invariant sets Om exist: Om1= [0 3], Om2= [-2 -1 1 2 4 5]; Furthermore, JSR(T(Om1)|V_0)=1/2, JSR(T(Om2)|V_0)=1
        [val]=parsem( '1_daubechies', varargin ); if( val ); [N,varargin] = parsem( '1_daubechies',varargin,0); S{end+1,2} = 2; S{end,1} =double(daubechiesmask(N)).'; S{end,3} = [ 0 1]; end; %Daubechies wavelet
        [val,varargin] = parsem( '1_unbounded_tile', varargin ); if( val );     for i = 1:9;  S{end+1,2} = 3; S{end,3} = [0 3*i^i-1]; S{end,1} = [1 1 1]'; S{end,end} = ['1_unbounded_tile (' num2str(i) ')']; end; end; %unbounded tile, no digit set
        [val,varargin] = parsem( '1_exponential', varargin ); if( val );        for i = 0:15; S{end+1,2} = 2; S{end,1} = [1 exp(2^(-i-1))]'; S{end,end} = ['1_exponential (' num2str(i) ')']; end; end; %blf(1:15,'1_exp','nocheck','start',[exp(0) exp(1) -exp(2)]'); %nice example %asymptotic sum-rules
        [val,varargin] = parsem( '1_cantor', varargin ); if( val );             S{end+1,2} = 3; S{end,3} = [0 2]; S{end,end} = '1_cantor'; end;   %Tile: Cantor set %no full masks sets
        
        %no sum-rules
        [val,varargin] = parsem( '1_cex_uniquesol', varargin ); if( val );  %#ok<ALIGN>
            %two non-zero limit functions. the second consists again of two functions, which converge in the C_1 Norm. 
            %c=blf(15,'1_cex_uniquesol'); plot(c(1:2:end)); plot(c(2:2:end)); 
            %plot(diff(c(1:2:end),1)); plot(diff(c(2:2:end),1)); 
            %plot(diff(c(2:4:end),1)); plot(diff(c(4:4:end),1));             
            S{end+1,2} = 2;  S{end,1} = 2/12*[1 1 2 1 2 1 2 1 1]'; S{end,3} = [0 1]; S{end,end} = '1_cex_uniquesol';
        end;   
        
        %=====================================================
        %Dimension 2
        %=====================================================
        [~,varargin] = parsem( 'parse_one', varargin );   %remove parse_one
        [val,varargin] = parsem( '2_all', varargin ); if( val ); [~,varargin] = parsem( 'parse_one',varargin,1,'set'); end;
        [val,varargin] = parsem( '2_rand', varargin ); if( val );   S(end+1,:)=randS(2); end;
        
        [val,varargin] = parsem( '2_exlang', varargin ); if( val ); S{end+1,2} = [2 0; 0 3]; S{end,1} = 1/6*[1 1 3 2 1; 2 4 6 4 2; 1 3 3 2 1]; S{end,end} = '2_exlang'; end;
        
        [val,varargin] = parsem( {'2_sierp','2_sierpinski'}, varargin ); if( val );  S{end+1,2} = [2 0; 0 2]; S{end,1} = [0 1 0;1 1 0;0 0 1]; S{end,3} = [0 0; 1 0; 0 1; 1 1]';  S{end,end} = '2_sierp'; end; %sum-rules, but irreducible symbol %blf is Sierpinski Triangle
        [val,varargin] = parsem( '2_rqj43', varargin ); if( val );  S{end+1,2} = [1 -1; 1 1]; S{end,1} = 1/32*[0 -1 0 -1 0; -1 0 10 0 -1; 0 10 32 10 0; -1 0 10 0 -1; 0 -1 0 -1 0]; S{end,end} = '2_rqj43'; end;%Rong Qing Jia p656, Example 4.3
        [val,varargin] = parsem( '2_maria', varargin ); if( val );  S{end+1,2} = [1 3; 1 0]; S{end,3} = [0 0; 1 0; 2 0]'; S{end,end} = '2_maria'; end; %Tile: von Example von Maria
        [val,varargin] = parsem( '2_M2102', varargin ); if( val );  S{end+1,2} = [2  1;  0 2]; S{end,1} = 1/4*[0 0 0 0; -2 -1 4 3; -4 -2 8 6; -2 -1 4 3]; S{end,end} = '2_M2102'; end; %support phi: for example in appendix
        [val,varargin] = parsem( '2_not_jointly_expanding', varargin ); 
        if( val );
            S{end+1,2} = [2 -1; 2 0]; S{end,1} = 1/6*[0 1;0 3;2 2;3 0;1 0]; S{end,end} = '2_not_jointly_expanding (1)'; 
            S{end+1,2} = [0 -1;2 2];  S{end,1} = 1/4*[1 1 0 0;0 2 2 0;0 0 1 1];  S{end,end} = '2_not_jointly_expanding (2)'; 
        end;  %two matrices which are not jointly expanding
        [val,varargin] = parsem( '2_2_times_2', varargin ); if( val );      S{end+1,2} = [2  0;  0  2]; S{end,3} = [0 0; 1 0; 0 1; 1 1]'; S{end,end} = '2_2_times_2'; end;   %Tile: 2x2 
        [val,varargin] = parsem( '2_3_times_3', varargin ); if( val );      S{end+1,2} = [3  0;  0  3]; S{end,3} = [0 0; 1 0; 2 0; 0 1; 1 1; 2 1 ; 0 2; 1 2; 2 2;]'; S{end,end} = '2_3_times_3'; end;  %Tile: 3x3
        [val,varargin] = parsem( '2_Cha2017', varargin );   if( val );      S{end+1,2} = [2 1; 1 -1]; S{end,3} = [0 0; 1 0; 2 0]'; S{end,1} = 1/4*[0 2 0; 2 4 1; 0 3 0]; S{end,end} = '2_Cha2017'; end; %Example from Cha2017
        [val,varargin] = parsem( '2_frayed_squares', varargin ); if( val ); i=1; S{end+1,2} = [3 0; 0 3]; S{end,3} = [0 0; 1 0; 2 0; 0 1; 0 2; 1 2; 2 1; 2 2; 3*i+1 3*i+1]'; S{end,end} = '2_frayed_squares';   end; %Ausgefranste Quadrate, unbounded
        [val,varargin] = parsem( '2_frayed_squares_5', varargin ); if( val );    i=5; S{end+1,2} = [3 0; 0 3]; S{end,3} = [0 0; 1 0; 2 0; 0 1; 0 2; 1 2; 2 1; 2 2; 3*i+1 3*i+1]'; S{end,end} = '2_frayed_squares_5'; end; %Ausgefranste Quadrate, unbounded
        [val]=parsem( '2_frayed_squares_n', varargin ); 
        if( val );    
            [n,varargin] = parsem( '2_frayed_squares_n',varargin,0); 
            for i=1:n; 
                S{end+1,2} = [3 0; 0 3]; S{end,3} = [0 0; 1 0; 2 0; 0 1; 0 2; 1 2; 2 1; 2 2; 3*i+1 3*i+1]'; S{end,end} = ['2_frayed_squares_n (' num2str(i) ')']; 
            end;  
        end; %Ausgefranste Quadrate, unbounded
        [val,varargin] = parsem( '2_unbd_dilation', varargin ); 
        if( val ); 
            %Example for unbounded matrices with bounded digits
            %============================================================
            %Diese Matrizen sind nicht jointly expanding. Alle Matrizen haben DET=2, und das gleiche Digit set, Eigenwerte fuer alle Matrizen: I/Sqrt[DET], 
            %Es scheint als wuerden sie zumindest in der Reihenfolge M1,M2,M3,M4,... expanding sein, %fuer DET=2 
            %Attractor is immer unabhaengig von DET
            %Digit set sollte fuer DET=2 stimmen, fuer hoehere DET wahrhscheinlich
            %auskommentierter Code dient zum Testen der Matrizen
                  MAX = 10; 
                  DET = 2; 
                  for k = 0:MAX; 
                      %k=k*0.1;
                      S{end+1,2} = [k -k^2-DET; 1 -k]; 
                      if(DET>0); 
                          S{end,3} = [0:-1:-DET+1;zeros(1,DET)]; 
                      else; 
                          S{end,3} = [0:-DET-1;zeros(1,-DET)]; end;
                      %S{end,1} = 1/2*[1 1; 1 1];
                      S{end,end} = ['2_unbd_dilation (' num2str(k) ')'];
                  end; %
            %for k=1:MAX; M{k}=inv(M{k}); end;
            %BEGIN=4; PROD=eye(2);
            %for k=BEGIN:MAX; PROD=PROD*M{k}; end;
            %max(abs(eig(PROD)))
            %clear MAX DET BEGIN PROD;
        end;
        
        [val,varargin] = parsem( '2_loop', varargin ); if( val );                S{end+1,2} = [1 2; -2 -1];  S{end,3} = [0 0; 1 0; 1 -1]'; S{end,end} = '2_loop';  end;  %Tile: loop
        [val,varargin] = parsem( '2_butterfly', varargin ); if( val );           S{end+1,2} = [2 0; 0 2]; S{end,1} = 1/16*[0 0 0 0 -1 -1 0;0 0 -1 0 2 0 -1;0 -1 2 8 8 2 -1;0 0 8 16 8 0 0;-1 2 8 8 2 -1 0; -1 0 2 0 -1 0 0;0 -1 -1 0 0 0 0];  S{end,end} = '2_butterfly'; end;%butterfly
        [val,varargin] = parsem( {'2_twindragon','2_dragon','2_doubledragon'}, varargin); 
                                                        if( val );               S{end+1,2} = [1  1; -1  1]; S{end,3} = [0 0; 1 0]';  S{end,end} = '2_twindragon'; end;  %Tile: TwinDragon
        [val,varargin] = parsem( '2_flash', varargin ); if( val );               S{end+1,2} = [1  3;  1  0]; S{end,3} = [0 0; 1 0; 2 0]';  S{end,end} = '2_flash'; end;   %Tile: Flash 
        [val,varargin] = parsem( '2_flflash_ru', varargin ); if( val );          S{end+1,2} = [1  2;  1 -1]; S{end,3} = [0 0; 1 0; 2 0]';  S{end,end} = '2_flflash_ru';end;   %Tile: fractalless Flash/Paralellogram: right up
        [val,varargin] = parsem( '2_fltwindragon_ru', varargin ); if( val );     S{end+1,2} = [1  1;  1 -1]; S{end,3} = [0 0; 1 0]';  S{end,end} = '2_fltwindragon_ru'; end;   %Tile: fractalless TwinDragon/Paralellogram: right up
        [val,varargin] = parsem( '2_fldragon_rd', varargin ); if( val );         S{end+1,2} = [1  1;  1 -1]; S{end,3} = [0 0; 0 1]';  S{end,end} = '2_fldragon_ru'; end;   %Tile: Paralellogram: right down
        [val,varargin] = parsem( '2_fldragon_ld', varargin ); if( val );         S{end+1,2} = [1  1;  1 -1]; S{end,3} = [0 0; -1 0]'; S{end,end} = '2_fldragon_ld'; end;   %Tile: Paralellogram: left down 
        [val,varargin] = parsem( '2_spline_binary_0', varargin ); if( val );     S{end+1,2} = [2 0; 0 2]; S{end,1} = [1 1;1 1];  S{end,end} = '2_spline_binary_0'; end;
        [val,varargin] = parsem( '2_spline_binary_1', varargin ); if( val );     S{end+1,2} = [2 0; 0 2]; S{end,1} = 1/2*[0 1 1;1 2 1; 1 1 0];  S{end,end} = '2_spline_binary_1'; end;
        [val,varargin] = parsem( '2_spline_binary_2', varargin ); if( val );     S{end+1,2} = [2 0; 0 2]; S{end,1} = 1/8*[0 1 2 1; 1 4 5 2; 2 5 4 1; 1 2 1 0];  S{end,end} = '2_spline_binary_2'; end;
        [val,varargin] = parsem( '2_spline_binary_3', varargin ); if( val );     S{end+1,2} = [2 0; 0 2]; S{end,1} = 1/16*[0 0 1 2 1;0 2 6 6 2;1 6 10 6 1;2 6 6 2 0;1 2 1 0 0];  S{end,end} = '2_spline_binary_3'; end;
        [val,varargin] = parsem( '2_spline_binary_4', varargin ); if( val );     S{end+1,2} = [2 0; 0 2]; S{end,1} = 1/1024*[0 0 0 0 1 4 6 4 1;0 0 0 4 20 40 40 20 4;0 0 6 40 106 144 106 40 6;0 4 40 144 260 260 144 40 4;1 20 106 260 346 260 106 20 1;4 40 144 260 260 144 40 4 0;6 40 106 144 106 40 6 0 0;4 20 40 40 20 4 0 0 0;1 4 6 4 1 0 0 0 0];    S{end,end} = '2_spline_binary_4'; end;  
        [val,varargin] = parsem( '2_cex_tilearea', varargin ); if( val );        S{end+1,2} = [-4 -2;5 1]; S{end,3} = [-5 5;-4 4;-3 3;-2 2;-1 1;0 0]';   S{end,end} = '2_cex_tilearea'; end; %tile with standard digit set, but tile has area 2???
        
        %examples from papers
        %======================================
        [val,varargin] = parsem( '2_V0neqV0bar_1', varargin ); %convergent scheme, where V0 ~= V0bar
        if( val ); 
            S{end+1,2} = [-3 -4; 4 4];
            S{end,1} = [0 0 0 1 0; 1/3 0 1 0 1/3; 0 1/2 0 0 0; 0 1/2 0 0 1/3];
            S{end,3} = [-3 3; -2 2 ; -1 1 ; 0 0]';
            S{end,end} = '2_V0neqV0bar_1';
        end;
        
        [val,varargin] = parsem( '2_superomega', varargin );%example from the thesis     
        if( val ); 
            S{end+1,2} = [1 -1; 1 1];
            S{end,1} = [0 1/4 0;1/4 1 1/4; 0 1/4 0]; S{end,3} = [0 0; 0 1]';             
            S{end,end} = '2_superomega (1)';
            
            S{end+1,2} = [2 0; 0 2];            
            S{end,1} = [1/2 1/2 ;1 1; 1/2 1/2 ]; S{end,3} = [0 0; 0 1; 1 0; 1 1]'; 
            S{end,end} = '2_superomega (2)';
        end;

        [val,varargin] = parsem( '2_CDRT_ex531', varargin );
        if( val ); 
            S{end+1,2} = [2 0; 0 3]; 
            S{end,1} = [1/6 1/3 1/2 1/3 1/6; 1/3 2/3 1 2/3 1/3; 1/6 1/3 1/2 1/3 1/6]; 
            S{end,end} = '2_CDRT_ex531'; 
        end;
        
        [val,varargin] = parsem( '2_CDRT_ex532', varargin );
        if( val );         
            S{end+1,2} = [2 0; 0 3]; 
            S{end,1} = [0 0 0 -1/48 -1/24 -1/16 -1/24 -1/48 0 0 0; 0 0 0 0 0 0 0 0 0 0 0;-2/81 -5/162 0 89/432 89/216 9/16 89/216 89/432 0 -5/162 -2/81;-4/81 -5/81 0 10/27 20/27 1 20/27 10/27 0 -5/81 -4/81;-2/81 -5/162 0 89/432 89/216 9/16 89/216 89/432 0 -5/162 -2/81;0 0 0 0 0 0 0 0 0 0 0;0 0 0 -1/48 -1/24 -1/16 -1/24 -1/48 0 0 0]; 
            S{end,end} = '2_CDRT_ex532'; 
        end;
        
        [val,varargin] = parsem( '2_CDRT_ex533', varargin );
        if( val );         
        	S{end+1,2} = [2 0; 0 3]; 
            AA= [0 0 0 0 0 0 1/256 1/128; 0 0 0 0 0 0 0 0;0 0 0 1/324 5/1296 0 -241/6912 -241/3456;0 0 0 0 0 0 0 0;7/1458 4/729 0 -121/2916 -605/11664 0 20809/93312 20809/46656];
            BB= [7/729 8/729 0 -56/729 -70/729 0 280/729 560/729];
            CC= [3/256 0 -25/256 0 75/128 1 75/128 0 -25/256 0 3/256]';
            AA= [AA; BB; flipud(AA)];
            AA= [AA CC fliplr(AA)];
            S{end,1} = AA;
            S{end,end} = '2_CDRT_ex533'; 
        end;    
        
        [val,varargin] = parsem( '2_CDRT_ex541', varargin );
        if( val );         
            S{end+1,2} = [2 0; 0 5]; 
            S{end,1} = 1/10*[1 2 3 4 5 4 3 2 1; 2 4 6 8 10 8 6 4 2; 1 2 3 4 5 4 3 2 1];  
            S{end,end} = '2_CDRT_ex541'; 
        end;
        
        [val,varargin] = parsem( '2_CGRS', varargin );
        if( val );         
            S{end+1,2} = [1 1; 1 -2];  S{end,1} = sequence(1/3*[1 2 3 2 1],[0; -2]); S{end,end} = '2_CGRS (1)'; 
            S{end+1,2} = [2 -1; 1 -2]; S{end,1} = sequence(1/3*[1 2 3 2 1],[0; -2]); S{end,end} = '2_CGRS (2)';             
        end;
        
        [val,varargin] = parsem( '2_cool_1', varargin ); if( val );     S{end+1,2} = [-1 -1;4 -1];  S{end,end} = '2_cool_1'; end; %cool looking stuff
        [val,varargin] = parsem( '2_identity', varargin ); if( val );   S{end+1,2} = eye(2); S{end,end} = '2_identity';end; %no dilation matrices
        [val,varargin] = parsem( '2_cantor', varargin ); if( val );     S{end+1,2} = [3 0; 0 3]; S{end,3} = [0 0 ; 2 0; 0 2; 2 2]'; end;   %Tile: Cantor set %no full masks sets
        [val,varargin] = parsem( '2_McLure', varargin ); if( val );     S{end+1,2} = [1 -sqrt(3); sqrt(3) 1]; S{end,3} = [0 0; 1 0; -1/2 sqrt(3)/2; -1/2 -sqrt(3)/2]';  S{end,end} = '2_McLure';      end; %%non-integer scheme, self affine 4-tile, Mark McLure - Self affine tiles, p1
        [val,varargin] = parsem( '2_tridragon', varargin ); if( val );  S{end+1,2} = [3/2 -sqrt(3)/2; sqrt(3)/2 3/2]; S{end,3} = [0 0; 1 0; 1/2 sqrt(3)/2]';  S{end,end} = '2_tridragon';  end;   %%non-integer scheme, Tridragon for hexagonal lattice
        
        
        %=====================================================
        %Dimension 3
        %=====================================================
        [~,varargin] = parsem( 'parse_one', varargin );  %remove parse_one
        [val,varargin] = parsem( '3_all', varargin ); if( val ); [~,varargin] = parsem( 'parse_one',varargin,1,'set'); end;
        [val,varargin] = parsem( '3_rand', varargin ); if( val ); S(end+1,:) = randS( 3 ); end;
        [val,varargin] = parsem( '3_spline', varargin ); if( val ); val= convm([1 1],[1 1],[1 1],'outer'); val = convm( val, val ); S{end+1,1}=val; S{end,2} = 2*eye( 3 ); S{end,end} = '3_spline'; end;
        
        
 
        %=====================================================
        %Dimension 4
        %=====================================================
        [~,varargin] = parsem( 'parse_one', varargin );  %remove parse_one
        [val,varargin] = parsem( '4_all', varargin ); if( val ); [~,varargin] = parsem( 'parse_one', varargin, 1, 'set' ); end;
        [val,varargin] = parsem( '4_rand', varargin ); if( val ); S(end+1,:) = randS(4); end;
        [val,varargin] = parsem( '4_cex_Pot97', varargin ); if( val ); S{end+1,2} = [0 1 0 0; 0 0 1 0; 0 0 -1 2;-1 0 -1 1];  S{end,end} = '4_cex_Pot97'; end; %dilation matrix which has no digit set which generates a tile [Pot97]

end

function Srand = randS( dim )
    s = warning( 'off', 'normalizeS:support' );
    MM = 0;
    while( true );
        while( min(abs(eig(MM)))<=1.1 || abs(det(MM))>6 || abs(det(MM))<=1.5 ); 
            MM = randi( 6, dim )-3; end; 
        aa = randi( 5, [5*ones(1,dim) 1] )-1;

        Srand = normalizeS( {aa,MM,[],[num2str(dim) '_rand']}, 'verbose',0 );             
        if( ~isempty(Srand) ); 
            break; end; end;     
    warning(s);
end

function [S,varargin] = getS_namevalue( S, varargin )
    [S{end+1,1},varargin] = parsem( {'mask','a'}, varargin, [] );      
    [S{end,2},varargin] = parsem( {'dilation','M'}, varargin, [] );  
    [S{end,3},varargin] = parsem( {'digit','D'}, varargin, [] );     
    [S{end,end},varargin] = parsem( {'n','name'}, varargin,[]);      
end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 