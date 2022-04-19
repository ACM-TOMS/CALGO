function [ X, typeout ] = sph2cartm( sco, type )
% [ X, typeout ] = sph2cartm( sco, [type] )
% Transforms hypershperical to cartesian coordinates.
%
% See also: cart2sphm
%
% Written by: tommsch, 2018

%parse input
dim = size( sco, 1 );

if( nargin<=1 || isempty(type) );
    type=0; end;   

if( ischar(type) )
    switch type
        case {'hsph','h',''};       type =    0;
        case {'hsphdeg','hd'};      type =   -1;
        case {'azel','ae'};         type =    2;
        case {'azeldeg','aed'};     type =   -2;
        case {'azinc','ai'};        type =    3;
        case {'azincdeg','aid'};    type =   -3;
        case {'ant'};               type =    4;
        case {'antdeg','antd'};     type =   -4;
        case {'hcyl','hc'};         type =  100;
        case {'hcyldeg','hcd'};     type = -100;
        case {'hant'};              type =  104;
        case {'hantdeg','hantd'};   type = -104;
        otherwise; error( 'Wrong string for type.' ); end; end;
if( type==9 || type==10 || type==11 || type==12 )
    fprintf( 'Old type used. See the help of the function cart2sphm for the new numbers of the types.\nFix calling arguments.\n' ); 
    error( 'Error in sph2cartm' ); end;    
 
typeout = type;


%handle abnormal cases
if( nargin==0 || isempty(sco) ); 
    X=[]; 
    return; end;

%check input
switch type
    case {0,1,-1};
        %do nothing
    case {2,3,4,-2,-3,-4};
        if( dim<2 || dim>3);
            error( '''Type +-2, +-3, +-4'' only implemented for dim==2,3 yet.' ); end;       
    case {100,101,-100,-101}
        if( dim==1 );
            error( '''Type +-100, +-101'' (Hyperspherical coordinates) do not make sense for dim==1.' ); end;
    case {104,-104};
        if( dim~=3 );
            error( '''Type +-104'' only implement for dim=3 yet.' ); end;        
    case {105,-105};
        if( dim==2 );
            error( '''Type +-105'' (Spherinder coordinates) do not make sense for dim==2.' ); end;
    otherwise
        error( 'Unkown type.' ); end;

%transform coordinates to radians
switch type
    case {-1,-2,-3,-4};
        sco(1:end-1,:) = sco(1:end-1,:)/180*pi;
        type = -type;
    case {-100,-101,-104};
        sco(1,:) = sco(1,:)/180*pi;
        type = -type;
    case {-105};
        if( dim==1 );
            %do nothing
        else %dim>2
            sco(1:end-2,:) = sco(1:end-2,:)/180*pi; end;
    otherwise
        end; %fallthrough;
        

%pre-process coordinates
switch type
    case 3;
        %transforms type 3 to type 2
        if( dim==3 ) 
            sco(2,:) = pi/2-sco(2,:); end;
        type = 2;    
    case 4;
        %transforms type 4 to type 2
        sco(1,:) = pi/2-sco(1,:);
        if( dim==3); 
            sco(2,:) = -sco(2,:); end;
        type = 2;
    case 104;
        %transforms type 104 to type 100
        sco(1,:) = pi/2-sco(1,:);
        type = 100;
    otherwise;
        end; %fallthrough

%transform coordinates
switch type
    case {0,1};
        X = zeros( dim,size(sco,2) );
        r = sco(end,:);
        sinval = sin(sco);
        cosval = cos(sco);
        cosval(end,:) = 1;
        for i=1:dim
            X(i,:) = prod( [sinval(1:i-1,:); cosval(i,:)], 1 ); end;
        X=X.*r;

    case {2};
        switch dim;
            case 2;
                [x,y] = pol2cart( sco(1,:), sco(2,:) );
                X = [x;y];
            case 3;
                [x,y,z] = sph2cart( sco(1,:), sco(2,:), sco(3,:) );
                X = [x;y;z]; 
            otherwise
                error('Fatal error. Check Code.'); end;
            
    case {100, 101};
        [x,y] = pol2cart( sco(1,:), sco(end,:) );
        X = [x;y; sco(2:end-1,:)];
        
    case {105};
        error( 'Case +-105 not implemented yet.' );
        
    otherwise
        error('Fatal error. Check Code.'); end;

end

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing.   