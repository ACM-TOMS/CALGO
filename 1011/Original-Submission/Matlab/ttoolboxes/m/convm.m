function d=convm(varargin)
% d = convm( a1 ,a2, a3, ..., an, [options])
% d = ( ... ((a1 * a2) * a3) * ... * an )
% Consistent behaviour of conv, with regards to multi-dimensional applications.
%
% Input:
%   ai                  Array or size compatible cell-arrays
%                       
% Options:
%   'mult',val          vector of length equal to the number of ai's
%                       makes the conv (a1*..*a1 ) * (a2*..*a2) * ... * (an*...*an) each of which mult(i) often
%   'outer'             computes the outer or tensor product of arrays instead
%   'normalize'         divides through (product of all sums)
%
% E.g.: convm([1 1],[1 1],[1 1],'normalize')
%       convm([1 2 1],[1; 3; 1],'mult',[2 1])
%       convm([1 2 3],[2 3 4],'outer')
%
% See also conv, conv2, convn
%
% Written by: tommsch, 2017
% Uses code from: https://de.mathworks.com/matlabcentral/profile/authors/11-us

 %#ok<*ALIGN>

[mult,varargin]=parsem('mult',varargin,[]);
[outer,varargin]=parsem('outer',varargin);
[normalize,varargin]=parsem('normalize',varargin);
    
if(~isempty(mult));
    varargin=rude(mult, varargin);
end

na=size(varargin,2);

if(outer)
    for i=1:na; varargin{i}=squeezem(varargin{i}); end;
end

if(na==1); 
    d=varargin{1}; 
else
    if(outer)
        ndim=cumsum(cellfun(@ndimsm,varargin));
        for i=1:na; varargin{i}=permute(varargin{i}, circshift(1:ndim(end),ndim(1)-ndim(i))); end
    end
    
    
    d=convm_worker(varargin{1},varargin{2}); %initialize d
    for i=3:na
        d=convm_worker(d,varargin{i});
    end
end

if(normalize);
    d=d/prod(cellfun(@(x) sum(x(:)),varargin));    
end

end

function d=convm_worker(a, c)
    if(isempty(a)); 
        d=a; 
        return; end;
    if(isempty(c)); 
        d=c; 
        return; end;
    
    try
        issyma=issym(a);
        issymc=issym(c);
    catch
        issyma=0;
        issymc=0;
    end
    
    
    if(~issyma && ~issymc && ~anym(a) && ~anym(c)); 
        d=0; 
        return; end;
    

    dim=max(ndimsm(a),ndimsm(c));
    if(issyma || issymc)
        d=symbol2mask(mask2symbol(a,'dim',dim,'var','z1')*mask2symbol(c,'dim',dim,'var','z1'),'dim',dim,'sym','var','z1');
    elseif(dim==1);
        d=conv(a,c);
    elseif(dim==2);
        d=conv2(a,c);
    else
        d=convn(a,c);
    end
end

function	[p1,p2]=rude(varargin)
% [ vec ] = rude( len, val )
% Run-length DEcoding
%
% [ len, val ] = rude( vec )
% Run-length ENcoding
%
% P = rude
% Retrieve subroutine handles in structure P
%   P.d		for run-length DEcoding
%			vec = P.d( len, val )
%   P.e		for run-length ENcoding
%			[ len, val ] = P.e( vec )
%
% len   : repeat each val corresponding len times to create vec
% val	: 1xN template array of vals to be repeated len times
%           - numericals
%           - strings
%           - cells (any contents)
% vec   : DEcode = reconstruced output vector from len/val
%           ENcode = input vector to be encoded into len/val
%
% Notes:
% 	len <= 0 will remove corresponding VALs
% 	<NaN>s and <+Inf>s are treated as not-equal (expected behavior)
%
% E.g.: vec = rude([1 2 3 0 4],[10 inf nan pi 20])
%       vec = 10 Inf Inf NaN NaN NaN 20 20 20 20
%       [len,val] = rude(vec)
%       % len = 1 1 1 1 1 1 4 % note nan~=nan / inf~=inf!
%       % val = 10 Inf Inf NaN NaN NaN 20
%
%       s = rude;
%       v.x = pi;
%       w.x = pi; % note <v> and <w> are equal!
%       vec = s.d([1 0 3 2 2 2 3],{'a' 'b' 'cd' 1:3 v w magic(3)})
%       % vec = 'a' 'cd' 'cd' 'cd' 1x3D 1x3D 1x1S 1x1S 1x1S 1x1S 3x3D
%       [len,val] = s.e(vec)
%       % len = 1 3 2 4 3
%       % val = 'a' 'cd' 1x3D 1x1S 3x3D
%
% Written by : us, 18-Nov-2004
% Modified: us,	30-Nov-2004 20:42:03	/ TMW FEX
    p2=[];
	if( ~nargin & nargout )
        p1.d=@rl_decode;
        p1.e=@rl_encode;
	elseif	~nargin
        help(mfilename);
        return;
	else
        if( nargin == 1 )
            [p1,p2]=rl_encode(varargin{1});
        elseif( nargin >= 2 )
            p1=rl_decode(varargin{1:2});
        end
    end
end
%--------------------------------------------------------------------------------
function	vec=rl_decode(len,val)
% run-length decoder
    lx=len>0 & ~(len==inf);
	if( ~any(lx) );
		vec=[];
		return;
	end
	if( numel(len) ~= numel(val) );
		error('rl-decoder: length mismatch\nlen = %-1d\nval = %-1d',numel(len),numel(val));
	end
    len=len(lx);
    val=val(lx);
    val=val(:).';
    len=len(:);
    lc=cumsum(len);
    lx=zeros(1,lc(end));
    lx([1;lc(1:end-1)+1])=1;
    lc=cumsum(lx);
    vec=val(lc);
end
%--------------------------------------------------------------------------------
function	[len,val]=rl_encode(vec)
% run-length encoder
	switch	class(vec)
        case 'cell'; 
            [len,val]=rl_encodec(vec);
        case 'char'; 
            [len,val]=rl_encoden(double(vec)); val=char(val);
        otherwise; 
            [len,val]=rl_encoden(vec);
	end
end
%--------------------------------------------------------------------------------
function	[len,val]=rl_encoden(vec)
% run-length encode doubles
	if( isempty(vec) )
		len=0;
		val=[];
		return;	end;
    vec=vec(:).';
    vx=[1 diff(double(vec))];
    vx=vx~=0;
    val=vec(vx);
    vc=cumsum(vx).';
    len=accumarray(vc,ones(size(vc))).';
end
%--------------------------------------------------------------------------------
% run-length encode cells
function [len,val]=rl_encodec(vec)

    len=0;
    val={};
    vl=length(vec)+1;
    cl=cellfun('length',vec);
    tix=0;
    tlen=1;
    for	i=1:length(cl)
        if(	cl(i) )
            tmpl=vec{i};
            for	j=i+1:vl
                if	( j>length(cl) || ~isequalwithequalnans(tmpl,vec{j}) ) %#ok<DISEQN>
                    tix=tix+1;
                    val{tix}=tmpl; %#ok<AGROW>
                    len(tix)=tlen; %#ok<AGROW>
                    cl(i:j-1)=0;
                    tlen=1;
                    break;
                else
                    tlen=tlen+1;
                end
            end
        end
    end
end
%--------------------------------------------------------------------------------

function dummy; end %#ok<DEFNU> %Generates an error, if the 'end' of a function is missing. 