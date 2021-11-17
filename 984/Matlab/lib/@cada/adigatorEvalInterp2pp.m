function zi = adigatorEvalInterp2pp(pp,xi,yi)
% CADA overloaded version of adigatorEvalInterp2pp
% 
% Syntax: ZI = adigatorEvalInterp2pp(pp,XI,YI)
%
% This function is used to evaluate the 2-D piecewise polynomials generated
% using the adigatorGenInterp2pp function. Since MATLAB does not have a
% built in 2-D polynomial evaluation command, this file will need to be in
% the user's path if they create a derivative file by using the adigator
% command on a function file which contains an interp2 command. The inputs
% to this are the pp generated by pp = adigatorGenInterp2pp(X,Y,Z,method),
% as well as the inputs XI, YI. 
%
% The following two lines of code:
%
%     pp = adigatorGenInterp2pp(X,Y,Z,method);
%     ZI = adigatorEvalInterp2pp(pp,XI,YI);
%
% Are (roughly) equivalent to the single line of code:
%
%     ZI = interp2(X,Y,Z,XI,YI,method,'extrap');
%
% Where the difference is due to rounding errors since when using the 
% adigatorGenInterp2pp and adigatorEvalInterp2pp commands, the 2-D
% coefficients are calculated, whereas MATLAB interp2 essentially performs
% interp1 twice, thus generating 1-D coefficients twice.
%
% Copyright 2011-2014 Matthew J. Weinstein and Anil V. Rao
% Distributed under the GNU General Public License version 3.0
%
% See also interp2, interp1, ppval, adigatorGenInterp2pp
global ADIGATOR
NUMvod  = ADIGATOR.NVAROFDIFF;
fid     = ADIGATOR.PRINT.FID;
PFLAG   = ADIGATOR.PRINT.FLAG;
indent  = ADIGATOR.PRINT.INDENT;
NDstr   = sprintf('%1.0f',ADIGATOR.DERNUMBER);

% -------------------------- Get pp Data -------------------------------- %

if isa(pp,'cada')
  if ADIGATOR.EMPTYFLAG
    zi = cadaEmptyEval(pp,xi,yi);
    return
  end
  ADIGATOR.VARINFO.LASTOCC(pp.id,1) = ADIGATOR.VARINFO.COUNT;
  ppknown = pp.func.pp;
  fpp = pp.func.value;
  fppStr = pp.func.name;
elseif isstruct(pp)
  if ADIGATOR.EMPTYFLAG
    zi = cadaEmptyEval(pp.xorder,pp.yorder,pp.coefs,pp.xbreaks,pp.ybreaks,xi,yi);
    return
  end
  ppknown = 1;
  fpp = pp;
  if ADIGATOR.FORINFO.FLAG && ~ADIGATOR.OPTIONS.UNROLL
    OverID = nonzeros(ADIGATOR.VARINFO.OVERMAP.FOR(fpp.xbreaks.id,:));
    if length(OverID) > 1; OverID = OverID(1); end
    if ~PFLAG
      if ~isempty(OverID)
        ADIGATOR.VARINFO.OVERMAP.FOR(ADIGATOR.VARINFO.OVERMAP.FOR == OverID) = 0;
      end
      OverID = nonzeros(ADIGATOR.VARINFO.OVERMAP.FOR(fpp.ybreaks.id,:));
      if ~isempty(OverID)
        if length(OverID) > 1; OverID = OverID(1); end
        ADIGATOR.VARINFO.OVERMAP.FOR(ADIGATOR.VARINFO.OVERMAP.FOR == OverID) = 0;
      end
        for I = 1:numel(fpp.coefs)
          OverID = nonzeros(ADIGATOR.VARINFO.OVERMAP.FOR(fpp.coefs{I}.id,:));
          if ~isempty(OverID)
            if length(OverID) > 1; OverID = OverID(1); end
            ADIGATOR.VARINFO.OVERMAP.FOR(ADIGATOR.VARINFO.OVERMAP.FOR == OverID) = 0;
          end
        end
    elseif PFLAG && isempty(OverID)
      ppknown = 0;
    end
  end
  fpp.xorder = fpp.xorder.func.value;
  fpp.yorder = fpp.yorder.func.value;
  if (isempty(fpp.xorder) && prod(pp.xorder.func.size)>1) || ...
    (isempty(fpp.yorder) && prod(pp.yorder.func.size)>1)
    error(['If using adigatorEvalInterp2pp within a loop which you wish to keep rolled, ',...
            'then the pp xorder and yorder must stay consistent across all calls']);
  end
  fpp.xbreaks = fpp.xbreaks.func.value;
  if isempty(fpp.xbreaks) && prod(pp.xbreaks.func.size) > 1; ppknown = 0; end
  fpp.ybreaks = fpp.ybreaks.func.value;
  if isempty(fpp.xbreaks) && prod(pp.xbreaks.func.size) > 1; ppknown = 0; end
  ADIGATOR.VARINFO.LASTOCC([pp.xbreaks.id pp.ybreaks.id pp.xorder.id pp.yorder.id]) = ADIGATOR.VARINFO.COUNT;
  for I = 1:numel(fpp.coefs)
    fpp.coefs{I} = fpp.coefs{I}.func.value;
    if isempty(fpp.coefs{I}) && prod(pp.coefs{I}.func.size) > 1; ppknown = 0; end
    ADIGATOR.VARINFO.LASTOCC(pp.coefs{I}.id,1) = ADIGATOR.VARINFO.COUNT;
  end
  if ~ppknown && ~(ADIGATOR.FORINFO.FLAG && ~ADIGATOR.OPTIONS.UNROLL && PFLAG)
    error('all values of adigatorInterp2 polys must be known numerically')
  end
  if PFLAG
    nameloc = strfind(pp.xbreaks.func.name,'.xbreaks');
    fppStr = pp.xbreaks.func.name(1:nameloc(end)-1);
  end
else
  error('pp input must be generated by adigatorGenInterp2pp')
end

xorder = fpp.xorder;
yorder = fpp.yorder;

% ------------------------ Parse XI and YI ------------------------------ %
[xi, xiMrow, xiNcol] = parseinput(xi);
[yi, yiMrow, yiNcol] = parseinput(yi);

if xi.func.size == yi.func.size
  % X, Y, Z same size
  ziMrow = xi.func.size(1); ziNcol = xi.func.size(2);
  repflag = 0;
elseif (xiMrow == 1 && yiNcol == 1) || (xiNcol == 1 && yiMrow == 1)
  % Z takes on row dimension of length of Y and col dimension of length of
  % X
  ziMrow = yiMrow*yiNcol;   ziNcol = xiMrow*xiNcol;
  if isinf(ziMrow) && isinf(ziNcol)
    error('Result of vectorized interp2 may only have 1 vectorized dimension')
  elseif isinf(ziMrow)
    yrep = ones(1,ziNcol);
    if cadaCheckForDerivs(xi)
      error('Cannot perform interp2 when YI vectorized, XI non-vectorized and has derivative information')
    end
  elseif isinf(ziNcol)
    xrep = ones(ziMrow,1);
    if cadaCheckForDerivs(yi)
      error('Cannot perform interp2 when XI vectorized, YI non-vectorized and has derivative information')
    end
  else
    xrep = repmat(1:ziNcol,[ziMrow 1]);     xrep = xrep(:);
    yrep = repmat((1:ziMrow).',[1 ziNcol]); yrep = yrep(:);
  end
  repflag = 1; 
else
  error('XI and YI must be the same size or vectors of different orientations.')
end

% Build ZI func
zi.id = ADIGATOR.VARINFO.COUNT;
[funcstr,DPFLAG] = cadafuncname();
zi.func = struct('name',funcstr,'size',[ziMrow ziNcol],'zerolocs',...
  [],'value',[]);
if isempty(xi.func.value) && isempty(yi.func.value) && ppknown
  zi.func.value = adigatorEvalInterp2pp(fpp,xi.func.value,yi.func.value);  
end
if isinf(ziMrow); 
  ziMrow = 1; zvec = 1; 
elseif isinf(ziNcol); 
  ziNcol = 1; zvec = 2;
else
  zvec = 0;
end
% Check Logical Reference
if repflag
  if isfield(xi.func,'logicref') || isfield(yi.func,'logicref')
    zi.func.logicref = zeros(1,2);
    if isfield(xi.func,'logicref')
      zi.func.logicref(2) = nonzeros(xi.func.logicref);
    end
    if isfield(yi.func,'logicref')
      zi.func.logicref(1) = nonzeros(yi.func.logicref);
    end
  end
elseif isfield(xi.func,'logicref') && isfield(yi.func,'logicref') && ...
    isequal(xi.func.logicref,yi.func.logicref)
  zi.func.logicref = xi.func.logicref;
elseif isfield(xi.func,'logicref') || isfield(yi.func,'logicref')
  error('Invalid binary operation on result of an unknown logical reference')
end

% Build ZI deriv
ppxflag = 0; ppyflag = 0;
zi.deriv = struct('name',cell(NUMvod,1),'nzlocs',cell(NUMvod,1));
for Vcount = 1:NUMvod
  nv = ADIGATOR.VAROFDIFF(Vcount).usize;
  if xorder < 2; xi.deriv(Vcount).nzlocs = []; end
  if yorder < 2; xi.deriv(Vcount).nzlocs = []; end
  % Repmat derivs if needed
  if repflag && ~isempty(xi.deriv(Vcount).nzlocs)
    xrows = xi.deriv(Vcount).nzlocs(:,1);
    xcols = xi.deriv(Vcount).nzlocs(:,2);
    dx = sparse(xrows,xcols,1:length(xrows),ziNcol,nv);
    dx = dx(xrep,:);
    [xrows, xcols,xrepind] = find(dx);
    if size(xrows,2) > 1; xrows = xrows.'; xcols = xcols.'; end
    xi.deriv(Vcount).nzlocs = [xrows, xcols];
    if DPFLAG
      TD1 = ['cada',NDstr,'td1'];
      TDind1 = cadaindprint(xrepind);
      if zvec
        fprintf(fid,[indent,TD1,' = ',xi.deriv(Vcount).name,'(:,',TDind1,');\n']);
      else
        fprintf(fid,[indent,TD1,' = ',xi.deriv(Vcount).name,'(',TDind1,');\n']);
      end
      xi.deriv(Vcount).name = TD1;
    end
  end
  if repflag && ~isempty(yi.deriv(Vcount).nzlocs)
    yrows = yi.deriv(Vcount).nzlocs(:,1);
    ycols = yi.deriv(Vcount).nzlocs(:,2);
    dy = sparse(yrows,ycols,1:length(yrows),ziMrow,nv);
    dy = dy(yrep,:);
    [yrows, ycols,yrepind] = find(dy);
    if size(yrows,2) > 1; yrows = yrows.'; ycols = ycols.'; end
    yi.deriv(Vcount).nzlocs = [yrows, ycols];
    if DPFLAG
      TD2 = ['cada',NDstr,'td2'];
      TDind1 = cadaindprint(yrepind);
      if zvec
        fprintf(fid,[indent,TD2,' = ',yi.deriv(Vcount).name,'(:,',TDind1,');\n']);
      else
        fprintf(fid,[indent,TD2,' = ',yi.deriv(Vcount).name,'(',TDind1,');\n']);
      end
      yi.deriv(Vcount).name = TD2;
    end
  end
  % Print out dZdX and dZdY if needed
  if DPFLAG && ~ppxflag && ~isempty(xi.deriv(Vcount).nzlocs)
    % Get pp which is for dzi/dxi
    if ppknown
      ppxStr     = getppderiv(fpp,'x');
    else
      ppxStr     = getppderiv(fpp,'x',fppStr,fid,indent,NDstr);
    end
    dZdXstr = ['cada',NDstr,'tf1'];
    fprintf(fid,[indent,dZdXstr,' = adigatorEvalInterp2pp(',ppxStr,',',xi.func.name,',',yi.func.name,');\n']);
  end
  if DPFLAG && ~ppyflag && ~isempty(yi.deriv(Vcount).nzlocs)
    if ppknown
      ppyStr     = getppderiv(fpp,'y');
    else
      ppyStr     = getppderiv(fpp,'y',fppStr,fid,indent,NDstr);
    end
    dZdYstr = ['cada',NDstr,'tf2'];
    fprintf(fid,[indent,dZdYstr,' = adigatorEvalInterp2pp(',ppyStr,',',xi.func.name,',',yi.func.name,');\n']);
  end
  if ~isempty(xi.deriv(Vcount).nzlocs) && ~isempty(yi.deriv(Vcount).nzlocs)
    % XI and YI have derivs
    derivstr = cadadername(funcstr,Vcount);
    zi.deriv(Vcount).name = derivstr;
    [zi.deriv(Vcount).nzlocs,dzxind,dzyind,xindflag,yindflag] = ...
      cadaunion(xi.deriv(Vcount).nzlocs,yi.deriv(Vcount).nzlocs,...
      ziMrow*ziNcol,ADIGATOR.VAROFDIFF(Vcount).usize);
    xrows = xi.deriv(Vcount).nzlocs(:,1);
    yrows = yi.deriv(Vcount).nzlocs(:,1);
    nzz   = size(zi.deriv(Vcount).nzlocs,1);
    if DPFLAG
      TD3 = ['cada',NDstr,'td3'];
      % Print dZdX
      TF3 = ['cada',NDstr,'tf3'];
      TFind = cadaindprint(xrows);
      if zvec == 1
        fprintf(fid,[indent,TF3,' = ',dZdXstr,'(:,',TFind,');\n']);
      elseif zvec == 2
        fprintf(fid,[indent,TF3,' = ',dZdXstr,'(',TFind,',:).'';\n']);
      else
        fprintf(fid,[indent,TF3,' = ',dZdXstr,'(',TFind,');\n']);
      end
      if xindflag
        DZXstr = TD3;
      elseif zvec
        fprintf(fid,[indent,TD3,' = zeros(size(',xi.deriv(Vcount).name,',1),%1.0f);\n'],nzz);
        TDind = cadaindprint(dzxind);
        DZXstr = [TD3,'(:,',TDind,')'];
      else
        fprintf(fid,[indent,TD3,' = zeros(%1.0f,1);\n'],nzz);
        TDind = cadaindprint(dzxind);
        DZXstr = [TD3,'(',TDind,')'];
      end
      if zvec
        fprintf(fid,[indent,DZXstr,' = ',TF3,'.*',xi.deriv(Vcount).name,';\n']);
      else
        fprintf(fid,[indent,DZXstr,' = ',TF3,'(:).*',xi.deriv(Vcount).name,';\n']);
      end
      
      % Print dZdY
      TFind = cadaindprint(yrows);
      if zvec == 1
        fprintf(fid,[indent,TF3,' = ',dZdYstr,'(:,',TFind,');\n']);
      elseif zvec == 2
        fprintf(fid,[indent,TF3,' = ',dZdYstr,'(',TFind,',:).'';\n']);
      else
        fprintf(fid,[indent,TF3,' = ',dZdYstr,'(',TFind,');\n']);
      end
      if yindflag
        DZYstr = TD3;
      elseif zvec
        TDind = cadaindprint(dzyind);
        DZYstr = [TD3,'(:,',TDind,')'];
      else
        TDind = cadaindprint(dzyind);
        DZYstr = [TD3,'(',TDind,')'];
      end
      if zvec
        fprintf(fid,[indent,DZYstr,' = ',DZYstr,' + ',TF3,'.*',yi.deriv(Vcount).name,';\n']);
      else
        fprintf(fid,[indent,DZYstr,' = ',DZYstr,' + ',TF3,'(:).*',yi.deriv(Vcount).name,';\n']);
      end
      fprintf(fid,[indent,derivstr,' = ',TD3,';\n']);
    end
  elseif ~isempty(xi.deriv(Vcount).nzlocs)
    % XI has derivs
    derivstr = cadadername(funcstr,Vcount);
    zi.deriv(Vcount).name = derivstr;
    zi.deriv(Vcount).nzlocs = xi.deriv(Vcount).nzlocs;
    if DPFLAG
      % Print dZdX
      xrows = xi.deriv(Vcount).nzlocs(:,1);
      TFind = cadaindprint(xrows);
      if zvec == 1
        fprintf(fid,[indent,derivstr,' = ',dZdXstr,'(:,',TFind,').*',xi.deriv(Vcount).name,';\n']);
      elseif zvec == 2
        fprintf(fid,[indent,derivstr,' = ',dZdXstr,'(',TFind,',:).''.*',xi.deriv(Vcount).name,';\n']);
      else
        TF3 = ['cada',NDstr,'tf3'];
        fprintf(fid,[indent,TF3,' = ',dZdXstr,'(',TFind,');\n']);
        fprintf(fid,[indent,derivstr,' = ',TF3,'(:).*',xi.deriv(Vcount).name,';\n']);
      end
    end
  elseif ~isempty(yi.deriv(Vcount).nzlocs)
    % YI has derivs
    derivstr = cadadername(funcstr,Vcount);
    zi.deriv(Vcount).name = derivstr;
    zi.deriv(Vcount).nzlocs = yi.deriv(Vcount).nzlocs;
    if DPFLAG
      % Print dZdY
      yrows = yi.deriv(Vcount).nzlocs(:,1);
      TFind = cadaindprint(yrows);
      if zvec == 1
        fprintf(fid,[indent,derivstr,' = ',dZdYstr,'(:,',TFind,').*',yi.deriv(Vcount).name,';\n']);
      elseif zvec == 2
        fprintf(fid,[indent,derivstr,' = ',dZdYstr,'(',TFind,',:).''.*',yi.deriv(Vcount).name,';\n']);
      else
        TF3 = ['cada',NDstr,'tf3'];
        fprintf(fid,[indent,TF3,' = ',dZdYstr,'(',TFind,');\n']);
        fprintf(fid,[indent,derivstr,' = ',TF3,'(:).*',yi.deriv(Vcount).name,';\n']);
      end
    end
  end
end

% Print Function
if PFLAG
  fprintf(fid,[indent,funcstr,' = adigatorEvalInterp2pp(',fppStr,',',...
    xi.func.name,',',yi.func.name,');\n']);
end
% if isstruct(pp)
ADIGATOR.VARINFO.LASTOCC([xi.id yi.id],1) = ADIGATOR.VARINFO.COUNT;
% else
%   ADIGATOR.VARINFO.LASTOCC([pp.id xi.id yi.id],1) = ADIGATOR.VARINFO.COUNT;
% end
ADIGATOR.VARINFO.COUNT = ADIGATOR.VARINFO.COUNT+1;
zi = cada(zi);
return
end

function [y, yMrow, yNcol] = parseinput(x)
if isa(x,'cada')
  y = x; yMrow = x.func.size(1); yNcol = x.func.size(2);
elseif isnumeric(x)
  [yMrow, yNcol] = size(x);
  y.id = [];
  y.func = struct('name',[],'size',[yMrow, yNcol],'zerolocs',[],'value',x);
else
  error('Invalid input to interp1')
end
end

function dppStr = getppderiv(pp,dim,varargin)
yorder    = pp.yorder;
xorder    = pp.xorder;
if nargin == 2
  switch dim
    case 'y'
      dyorder   = yorder-1;
      coefs     = pp.coefs;
      % Remove Constant Elements
      coefs(yorder,:) = [];
      % Multiply
      for I = 1:dyorder
        for J = 1:xorder
          coefs{I,J} = coefs{I,J}.*(yorder-I);
        end
      end
      % Assign dpp
      dpp = pp;
      dpp.coefs  = coefs;
      dpp.yorder = dyorder;
    case 'x'
      dxorder   = xorder-1;
      coefs     = pp.coefs;
      % Remove Constant Elements
      coefs(:,xorder) = [];
      % Multiply
      for J = 1:dxorder
        for I = 1:yorder
          coefs{I,J} = coefs{I,J}.*(xorder-J);
        end
      end
      % Assign dpp
      dpp = pp;
      dpp.coefs  = coefs;
      dpp.xorder = dxorder;
  end
  dppStr  = cadamatprint(dpp);
else
  fppStr = varargin{1}; fid = varargin{2}; 
  indent = varargin{3}; NDstr = varargin{4};
  switch dim
    case 'y'
      dyorder = yorder-1;
      dppStr = ['cada',NDstr,'dpp2y'];
      % Initialize and set order
      fprintf(fid,[indent,dppStr,'.form = ''adigatorpp2'';\n']);
      fprintf(fid,[indent,dppStr,'.xbreaks = ',fppStr,'.xbreaks;\n']);
      fprintf(fid,[indent,dppStr,'.ybreaks = ',fppStr,'.ybreaks;\n']);
      fprintf(fid,[indent,dppStr,'.yorder = %1.0f;\n'],dyorder);
      fprintf(fid,[indent,dppStr,'.xorder = %1.0f;\n'],xorder);
      % Remover Constant Elements
      fprintf(fid,[indent,dppStr,'.coefs = cell(%1.0f,%1.0f);\n'],dyorder,xorder);
      % Multiply
      for I = 1:dyorder
        for J = 1:xorder
          fprintf(fid,[indent,dppStr,'.coefs{%1.0d,%1.0d} = ',...
            fppStr,'.coefs{%1.0d,%1.0d}.*%1.0d;\n'],I,J,I,J,yorder-I);
        end
      end
    case 'x'
      dxorder = xorder-1;
      dppStr = ['cada',NDstr,'dpp2x'];
      % Initialize and set order
      fprintf(fid,[indent,dppStr,'.form = ''adigatorpp2'';\n']);
      fprintf(fid,[indent,dppStr,'.xbreaks = ',fppStr,'.xbreaks;\n']);
      fprintf(fid,[indent,dppStr,'.ybreaks = ',fppStr,'.ybreaks;\n']);
      fprintf(fid,[indent,dppStr,'.yorder = %1.0f;\n'],yorder);
      fprintf(fid,[indent,dppStr,'.xorder = %1.0f;\n'],dxorder);
      % Remover Constant Elements
      fprintf(fid,[indent,dppStr,'.coefs = cell(%1.0f,%1.0f);\n'],yorder,dxorder);
      % Multiply
      for J = 1:dxorder
        for I = 1:yorder
          fprintf(fid,[indent,dppStr,'.coefs{%1.0d,%1.0d} = ',...
            fppStr,'.coefs{%1.0d,%1.0d}.*%1.0d;\n'],I,J,I,J,xorder-J);
        end
      end
  end
end
end