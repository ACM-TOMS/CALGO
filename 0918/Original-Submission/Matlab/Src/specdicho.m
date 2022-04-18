function [P, H] = specdicho(varargin)
%                                                                  %
%------------------------------------------------------------------%
% SPECTRAL DICHOTOMY OF A REGULAR MATRIX PENCIL WITH RESPECT TO    %
%     THE CIRCLE, ELLIPSE, IMAGINARY AXIS AND PARABOLA             %
%------------------------------------------------------------------%
% SPECDICHO computes:                                              %
% -- P: THE SPECTRAL PROJECTOR ONTO THE RIGHT DEFLATING SUBSPACE   %
%       OF (\lambda*B - A) CORRESPONDING TO THE EIGENVALUES INSIDE %
%       THE CIRCLE, ELLIPSE, PARABOLA AND WITH REAL PARTS IN THE   %
%       CASE OF IMAGINARY AXIS                                     %
% -- H: THE MATRIX INTEGRAL WHOSE NORM ||H||_2 INDICATES THE       %
%       NUMERICAL QUALITY OF P.                                    %
%       ||H||_2 IS CALLED THE DICHOTOMY CONDITION NUMBER.          %
%                                                                  %
%                                                                  %
% THE SPECTRAL DICHOTOMY OF (\lambda*B - A) WITH RESPECT TO THE    %
% CIRCLE OF CENTER C and RADIUS R                                  %
%                                                                  %
%                                                                  %
% THE SPECTRAL DICHOTOMY OF (\lambda*B - A) WITH RESPECT TO THE    %
% ELLIPSE CENTERED AT THE ORIGIN AND WITH SEMI-MAJOR AXIS a AND    %
% SEMI-MINOR AXIS b, IS GIVEN BY APPLYING THE UNIT CIRCULAR        %
% SPECTRAL DICHOTOMY TO (\lambda*BB - AA) WHERE                    %
% BB = [a1*B , -A ; 0 , a1*B] and  AA = [-b1*B , 0 ; A , -b1*B]    %
% a1 = (a+b)/2                and  b1 = (a-b)/2                    %
%                                                                  %
%                                                                  %
%                                                                  %
% THE SPECTRAL DICHOTOMY OF (\lambda*B - A) WITH RESPECT TO THE    %
% IMAGINARY AXIS IS OBTAINED BY APPLYING THE SPECTRAL DICHOTOMY TO %
% (\lambda*(A-B) - (A+B)) WITH RESPECT TO THE UNIT CIRCLE.         %
% WE ASSUME THAT THE MATRIX B IS NONSINGULAR.                      %
%                                                                  %
%                                                                  %
%                                                                  %
% THE SPECTRAL DICHOTOMY OF (\lambda*B - A) WITH RESPECT TO THE    %
% PARABOLA Y^2 = 2*P*(P/2 - X) IS OBTAINED BY APPLYING THE SPECTRAL%
% DICHOTOMY BY THE IMAGINARY AXIS TO (\lambda*BB - AA).            %
% WHERE                                                            %
% BB = [B , 0 ; 0 , I]  and                                        %
% AA = [-\SQRT{P/2}*B , A ; I , -\SQRT{P/2}*I]                     %
%                                                                  %
%------------------------------------------------------------------%
%------------------------------------------------------------------%
% OutPut Options:                                                  %
% 0. SPECDICHO(A) DISPLAYS THE SPECTRAL PROJECTOR P.               %
%    DEFAULT VALUE B = EYE(SIZE(A))                                %
%                                                                  %
% 1. SPECDICHO(A,B)  DISPLAYS THE SPECTRAL PROJECTOR P.            %
%    A MUST BE NUMERIC SQUARE AND DENSE                            %
%    B MUST BE NUMERIC SQUARE AND DENSE                            %
%    B MUST BE OF THE SAME SIZE AS A                               %
%                                                                  %
%                                                                  %
% 2. [P] = SPECDICHO(A,B) RETURNS THE MATRIX P                     %
%                                                                  %
% 3. [P,H] = SPECDICHO(A,B) RETURNS THE MATRICES P and H           %
%------------------------------------------------------------------%
%------------------------------------------------------------------%
% InPut Options:                                                   %
%     ... = SPECDICHO(A,OPTS) OR  ... = SPECDICHO(A,B,OPTS)        %
%  OPTS IS A STRUCTURE CONTAINING INPUT PARAMETERS.                %
%                                                                  %
%  THE STRUCTURE OPTS MAY CONTAIN SOME OR ALL OF THE FOLLOWING     %
%  INPUT PARAMETERS.                                               %
%  THE STRING FOR THE INPUT PARAMETERS CAN CONTAIN UPPER OR LOWER  %
%  CASE CHARACTERS.                                                %
%                                                                  %
%  IF THE INPUT ARGUMENT B IS OMITTED, THEN ITS DEFAULT VALUE IS   %
%  THE IDENTITY MATRIX B = EYE(SIZE(A))                            %
%                                                                  %
% INPUT PARAMETER:                                                 %
%                                                                  %
%    opts.geom   positively oriented contour in the complex plane. %
%                This allows the user to specify the geometry      %
%                of work. It must be equal to 'C': circle or       %
%                'E': ellipse or 'I': imaginary axis or            %
%                'P': parabola.                                    % 
%                Its default value is set to 'C'.                  %
%                                                                  %
%    opts.c      center of the circle.                             %
%                This  allows the user to specify the center       %
%                of the circle when geom is equal to 'C'.          %
%                It must be a complex number.                      %
%                Its default value is set to 0.                    %
%                                                                  %
%   opts.r       radius of the circle.                             %
%                This allows the user to specify the radius        %
%                of the circle when geom is equal to 'C'.          %
%                It must be a positive real number.                %
%                Its default value is set to  1.                   %
%                                                                  %
%   opts.a      semi-major axis of the ellipse                     %
%                   (X/a)^2 + (Y/b)^2 = 1 .                        %
%               This allows the user to specify the                %
%               semi-major axis of the ellipse when                %
%               geom is equal to 'E'.                              %
%               It must be a positive real number.                 %
%               Its default value is set to 5.                     %
%                                                                  %
%   opts.b      semi-minor axis of the ellipse                     %
%                   (X/a)^2 + (Y/b)^2 = 1 .                        %
%               This allows the user to specify the                %
%               semi-minor axis of the ellipse  when               %
%               geom is equal to 'E'.                              %
%               It must be a positive real number and  a >= b >0.  %
%               Its default value is set to 1.                     %
%                                                                  %
%   opts.p      positive real parameter of the parabola            %
%                    Y^2 = 2*p*(p/2 - X).                          %     
%               This allows the user to specify the                %
%               parameter of the parabola given by when            %
%               geom is equal to 'P'.                              %
%               It must be a positive real number.                 %
%               Its default value is set to 1.                     %
%                                                                  %
%  opts.mxiter maximal number of iterations to perform.            %
%              This allows the user to specify the maximal         %
%              number of iterations SPECDICHO will perform.        %
%              The default value is set to 10.                     %
%              The user can set a larger value for particularly    %
%              difficult problems.                                 %
%                                                                  %
%  opts.tol    tolerance used for convergence check.               %
%              This allows the user to specify the value of        %
%              tolerance.                                          %
%              Its default value is set to 1.0e-10.                %
%              It can be much smaller (e.g., tol = eps) for        % 
%              difficult problems.                                 %
%                                                                  %
%  opts.Ho     Hermitian positive definite matrix used             %
%              for scaling purposes.                               %
%              This option  must be of the same size as A          % 
%              when geom = ['C'|'I'] and of twice the order of A   %
%              when geom = ['E'|'P'].                              %
%              Its default value is set at eye(size(A))            %
%              when geom = ['C'|'I'] and set at eye(2*size(A))     %
%              when geom = ['E'|'P'].                              %
%                                                                  %
%------------------------------------------------------------------%
%------------------------------------------------------------------%
% DATE MODIFIED: 2011/06/05                                        %
% VERSION:  1.0                                                    %
%                                                                  %
% AUTHORS:                                                         %
% Miloud Sadkane  D\'epartement de Math\'ematiques,                %
%                 Universit\'e de Brest,                           %
%                 Brest, France.                                   %
%                 E-mail: Miloud.Sadkane@univ-brest.fr             %
%                                                                  %
% Ahmed Touhami   Mathematics and Computer Science Department      %
%                 Hassan 1st University,                           %
%                 Faculty of Sciences and Technologies,            %
%                 Settat, Morocco.                                 %
%                 E-mail: Ahmed.Touhami@gmail.com                  %
%                                                                  %
%------------------------------------------------------------------%
%------------------------------------------------------------------%
% REFERENCES:                                                      %
% M. Sadkane and A. Touhami. Modification of the Spectral Dichotomy%
%  Algorithm by the Circle. (2009)                                 %
%------------------------------------------------------------------%
%------------------------------------------------------------------%
% Example:
% A = gallery('frank',6);
% specdicho(A)
% [P]   = specdicho(A);
% [P,H] = specdicho(A);
% opts.mxiter = 20;
% opts.tol    = 1e-12;
% opts.r      = 5.38;
% [P,H] = specdicho(A,opts);
%
%------------------------------------------------------------------%
% PARSE INPUT VALUES
%------------------------------------------------------------------%
%
% 1. No input arguments, return help.
%
if (nargin < 1),
  error('Not enough input arguments.');
else
  if (nargin > 3),
    error('Too many input arguments.');
  end
end
%
%------------------------------------------------------------------%
% 2. Get matrix A. Check type (numeric or character)
%    and dimensions.
%
A = varargin{1};
if (~isnumeric(A)) | (ndims(A) ~= 2) | (isempty(A)),
  error('The first argument must be a nonempty matrix.');
 end 
[m,n] = size(A);
if m~= n
  error('Matrix A is not square.');
 
end
%------------------------------------------------------------------%
Bpresent= []; 
optsInd = [];
if (nargin >= 2),
  Bpresent = (nargin > 2) | ((nargin == 2) & (~isstruct(varargin{2})));
end

if Bpresent,
  B = varargin{2};
  if (~isnumeric(B)) | (ndims(B) ~= 2) | (~isequal(size(B),[n,n])),
    error('B must be a square matrix of the same size as A.');
  end
else
  B = eye(n);
end

if (nargin == 2 & ~Bpresent) | (nargin == 3 & Bpresent),
  optsInd = 3 - (~Bpresent);
end
%
%------------------------------------------------------------------%
% 3. Set all input options to default values.
%
geom   = 'C'         ;
c      = 0           ;
r      = 1           ;
a      = 5           ;
b      = 1           ;
p      = 1           ;
mxiter = 10          ;
tol    = 1e-10       ;
Ho     = []          ;
%
%
if ~isempty(optsInd),
  opts = varargin{optsInd};
  if ~isstruct(opts),
    error(sprintf(['The last argument must be a structure.']));
  end
  names = fieldnames(opts);
  J = strmatch('GEOM',upper(names),'exact');
  if ~isempty(J), geom = upper(getfield(opts,names{J})); end
  J = strmatch('C',upper(names),'exact');
  if ~isempty(J), c = getfield(opts,names{J}); end 
  J = strmatch('R',upper(names),'exact');
  if ~isempty(J), r = getfield(opts,names{J}); end  
  J = strmatch('A',upper(names),'exact');
  if ~isempty(J), a = getfield(opts,names{J}); end 
  J = strmatch('B',upper(names),'exact');
  if ~isempty(J), b = getfield(opts,names{J}); end
  J = strmatch('P',upper(names),'exact');
  if ~isempty(J), p = getfield(opts,names{J}); end 
  J = strmatch('MXITER',upper(names),'exact');
  if ~isempty(J), mxiter = getfield(opts,names{J}); end
  J = strmatch('TOL',upper(names),'exact');
  if ~isempty(J), tol = getfield(opts,names{J}); end
  J = strmatch('HO',upper(names),'exact');
  if ~isempty(J), Ho = getfield(opts,names{J}); end
end
%
% 4. Check for input errors in the data structure.
%
%-------------------------------------------------------------------%
%        Check for input errors in the geometry geom                %
%-------------------------------------------------------------------%
if    ~strcmp(geom,'C') & ~strcmp(geom,'I') & ...
      ~strcmp(geom,'E') & ~strcmp(geom,'P') ,
  error(sprintf(['The optional argument GEOM must be equal to\n' ...
                 '[C|c | E|e | I|i | P|p]: unit circle, ellipse\n'...
                 'imaginary axis or parabola.']));
end
%
%-------------------------------------------------------------------%
%     Check for input errors in the center c and radius r           %
%     when the input optional argument geom is equal to 'C'         %
%-------------------------------------------------------------------%
if (geom == 'C'),
  if ~isnumeric(c) | length(c) ~=1,
    error(['The center c  must be a' ...
	   '  real/complex number.']);
  end
  %
  if  ~isreal(r) | ~isnumeric(r) | length(r) ~=1 | r <= 0,
    error(sprintf([' The radius r must be \n' ...
           ' a nonzero numeric positive value.']));
  end
end
%
%-------------------------------------------------------------------%
%       Check for input errors in the  semi-axes a and b            %
%     when the input optional argument geom is equal to 'E'         %
%-------------------------------------------------------------------%
if (geom == 'E'),
   if  ~isreal(a) | ~isnumeric(a) | length(a) ~=1 | a <= 0 | ...
         ~isreal(b) | ~isnumeric(b) | length(b) ~=1 | b <= 0,
    error(sprintf([' The Semi-Major and Minor axes must be \n' ...
           ' a nonzero numeric positive values.']));
  else
    if (a < b),
      error(sprintf([' The Semi-Major axis a must be larger than \n' ...
             ' the Semi-Minor axis b.']));
    end
  end
end
%
%-------------------------------------------------------------------%
%          Check for input errors in the  parameter p               %
%     when the input optional argument geom is equal to 'P'         %
%-------------------------------------------------------------------% 
if (geom == 'P'),
  if  ~isreal(p) | ~isnumeric(p) | length(p) ~=1 | p <= 0,
    error(sprintf([' The parameter p must be \n' ...
	   ' a nonzero numeric positive value.']));
  end
end
%
%-------------------------------------------------------------------%
%         Check for input errors in the parameter mxiter            %
%-------------------------------------------------------------------%  
  if ~isreal(mxiter) | ~isnumeric(mxiter) | length(mxiter) ~=1 | mxiter <= 0,
    error(sprintf([' The optional argument mxiter  must be \n' ...
           ' a not null numeric positive value.']));
end
%
%-------------------------------------------------------------------%
%          Check for input errors in the parameter tol              %
%-------------------------------------------------------------------% 
if  ~isreal(tol) | ~isnumeric(tol) | length(tol) ~=1 | tol <= 0,
 error(sprintf([' The optional argument tol  must be \n' ...
	 ' a not null numeric positive value.']));
end

if (tol < eps),
  tol = eps; 
  disp(sprintf([' WARNING: The value of tol is set to ',num2str(tol)]));
end
%
%-------------------------------------------------------------------%
%          Check for input errors in the parameter Ho               %
%-------------------------------------------------------------------% 
Hsize = n;
if(geom == 'E') | (geom == 'P'),
  Hsize = 2*n; 
end
if ~isempty(Ho),
  if ~isnumeric(Ho),
    error('Incorrect type for input value(s) in the structure.');
  end
  % Check size first
  if any(size(Ho) ~= Hsize),  
    error([' Matrix Ho must be of the same size as A '...
           ' when GEOM = [C|c| I|i] and of twice the order '... 
           ' of A when GEOM = [E|e | P|p].']);
  end 
  %
  if (nnz(Ho) == 0),
    error('Matrix Ho contains all zeros.');
  end
  % Check for SPD
  if ~isequal(Ho,Ho'),
    error('Matrix Ho must be symmetric.');
  end
  [RHo, pHo] = chol(Ho);
  if (pHo ~= 0),
    error('Matrix Ho must be positive definite');
  end
else
  Ho = eye(Hsize);
end
%
%-------------------------------------------------------------------%
%     Check for non-singularity condition of the matrix B           %
%                  when GEOM = ['I' | 'P']                          %
%-------------------------------------------------------------------%
%
if(geom == 'I') | (geom == 'P'),
  if (rcond(B) < eps),
    error(sprintf(['Matrix B is singular,\n'...
                   'results may be inaccurate.']));
  end
end
%
%------------------------------------------------------------------%
%
if (geom == 'C'),
  disp(sprintf([' Circle of center c = (' ,num2str(real(c)) ',' ...
		,num2str(imag(c)) ') and radius r = ' ,num2str(r)])); 
  % Build A0, B0
  A0 = (1/r)*(A - c*B);
  B0 = B ;
elseif (geom == 'E'), 
  disp(sprintf([' Ellipse (X/a)^2 + (Y/b)^2 = 1\n' ...
                ' with a = ' ,num2str(a) ' and b = ' ,num2str(b)']));
  % Build A0, B0
  a1 = (a+b)/2;
  b1 = (a-b)/2; 
  A0 = [-b1*B , zeros(n,n) ; A , -b1*B];
  B0 = [a1*B  ,-A ; zeros(n,n) ,  a1*B];
  n  = size(A0,1);
elseif (geom == 'I'), 
  disp(sprintf([' Imaginary Axis\n']));
  % Build A0, B0
  A0 = A + B;
  B0 = A - B;
elseif (geom == 'P'),
  disp(sprintf([' Parabola Y^2 = 2*p*(p/2 - X) with p =' ...
                ' ' ,num2str(p)']));
  % Build A0, B0
  a1  = -sqrt(p/2);
  A0  = [(a1+1)*B  A ; eye(n) (a1+1)*eye(n)];
  B0  = [(a1-1)*B  A ; eye(n) (a1-1)*eye(n)];
  n   = size(A0,1);
end
%------------------------------------------------------------------%
%--------------------------%
%    Initialization        %
%--------------------------%
iter    = 1;
M1      = B0 - A0;
M2      = A0 + B0;
if (rcond(M1) < eps) |(rcond(M2) < eps),
  disp(['  At iteration ',num2str(iter)]);
  error(['  the used contour does not realize a spectral dichotomy.']);
 end
[L1,U1] = lu(M1);
[L2,U2] = lu(M2);
%----------------------------%
%    Solution for X, Y       %
%----------------------------%
X1      =  A0/U1;
X       =  X1/L1;
Y1      =  B0/U1 ;
Y       =  Y1/L1;
%----------------------------%
% Solution for Delta1,Nabla1 %
%----------------------------%
X2      = L2\X ;
Delta   = U2\X2;
Y2      = L2\Y ;
Nabla   = U2\Y2;
%----------------------------%
%    Computation of H1       %
%----------------------------%
Hj      = Delta'*Ho*Delta+Nabla'*Ho*Nabla;
Z1      = Delta;
Z2      = Nabla;
%---------------------------%
%        Iterations         %
%---------------------------%
while (iter<mxiter)
  iter   = iter +1;
  Hjold  = Hj;
  %---------------------------%
  %     Update  Aj            %
  %---------------------------%
  Aj    = -A0*Z1;
  %----------------------------%
  %   Solution for Deltaj      %
  %----------------------------%
  M = 2*Aj-eye(n);
  if (rcond(M)<eps),
    disp(['  At iteration ',num2str(iter)]);
    error(['  the used contour does not realize a spectral dichotomy.']);
  end
  Delta  = M\Aj;
  Nabla  = eye(n)-Delta;
  %-----------------------------%
  %  Computation of Hj, Z1, Z2j %
  %-----------------------------%
  Hj    = Delta'*Hj*Delta+Nabla'*Hj*Nabla;
  Z1    = Z1*Delta;
  Z2    = Z2*Nabla;
  %-----------------------------%
  %     Stopping criterion      %
  %-----------------------------%
  if ((norm(Hj - Hjold)/norm(Hj)) <= tol),
    disp([' At iteration ',num2str(iter)]);
    disp([' convergence to the desired tolerance tol = ',num2str(tol)]);
    break
  end
end
H = Hj;
P = Z2*B0;
%
% The spectral projector when GEOM = ['E' | 'P']
%
if (geom == 'E'),
  P = P(1:Hsize/2,1:Hsize/2) +  P((Hsize/2)+1:Hsize,(Hsize/2)+1:Hsize) - eye(Hsize/2,Hsize/2);
else
  if (geom == 'P'),
    P =  2*P(1:Hsize/2,1:Hsize/2) - eye(Hsize/2,Hsize/2);
  end
end
%------------------------------------------------------------------%
