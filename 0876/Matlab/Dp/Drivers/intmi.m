function [x,w] = intmi(a,b,n)
% Calculate the nodes and weights for the integration
%
%           B                       N
%       INT   F(X)*DX  APPROX=  SUM   W(I)*F(X(I))
%           A                      I=1
%
% The method is due to Iri, Moriguti, and Takasawa.  It uses a change of
% variable to render harmless most endpoint singularities, followed by use
% of the trapezoidal rule. This program is a translation of the original 
% Fortran program of Ian Robinson and Elise de Doncker, cited below.
%
% This function is called by both TestLogSing and TestAlgSing.  It would
% be convenient to place it in a private directory.

    x = zeros(n,1); w = zeros(n,1);

    % Initialize
    h = 2/(n+1);
    if n == 2*floor(n/2);
        m = n/2;
    else
        m = (n+1)/2;
    end
    fac = b-a;
    fac2 = h*(b-a);
    
    % Calculate nodes and weights
    for i=1:m
        [si,ti] = trnsfm(i*h-1.0);
        x(i) = a + ti*fac;
        x(n+1-i) = b - ti*fac;
        w(i) = fac2*si;
        w(n+1-i) = w(i);
    end
end %intmi

function [si,fi] = trnsfm(z)

% For given z, calculate
%
%        si(z) = con*exp(-4/(1-z*z))
%
%     and
%                        z
%        fi(z) = Integral  si(t)*dt
%                       -1
%
% where con is chosen to normalize fi(1)=1.
% fi(z) is computed using a truncated Chebyshev series.  
%
% This routine is taken from the following paper: 
% Ian Robinson and Elise de Doncker, "Automatic computation of improper
% integrals over a bounded or unbounded planar region, Computing (1981).

% Coefficients for some Chebyshev expansions.
q1008= [0.1307592530313668E-01   
0.8580903946529252E-02
0.1969565496511163E-02 
-0.6785908047534174E-04
0.5157989084218512E-05 
-0.3254283578286669E-06
0.3009165432294816E-07 
-0.2853960792992456E-08
0.3125316663674964E-09 
-0.369071915636499E-10
0.47226104267733E-11   
-0.6452090127750E-12
0.935240772542E-13     
-0.142880027692E-13
0.22888743282E-14      
-0.3827976860E-15
0.665889416E-16        
-0.120097537E-16
0.22395647E-17]';

q0806 = [0.7530091699129213E-01  
0.2215896190278894E-01
0.1517662771681320E-02 
-0.1374204580631322E-04
0.2501181491115358E-05 
-0.2491206397236787E-07
0.5430948732810670E-08 
-0.1406070367942514E-09
0.160397857131033E-10  
-0.7706914621139E-12
0.656131517151E-13     
-0.43257167630E-14
0.3459964512E-15       
-0.264283638E-16
0.21607480E-17]';

q0603 = [0.2285305746758029E00   
0.5670149715650661E-01
0.3881611773394649E-02  
0.1483758828946990E-03
0.1986416462810431E-04  
0.1166284710859293E-05
0.1048168134503124E-06  
0.6572636994171403E-08
0.5344588684897204E-09  
0.335115172128537E-10
0.26225503158527E-11    
0.1606252080762E-12
0.124685606769E-13      
0.73251000E-15
0.579829277E-16         
0.31817034E-17]';

q0300 = [0.5404466499651320E+00  
0.1034761954005156E+00
0.9090551216566596E-02  
0.9130257654553327E-03
0.1026685176623610E-03  
0.1035244509501565E-04
0.1034925154934048E-05  
0.1002050044392230E-06
0.9552837592432324E-08  
0.8941334072231475E-09
0.825727197051832E-10   
0.75255600668576E-11
0.6783483380205E-12     
0.605164239084E-13
0.53497536298E-14       
0.4689365198E-15
0.407909118E-16         
0.35230149E-17]';

q = [q1008, q0806, q0603, q0300];
% q(1:19) = q1008; q(20:34) = q0806; q(35:50) = q0603; q(51:68) = q0300;

% q contains the coefficients of the Chebyshev expansion of fi.  To
% achieve greater than 16-digit accuracy, it is necessary to increase the
% number of digits in these coefficients and also the number of
% coefficients used in the Chebyshev series.  This requires changes to the
% size of q and the coefficient arrays from which it is defined.  Refer to
% appendix 3 of the following for details:
%       A. Haegemans, Algorithm 34: An algorith for the automatics
%       integration over a triangle, Computing 19 (1977), pp. 179-187.

arg = [1,20,35,51,69];
uflow = realmin; eoflow = -log(realmin);
divsor = 7.0298584066096E-3; dv1 = 2*divsor; dv2 = 4*divsor;
a1 = [20,20,40/3,40/3]; a2 = [18,14,6,2];

    x = -abs(z);
    % Compute si, eliminating first the easier cases.    
    if x <= -1
        si = 0;
        if z > 0
            fi = 1;
        else
            fi = 0;
        end
        return
    else
        expont = 4/(1-x*x);
        if expont < eoflow
            si = exp(-expont);
        else
            si = 0;
            if z > 0
                fi = 1;
            else
                fi = 0;
            end
            return            
        end
    end

    % Determine the coefficients to be used in the Chebyshev series.
    if x >= -.3
        ctr = 4;
    elseif x >= -.6
        ctr = 3;
    elseif x >= -.8
        ctr = 2;
    else
        ctr = 1;
    end

    y = a1(ctr)*x + a2(ctr);
    kfirst = arg(ctr); klast = arg(ctr+1)-2;
    
    % Compute fi.
    t2 = 0; t3 = q(klast+1);
    ksum = kfirst + klast;
    for k=kfirst:klast
        krev = ksum - k;
        t1 = t2; t2 = t3;
        t3 = y*t2 - t1 + q(krev);
    end
    
    fac = (t3 - t1)/dv2;
    if fac > realmin/si
        fi = fac*si;
    else
        fi = 0;
    end
    si = si/dv1;
    if z > 0
        fi = 1 - fi;
    end
end % trnsfm
