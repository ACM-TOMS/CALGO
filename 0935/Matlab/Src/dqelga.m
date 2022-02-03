function [n,epstab,result,abserr,res3la,nres]=dqelga(n,epstab,res3la,nres)
%***BEGIN PROLOGUE  DQELG
%***SUBSIDIARY
%***PURPOSE  The routine determines the limit of a given sequence of
%            approximations, by means of the Epsilon algorithm of
%            P.Wynn. An estimate of the absolute error is also given.
%            The condensed Epsilon table is computed. Only those
%            elements needed for the computation of the next diagonal
%            are preserved.
%***LIBRARY   SLATEC
%***TYPE      doubleprecision (QELG-S, DQELG-D)
%***KEYWORDS  CONVERGENCE ACCELERATION, EPSILON ALGORITHM, EXTRAPOLATION
%***AUTHOR  Piessens, Robert
%             Applied Mathematics and Programming Division
%             K. U. Leuven
%           de Doncker, Elise
%             Applied Mathematics and Programming Division
%             K. U. Leuven
%***DESCRIPTION
%
%           Epsilon algorithm
%           Standard fortran subroutine
%           doubleprecision version
%
%           PARAMETERS
%              N      - Integer
%                       EPSTAB(N) contains the new element in the
%                       first column of the epsilon table.
%
%              EPSTAB - doubleprecision
%                       Vector of dimension 52 containing the elements
%                       of the two lower diagonals of the triangular
%                       epsilon table. The elements are numbered
%                       starting at the right-hand corner of the
%                       triangle.
%
%              RESULT - doubleprecision
%                       Resulting approximation to the integral
%
%              ABSERR - doubleprecision
%                       Estimate of the absolute error computed from
%                       RESULT and the 3 previous results
%
%              RES3LA - doubleprecision
%                       Vector of dimension 3 containing the last 3
%                       results
%
%              NRES   - Integer
%                       Number of calls to the routine
%                       (should be zero at first call)
%
%***REVISION HISTORY  (YYMMDD)
%   800101  DATE WRITTEN
%   890531  Changed all specific intrinsics to generic.  (WRB)
%   890531  REVISION DATE from Version 3.2
%   891214  Prologue converted to Version 4.0 format.  (BAB)
%   900328  Added TYPE section.  (WRB)
%   10/02/2010  Cleaned up for use in IIPBF. Simplified input/output 
%               variables (Jung Hun Kim)
%***end PROLOGUE  DQELG
%

% Reference:
% Ratnanather, J. T., Kim, J. H., Zhang, S., Davis, A. M. J., and Lucas, 
% S. K. 2013. Algorithm XXX: IIPBF, a MATLAB toolbox for infinite integral 
% of products of two Bessel functions. To appear in ACM Transactions on 
% Mathematical Software.

persistent delta1 delta2 delta3 e0 e1 e1abs e2 e3 epmach epsinf err1 err2 err3 error i ib ib2 ie indx k1 k2 k3 limexp newelm num oflow res ss tol1 tol2 tol3 ;

if isempty(delta1), delta1=0; end;
if isempty(delta2), delta2=0; end;
if isempty(delta3), delta3=0; end;
if isempty(epmach), epmach=0; end;
if isempty(epsinf), epsinf=0; end;
if isempty(error), error=0; end;
if isempty(err1), err1=0; end;
if isempty(err2), err2=0; end;
if isempty(err3), err3=0; end;
if isempty(e0), e0=0; end;
if isempty(e1), e1=0; end;
if isempty(e1abs), e1abs=0; end;
if isempty(e2), e2=0; end;
if isempty(e3), e3=0; end;
if isempty(oflow), oflow=0; end;
if isempty(res), res=0; end;
if isempty(ss), ss=0; end;
if isempty(tol1), tol1=0; end;
if isempty(tol2), tol2=0; end;
if isempty(tol3), tol3=0; end;
if isempty(i), i=0; end;
if isempty(ib), ib=0; end;
if isempty(ib2), ib2=0; end;
if isempty(ie), ie=0; end;
if isempty(indx), indx=0; end;
if isempty(k1), k1=0; end;
if isempty(k2), k2=0; end;
if isempty(k3), k3=0; end;
if isempty(limexp), limexp=0; end;
if isempty(newelm), newelm=0; end;
if isempty(num), num=0; end;
%
%           LIST OF MAJOR VARIABLES
%           -----------------------
%
%           E0     - THE 4 ELEMENTS ON WHICH THE COMPUTATION OF A NEW
%           E1       ELEMENT IN THE EPSILON TABLE IS BASED
%           E2
%           E3                 E0
%                        E3    E1    NEW
%                              E2
%           NEWELM - NUMBER OF ELEMENTS TO BE COMPUTED IN THE NEW
%                    DIAGONAL
%           ERROR  - ERROR = ABS(E1-E0)+ABS(E2-E1)+ABS(NEW-E2)
%           RESULT - THE ELEMENT IN THE NEW DIAGONAL WITH LEAST VALUE
%                    OF ERROR
%
%           MACHINE DEPENDENT CONSTANTS
%           ---------------------------
%
%           EPMACH IS THE LARGEST RELATIVE SPACING.
%           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
%           LIMEXP IS THE MAXIMUM NUMBER OF ELEMENTS THE EPSILON
%           TABLE CAN CONTAIN. IF THIS NUMBER IS REACHED, THE UPPER
%           DIAGONAL OF THE EPSILON TABLE IS DELETED.
%
%***FIRST EXECUTABLE STATEMENT  DQELG
% epmach is the largest relative spacing.
epmach = 1.0D-15;
% uflow is the smallest positive magnitude.
uflow = 1.0D-100;
% oflow is the largest magnitude.
oflow = 1.0D+100;

nres = fix(nres + 1);
abserr = oflow;
result = epstab(n);
if( n>=3 )
    limexp = 50;
    epstab(n+2) = epstab(n);
    newelm =fix(fix((n-1)./2));
    epstab(n) = oflow;
    num = fix(n);
    k1 = fix(n);
    for i = 1 : newelm;
        k2 = fix(k1 - 1);
        k3 = fix(k1 - 2);
        res = epstab(k1+2);
        e0 = epstab(k3);
        e1 = epstab(k2);
        e2 = res;
        e1abs = abs(e1);
        delta2 = e2 - e1;
        err2 = abs(delta2);
        tol2 = max(abs(e2),e1abs).*epmach;
        delta3 = e1 - e0;
        err3 = abs(delta3);
        tol3 = max(e1abs,abs(e0)).*epmach;
        if( err2>tol2 || err3>tol3 )
            e3 = epstab(k1);
            epstab(k1) = e1;
            delta1 = e1 - e3;
            err1 = abs(delta1);
            tol1 = max(e1abs,abs(e3)).*epmach;
            %
            %           IF TWO ELEMENTS ARE VERY CLOSE TO EACH OTHER, OMIT
            %           A PART OF THE TABLE BY ADJUSTING THE VALUE OF N
            %
            if( err1>tol1 && err2>tol2 && err3>tol3 )
                ss = 0.1d+01./delta1 + 0.1d+01./delta2 - 0.1d+01./delta3;
                epsinf = abs(ss.*e1);
                %
                %           TEST TO DETECT IRREGULAR BEHAVIOUR IN THE TABLE, AND
                %           EVENTUALLY OMIT A PART OF THE TABLE ADJUSTING THE VALUE
                %           OF N.
                %
                if( epsinf>0.1d-03 )
                    %
                    %           COMPUTE A NEW ELEMENT AND EVENTUALLY ADJUST
                    %           THE VALUE OF RESULT.
                    %
                    res = e1 + 0.1d+01./ss;
                    epstab(k1) = res;
                    k1 = fix(k1 - 2);
                    error = err2 + abs(res-e2) + err3;
                    if( error<=abserr )
                        abserr = error;
                        result = res;
                    end;
                    continue;
                end;
            end;
            n = fix(i + i - 1);
            % ***JUMP OUT OF DO-LOOP
            break;
        else
            %
            %           IF E0, E1 AND E2 ARE EQUAL TO WITHIN MACHINE
            %           ACCURACY, CONVERGENCE IS ASSUMED.
            %           RESULT = E2
            %           ABSERR = ABS(E1-E0)+ABS(E2-E1)
            %
            result = res;
            abserr = err2 + err3;
            % ***JUMP OUT OF DO-LOOP
            abserr = max(abserr,0.5d+01.*epmach.*abs(result));
            return;
        end;
    end;
    %
    %           SHIFT THE TABLE.
    %
    if( n==limexp )
        n = fix((2.*(fix(limexp./2))) - 1);
    end;
    ib = 1;
    if(((fix(num./2)).*2)==num )
        ib = 2;
    end;
    ie = fix(newelm + 1);
    for i = 1 : ie;
        ib2 = fix(ib + 2);
        epstab(ib) = epstab(ib2);
        ib = fix(ib2);
    end; i = fix(ie+1);
    if( num~=n )
        indx = fix(num - n + 1);
        for i = 1 : n;
            epstab(i) = epstab(indx);
            indx = fix(indx + 1);
        end; i = fix(n+1);
    end;
    if( nres>=4 )
        %
        %           COMPUTE ERROR ESTIMATE
        %
        abserr = abs(result-res3la(3)) + abs(result-res3la(2))+ abs(result-res3la(1));
        res3la(1) = res3la(2);
        res3la(2) = res3la(3);
        res3la(3) = result;
    else
        res3la(nres) = result;
        abserr = oflow;
    end;
end;
abserr = max(abserr,0.5d+01.*epmach.*abs(result));
end %subroutine dqelg
%DECK DQFORM
