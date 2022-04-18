%EXPERIMENTS to illustrate the use of BESSELINT

format long e

% Part I: testing the results
disp('Part I: explicitly known examples');

% Example 1. One factor.
% Example 1.1
disp(sprintf('\nEx. 1.1: a=[1]; nu=[1/3]; m=0; \nExact answer: 1'));
I=besselint(1,1/3,0)

% Example 1.2
q=1/3;nu=-1/4;I12=2^q*gamma((nu+q+1)/2)/gamma((nu-q+1)/2);
disp(sprintf('\nEx. 1.2: a=[1]; nu=[-1/4]; m=1/3; \nExact answer: %1.15e',I12));
I=besselint(1,-1/4,1/3)

% Example 1.3 Abel summability
disp(sprintf('\nEx. 1.3: a=[1]; nu=0; m=2; \nExact answer: -1'));
I=besselint(1,0,2)

% Example 1.4 Abel summability
disp(sprintf('\nEx. 1.4: a=[1]; nu=0; m=4; \nExact answer: 9'));
I=besselint(1,0,4)

% Example 2. Two factors.
% Example 2.1
disp(sprintf('\nEx. 2.1: a=[1 5]; nu=[0 1]; m=0; \nExact answer: 1/5'));
I=besselint([1 5],[0 1],0)

% Example 2.2
disp(sprintf('\nEx. 2.2: a=[5 1]; nu=[0 1]; m=0; \nExact answer: 0'));
I=besselint([5 1],[0 1],0)

% Example 2.3
a=3;mu=-1/2;nu=1/3;I23=((a^2-1)/2)^(nu-mu-1)/(a^nu*gamma(nu-mu));
disp(sprintf('\nEx. 2.3: a=[1 3]; nu=[-1/2 1/3]; m=1/6; \nExact answer: %1.15e',I23));
I=besselint([1 3],[-1/2 1/3],1/6)

% Example 2.4
disp(sprintf('\nEx. 2.4: a=[3 1]; nu=[-1/2 1/3]; m=1/6; \nExact answer: 0'));
I=besselint([3 1],[-1/2 1/3],1/6)

% Example 3. Three factors.
% Example 3.1
I31=1/(pi*sqrt(6));
disp(sprintf('\nEx. 3.1: a=sqrt([2 3 5]); nu=0; m=1; \nExact answer: %1.15e',I31));
I=besselint(sqrt([2 3 5]),0,1)

% Example 3.2
disp(sprintf('\nEx. 3.2: a=sqrt([2 3 11]); nu=0; m=1; \nExact answer: 0'));
I=besselint(sqrt([2 3 11]),0,1)

% Example 3.3
I33=1/(pi*sqrt(5));
disp(sprintf('\nEx. 3.3: a=sqrt([2 3 5]); nu=1; m=0; \nExact answer: %1.15e',I33));
I=besselint(sqrt([2 3 5]),1,0)

% Example 3.4
alfa1=pi-acos(3/sqrt(15));alfa2=pi-acos(2/sqrt(10));alfa3=pi/2;
delta=sqrt(6);
I34=(cos(alfa2-2*alfa1)+cos(2*alfa3+3*alfa2)+cos(3*alfa1+alfa3))/(3*pi*delta);
disp(sprintf('\nEx. 3.4: a=sqrt([2 3 5]); nu=[1 2 -3]; m=1; \nExact answer: %1.15e',I34));
I=besselint(sqrt([2 3 5]),[1 2 -3],1)

% Example 3.5
disp(sprintf('\nEx. 3.5: a=sqrt([2 3 11]); nu=[1 2 -3]; m=1; \nExact answer: 0'));
I=besselint(sqrt([2 3 11]),[1 2 -3],1)

% Example 4. Four factors.
% Example 4.1
I41=ellipke(1/16*(8+23*sqrt(5/42)))/(pi^2*sqrt(sqrt(210)));
disp(sprintf('\nEx. 4.1: a=sqrt([2 3 5 7]); nu=0; m=1; \nExact answer: %1.15e',I41));
I=besselint(sqrt([2 3 5 7]),0,1)

% Example 4.2
nu=5/4;
I42=gamma(2*nu)*gamma(nu)/(2*pi*gamma(3*nu)*gamma(nu+1/2)^2);
disp(sprintf('\nEx. 4.2: a=[1 1 1 1]; nu=5/4; m=-3/2; \nExact answer: %1.15e',I42));
I=besselint([1 1 1 1],5/4,-3/2)

% Example 5. Five factors.
% Example 5.1 Problem 8  from [2], 'exact' result known (thx to Dirk Laurie)
I51=0.06106434990872167061868472;
disp(sprintf('\nEx. 5.1: a=sqrt([2 3 5 7 11]); nu=0;  m=1; \nExact answer: %1.15e',I51));
I=besselint(sqrt([2 3 5 7 11]),0,1)

% Example 5.2 'exact' result known (thx to Dirk Laurie)
I52=1.7024879933914e-2;
disp(sprintf('\nEx. 5.2: a=sqrt([2 3 5 7 11]); nu=0;  m=2; \nExact answer: %1.13e',I52));
I=besselint(sqrt([2 3 5 7 11]),0,2)

% Example 6. k factors.
% Example 6.1
disp(sprintf('\nEx. 6.1: a=[8 2.5 2 1.5 1]; nu=1; m=-2; \nExact answer: 0'));
besselint([8 2.5 2 1.5 1],1,-2)

% Example 6.2
b=8;mu=1/3;nu=-1/4;a=[2.5 2 1.5 1 0.5];
I62=2^(mu-1)*gamma(mu)/b^mu*prod((a/2).^nu/(gamma(nu+1)));
disp(sprintf('\nEx. 6.2: a=[8 2.5 2 1.5 1 0.5]; nu=[1/3 -1/4 -1/4 -1/4 -1/4 -1/4]; m=7/12; \nExact answer: %1.15e',I62));
I=besselint([8 2.5 2 1.5 1 0.5],[1/3 -1/4 -1/4 -1/4 -1/4 -1/4],1/3+5/4-1)


% Part II: comparing requested, estimated and actual accuracy
disp(sprintf('\n----------------------------------------------\n'));
disp('Part II: how reliable are the error estimates?');

% Example 7. Requested relative precision
% Example 7.1
disp(sprintf('\nEx. 7.1: a=sqrt([2 3 5]); nu=0; m=1'));
disp('Requested relative precision: 1e-14');
[I,e]=besselint(sqrt([2 3 5]),0,1,1e-14);
disp(sprintf('Estimated relative precision: %1.2e',e(1)));
disp(sprintf('Actual relative precision: %1.2e',abs(I/I31-1)));

% Example 7.2
disp(sprintf('\nEx. 7.2: a=sqrt([2 3 5]); nu=[1 2 -3]; m=1'));
disp('Requested relative precision: 1e-6');
[I,e]=besselint(sqrt([2 3 5]),[1 2 -3],1,1e-6);
disp(sprintf('Estimated relative precision: %1.2e',e(1)));
disp(sprintf('Actual relative precision: %1.2e',abs(I/I34-1)));

% Example 7.3
disp(sprintf('\nEx. 7.3: a=[1]; nu=0; m=2'));
disp('Requested relative precision: 1e-10');
[I,e]=besselint(1,0,2,1e-10);
disp(sprintf('Estimated relative precision: %1.2e',e(1)));
disp(sprintf('Actual relative precision: %1.2e',abs(I+1)));

% Example 7.4
disp(sprintf('\nEx. 7.4: a=[1]; nu=0; m=4'));
disp('Requested relative precision: 1e-2');
[I,e]=besselint(1,0,4,1e-2);
disp(sprintf('Estimated relative precision: %1.2e',e(1)));
disp(sprintf('Actual relative precision: %1.2e',abs(I/9-1)));

% Example 8. Requested absolute precision
% Example 8.1 
disp(sprintf('\nEx. 8.1: a=[1]; nu=[1/3]; m=0'));
disp('Requested absolute precision: 1e-14');
[I,e]=besselint(1,1/3,0,[],1e-14);
disp(sprintf('Estimated absolute precision: %1.2e',e(1)));
disp(sprintf('Actual absolute precision: %1.2e',abs(I-1)));

% Example 8.2
disp(sprintf('\nEx. 8.2: a=sqrt([2 3 5 7]); nu=0; m=1'));
disp('Requested absolute precision: 1e-6');
[I,e]=besselint(sqrt([2 3 5 7]),0,1,[],1e-6);
disp(sprintf('Estimated absolute precision: %1.2e',e(1)));
disp(sprintf('Actual absolute precision: %1.2e',abs(I-I41)));

% Example 8.3
disp(sprintf('\nEx. 8.3: a=[1]; nu=[-1/4]; m=1/3'));
disp('Requested absolute precision: 1e-2');
[I,e]=besselint(1,-1/4,1/3,[],1e-2);
disp(sprintf('Estimated absolute precision: %1.2e',e(1)));
disp(sprintf('Actual absolute precision: %1.2e',abs(I-I12)));

% Example 8.4
disp(sprintf('\nEx. 8.4: a=sqrt([2 3 5 7 11]); nu=0; m=2'));
disp('Requested absolute precision: 1e-10');
[I,e]=besselint(sqrt([2 3 5 7 11]),0,2,[],1e-10);
disp(sprintf('Estimated absolute precision: %1.2e',e(1)));
disp(sprintf('Actual absolute precision: %1.2e',abs(I-I52)));

% Part III: testing some special cases
disp(sprintf('\n----------------------------------------------\n'));
disp('Part III: some unusual cases');

% Example 9. Discontinuous cases alpha = 0
% Example 9.1
disp(sprintf('\nEx. 9.1: a=[1 1]; nu=[1 4];  m=1; \nExact answer: -Inf'));
I=besselint([1 1],[1 4],1)

% Example 9.2
disp(sprintf('\nEx. 9.2: a=[1 1]; nu=[0 1];  m=0; \nExact answer: 1/2'));
I=besselint([1 1],[0 1],0)

% Example 9.3
disp(sprintf('\nEx. 9.3: a=[1 2 3]; nu=0; m=1; \nExact answer: Inf'));
I=besselint([1 2 3],0,1)

% Example 9.4
disp(sprintf('\nEx. 9.4: a=[1 2 3]; nu=[0 1 1/2]; m=1; \nExact answer: finite but unknown'));
I=besselint([1 2 3],[0 1 1/2],1)

% Example 9.5
disp(sprintf('\nEx. 9.5: a=[1 1]; nu=[-11/2 11/2]; \nExact answer: finite if m<3'));
m=2
I=besselint([1 1],[-11/2 11/2],m)
m=3-10*eps
I=besselint([1 1],[-11/2 11/2],m)
m=3
I=besselint([1 1],[-11/2 11/2],m)

% Example 9.6 Artificially constructed example
disp(sprintf('\nEx. 9.6: a=[1 2 3]; nu=[0 0.8027756377319946 0.3027756377319946]; \nExact answer: finite if m<3.5'));
m=2
I=besselint([1 2 3],[0 0.8027756377319946 0.3027756377319946],m)
m=3.5-10*eps
I=besselint([1 2 3],[0 0.8027756377319946 0.3027756377319946],m)
m=3.5
I=besselint([1 2 3],[0 0.8027756377319946 0.3027756377319946],m)

% Example 10. Continuous cases alpha = 0
% Example 10.1 'exact' result known (thx to Dirk Laurie)
I101 = 0.475270173593537344;
disp(sprintf('\nEx. 10.1: a=[1 2 3]; nu=0; m=0; \nExact answer: %1.15e',I101));
I=besselint([1 2 3],0,0)

% Example 10.2 'exact' result known (thx to Dirk Laurie)
I102 = 0.44371090379604393331635692;
disp(sprintf('\nEx. 10.2: a=[sqrt(2) sqrt(3) sqrt(2)+sqrt(3)]; nu=0; m=0; \nExact answer: %1.15e',I102));
I=besselint([sqrt(2),sqrt(3),sqrt(2)+sqrt(3)],0,0)

% Example 11. Large parameters. Computations may take much time.
% Example 11.1 Coefficients
disp(sprintf('\nEx. 11.1: a=[1 10 1000]; nu=0; m=1; \nExact answer: 0'));
[I,e]=besselint([1 10 1000],0,1)

% Example 11.2 Factors in integrand
disp(sprintf('\nEx. 11.2: a=sqrt([2 3 5 7 11 13 17 19 23]); nu=0; m=4; \nExact answer: unknown'));
[I,e]=besselint(sqrt([2 3 5 7 11 13 17 19 23]),0,4)

% Example 11.3 Order
disp(sprintf('\nEx. 11.3: a=1; nu=100; m=0; \nExact answer: 1'));
[I,e]=besselint(1,100,0)
