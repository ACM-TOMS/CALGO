function tGJ=gettGJ
% Computes relative efficiency of igamma compared to besselj

% ------------------------------------------------------------
% The parameter tGJ in the cost function, defined as the
% relative efficiency of the incomplete Gamma function
% compared to the Bessel function, reflects the fact that
% the cost is dominated by these two functions.
% 
% While the incomplete Gamma functions is always called with
% a scalar argument, the Bessel function is called with a
% vector argument.  And that has a serious impact in Matlab.
% This is taken into account in the second loop below.
% 
% There are other factors that cannot easily be taken into
% account.  E.g., The number of incomplete Gamma functions
% that is really evaluated is often lower than in eq. (10)
% because we cut the series.
% 
% Our advice is to set tGJ in besselint.m lower than the
% result of this script. E.g. this script suggested tGJ=19
% and we finally selected tGJ=15, considering the total
% execution time of experiments.m
% ------------------------------------------------------------
%   Authors: Joris Van Deun & Ronald Cools,
%            Dept. of Computer Science, K.U.Leuven, Belgium.
%
%   Software revision date: September 2, 2005
% ------------------------------------------------------------


% warm up 
igamma(-0.5,-i);
besselj(0,1);

total=17*120;
stepsize=17*8; % = 17 * unrolled levels

n=22;
for k=1:n,
    q=0;
    for j=1:17*8:total,
        % loop unrolling used to avoid q=0 on fast machines
        % will cause trouble in the future, if speed increases!
        t=cputime;
        besselj(0,j+[0:16]);
        besselj(0,j+[0:16]);
        besselj(0,j+[0:16]);
        besselj(0,j+[0:16]);
        besselj(0,j+[0:16]);
        besselj(0,j+[0:16]);
        besselj(0,j+[0:16]);
        besselj(0,j+[0:16]);
        q=q+cputime-t;
    end
    
    p=0;
    for j=1:total,
        t=cputime;igamma(-0.5,-i*j);p=p+cputime-t;
    end
    
    r(k)=p/q;
end
% take average, ignoring 2 extremes.
r=sort(r);
tGJ=sum(r(2:n-1))/(n-2);
disp('Suggested value for tGJ = ');

