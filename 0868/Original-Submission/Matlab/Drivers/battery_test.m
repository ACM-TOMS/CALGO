%battery_test
%This script offers to test the following four codes
%da0glob, da2glob, da3glob and quadl on
%the 23 problems in the battery test for twelve different
%accuracies: If you want more detailed information about each problem you
%may activate the last command (appears now as a comment) 
%before the end of the outer loop: Problem_number,work,extra_digits
%
%This code takes less than 50 seconds on an
%Dell Optiplex 260 (pentium 4) running Matlab version 7.0.
%
%Modified by T. O. Espelid 24/09/06

%Some of these problems imply division by zero for certain arguments
%and will therefore give a warning and this will happen for all codes and
%all accuracies. In order to reduce the amount of output we turn the 
%warning off
  format short
  warning('off')
for Problem_number=1:23,
% Get the correct answer and the interval of this problem:
   [exact,a,b]=correct(Problem_number);
   tol=1;
   for i=1:12
     tol=tol*10^-1;
%Here we use anonymous function calls and give the extra parameter
%Problem_number to the function problems(x,i).
    [q(i,1),work(i,1)]=da0glob(@(x)problems(x,Problem_number),[a,b],tol);
    [q(i,2),work(i,2)]=da2glob(@(x)problems(x,Problem_number),[a,b],tol);
    [q(i,3),work(i,3)]=da3glob(@(x)problems(x,Problem_number),[a,b],tol);
    [q(i,4),work(i,4)]=quadl(@(x)problems(x,Problem_number),a,b,tol);
   end
   true_error=abs(q-exact);
% the extra accuracy achieved: counted in extra digits
% A negative result implies a failure
   extra_digits= -log10(true_error+eps)-(1:12)'*ones(1,4);
   failures(Problem_number,:) = [Problem_number,sum((extra_digits<0))];
% Exclude those results that fail to achieve the given tolerance
   A=work;
   A(extra_digits<0)=100000;
   A=A';
% Then find the minimum number of function values
% and all codes achieving this minimum
   minvalues=(min(A))';minvalues(minvalues==100000)=0;
   B = (work==(minvalues*ones(1,4)));
% The number of times a code has achieved this minimum
   best_code(Problem_number,:)=[Problem_number,sum(B)];
% Activate the next line to get more information.
%  Problem_number,work,extra_digits
end

failures
sum_failures=sum(failures(:,2:5))
best_code
sum_best_code=sum(best_code(:,2:5))
%In order to illustrate the information available
%if you activate the last comment in the outer loop
%(giving us more information about problem 23)
Problem_number,work,extra_digits
warning('on')

%-------------------------------------------------------------------------
%%This codes gives the following output on a Dell Optiplex GX260
%%Running Matlab version 7.0
%%
%%We list, for each of the 23 different problems, (1) how many failures
%%out of the 12 accuracies each code has (2) how many times each code
%%uses the fewest number of function evaluations for these twelve accuracies
%%Finally we sum these numbers over all 23 problems for each of the four
%%codes.
%%
%%problem number da0glob da2glob da3glob quadl

%>>battery_test
%
%failures =

%     1     0     0     0     0
%     2     0     0     0     5
%     3     0     0     0     2
%     4     0     0     0     0
%     5     0     0     0     0
%     6     0     0     0     0
%     7     0     0     0     5
%     8     0     0     0     0
%     9     0     0     0     0
%    10     0     0     0     0
%    11     0     0     0     0
%    12     0     0     0     0
%    13     0     0     0     3
%    14     0     0     0     0
%    15     0     0     0     0
%    16     0     0     0     0
%    17     2     2     2     3
%    18     0     0     0     0
%    19     0     0     0     2
%    20     0     0     0     0
%    21     5     9     9     6
%    22     0     0     0     0
%    23     0     0     0     3


%sum_failures =

%     7    11    11    29


%best_code =

%     1    12    12    12     0
%     2    11    11    11     1
%     3    10    12    12     0
%     4    12    12    12     0
%     5     5     9     8     3
%     6    12     3     3     0
%     7     7    12     9     0
%     8     7    11    11     1
%     9     3    11     7     1
%    10     8    12    11     0
%    11    12    12    12     0
%    12    12    12    12     0
%    13     4    11     8     1
%    14     6    12    11     0
%    15     5    12    11     0
%    16     7    11     8     0
%    17     0     6     2     4
%    18     3    11     9     1
%    19     9    12    10     0
%    20     6     9     9     3
%    21     5     1     1     2
%    22     2    10     8     2
%    23     6    11     9     1


%sum_best_code =
%
%   164   235   206    20
%
%-------------------------------------------------------------------------
%%This is the first part of the output
%%COMMENTS to the first part:
%%Counting failures: quadl has a total of 29 failures out
%%of the 23*12=276 test cases, while da0glob has 7 failures.
%%Among these four codes da2glob is the best code totally ( counting the 
%%number of  function evaluations): in 235 out of 276 cases.
%%
%%------------------------------------------------------------------------
%%Here comes the second part of the output: for 12 different accuracies
%%and the four codes: da0glob da2glob da3glob quadl for the last problem: 23
%%------------------------------------------------------------------------

%Problem_number =

%    23


%work =

%    25    25    25    18
%    41    41    41    18
%    53    53    53    18
%    73    73    73    78
%    93    93    93   108
%   113   113   113   168
%   137   137   137   198
%   189   169   169   288
%   217   201   201   348
%   281   225   225   468
%   337   249   257   648
%   401   305   337   888


%extra_digits =

%    1.3901    1.3901    1.3901    0.9433
%    1.8967    1.8967    1.8967   -0.0567
%    3.2256    3.2256    3.2256   -1.0567
%    2.7929    2.7929    2.7929   -0.8062
%    2.8205    2.8205    2.8205    2.2487
%    3.0642    2.1908    2.1922    2.0989
%    2.4460    2.0907    2.0991    2.5361
%    1.7868    2.4685    2.4322    3.3440
%    3.0292    2.6688    3.0895    4.3691
%    1.7278    2.3128    2.1502    3.3161
%    2.2257    1.9650    1.7304    4.5704
%    1.6744    2.1792    1.6605    3.5103
%%Final comments:
%%This last part illustrates the information available for each of the 23
%%problems if you activate the last statement in the main loop.
%%work: is the number of function evaluations used by the codes.
%%extra_digits: tol implies a certain number of correct decimals.
%%These numbers indicate if the code got more correct digits (positive) 
%%or less correct digits (negative) than we asked for through tol.
%%
%%Thus we see that for problem 23 quadl has three negative numbers in
%%extra_digits: therefore these three cases do not compete in the best code
%%competition: thus quadl is the best code when tol=0.1 only: see the
%%best_code and the failures reports above. 
%%
%%Note: the values given in extra_digits depend on the actual floating
%%point processor used: you may experience differences to these numbers when
%%you run this code on your machine.
%%End of comments.
%%-------------------------------------------------------------------
