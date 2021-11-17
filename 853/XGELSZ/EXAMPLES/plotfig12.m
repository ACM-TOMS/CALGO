% plotfig12.m
% This Matlab script file plots the run times of DGELSZ, DGELSY,
% and DGELSD versus the matrix rank.  It produces figures similar
% to those in Figures 1 and 2 in the paper by Foster and Kommu.
% Before using plotfig12 run the Fortran program cfig12.f, dfig12.f
% sfig12.f or zfig12.f which produce the output file fig12.out. 
% Following this  run plotfig12.m from inside Matlab. If there are 
% runs with identical ranks then the times for these ranks are
% averaged before being plotted in figure(1) and all the times 
% are printed in figure(2).  For the plots to be meaningful the
% number of rows and the number of columns should be the same for
% every matrix.

% This code is part of a package for solving rank deficient least
% squares problems, written by:
% ==================================================================
% L. Foster                   and   R. Kommu
% Department of Mathematics         Department of Physics
% San Jose State University         San Jose State University
% San Jose, CA 95192                San Jose, CA 95192
% foster@math.sjsu.edu              rkommu@email.sjsu.edu
% ==================================================================
% 03/05/2004
%

load fig12.out

% determine distinct ranks and sort in increasing order
ranks = fig12(:,3);
[ranks,sort_index] = sort(ranks);
rankchange = diff(ranks);
rankdistinct = ranks(find(rankchange > 0)+1);
rankdistinct = [ranks(1); rankdistinct];
ndistinct = length(rankdistinct);

timed = fig12(:,4);
timey = fig12(:,5);
timez = fig12(:,6);

timed = timed(sort_index);
timey = timey(sort_index);
timez = timez(sort_index);

% average the times corresponding to the same rank
timed_ave = zeros(ndistinct,1);
timey_ave = zeros(ndistinct,1);
timez_ave = zeros(ndistinct,1);
timed_std = zeros(ndistinct,1);
timey_std = zeros(ndistinct,1);
timez_std = zeros(ndistinct,1);
timed_std_mean = zeros(ndistinct,1);
timey_std_mean = zeros(ndistinct,1);
timez_std_mean = zeros(ndistinct,1);
timez_y_std = zeros(ndistinct,1);

for i = 1:ndistinct
   rank_same_index = find( ranks==rankdistinct(i) );
   timed_same_rank = timed( rank_same_index );
   timey_same_rank = timey( rank_same_index );
   timez_same_rank = timez( rank_same_index );
   timed_ave(i) = mean(timed_same_rank);
   timey_ave(i) = mean(timey_same_rank);
   timez_ave(i) = mean(timez_same_rank);
   timed_std(i) = std( timed_same_rank );
   timey_std(i) = std( timey_same_rank );
   timez_std(i) = std( timez_same_rank );
   timed_std_mean(i) = timed_std(i) / sqrt(length(timed_same_rank));
   timey_std_mean(i) = timey_std(i) / sqrt(length(timey_same_rank));
   timez_std_mean(i) = timez_std(i) / sqrt(length(timez_same_rank));
   timez_y_std(i) = std(timez_same_rank-timey_same_rank)/sqrt(length(timez_same_rank));
end
   
%draw the plots
   
figure(2)
clf
% plot all the times
plot(ranks,timez,'x',ranks,timey,'o',ranks,timed,'+'),grid,shg
plot(ranks,timez,'x',ranks,timey,'o',ranks,timed,'+', ...
   rankdistinct,timez_ave,'-',rankdistinct,timey_ave,'-',  ...
   rankdistinct,timed_ave,'-'),grid,shg
ylabel('time in seconds')
xlabel('numerical rank')
legend('xGELSZ','xGELSY','xGELSD',2)
title('run times (solid lines are mean times)')

figure(1)
clf
% plot only the average times at each rank
plot(rankdistinct,timez_ave,'x-',rankdistinct,timey_ave,'o-',  ...
   rankdistinct,timed_ave,'+-'),grid,shg
ylabel('time in seconds')
xlabel('numerical rank')
legend('xGELSZ','xGELSY','xGELSD',2)
title('mean run times')

shg

