% .
% TestSummary displays the data collected during the testing
% of the xGESVDQ part (SVD computation) of the SIGMA library.
% It uses data colected in the files
% -> SGESVDQ.table and DGESDVQ.table 
%    (after running real tests rtestt or rtests)
% -> CGESVDQ.table and ZGESVDQ.table 
%    (after running complex tests ctestt or ctests)
% If the appropriate file does not exist in the working directory,
% the program aborts with an error message.
% After running the corresponding test program, for each of the new
% SVD subroutines (SGESVDQ, DGESVDQ, CGESVDQ, ZGESVDQ), TestSummary
% displays the norms of the residuals (assuming all singular values and
% vectors are computed), errors in the singular values and vectors (if 
% computed).
% Part of testing devices of the SIGMA library.
% Coded by Zlatko Drmac, drmac@math.hr, University of Zagreb, Croatia.
%
subroutine = menu('Choose the subroutine being tested', ...
             'SGESVDQ', 'DGESVDQ', 'CGESVDQ', 'ZGESVDQ') ;

switch subroutine
    case 1,
        sbrtn = 'SGESVDQ'      ;
        load     SGESVDQ.table ;
        A = SGESVDQ ;
        eps0 = (eps(single(1/2))) ;
        ref = 'DGESVJ' ;
    case 2,
        sbrtn = 'DGESVDQ'      ;
        load     DGESVDQ.table ;
        A = DGESVDQ    ;
        eps0 = eps     ; 
        ref = 'DGESVJ' ;
    case 3,
        sbrtn = 'CGESVDQ'      ;
        load     CGESVDQ.table ;
        A = CGESVDQ ; 
        eps0 = (eps(single(1/2))) ;
        ref = 'ZGESVJ' ;
    case 4,
        sbrtn = 'ZGESVDQ'      ;
        load     ZGESVDQ.table ;
        A = ZGESVDQ ;   
        eps0 = eps ; 
        ref = 'ZGESVJ' ; 
end
%..........................................................................
load dimensions.data
m = dimensions(1,1) ;
n = dimensions(1,2) ; 
mexp  = size(A,1) ; 
neps  = n*eps0*ones(mexp,1) ; 
meps  = m*eps0*ones(mexp,1) ; 
ifig  = 0 ; 
%..........................................................................
% residuals of the computed SVD
if ( A(1,3) > -1 )
ifig = ifig + 1 ;
figure(ifig), semilogy(1:mexp,A(:,3),'b.-',1:mexp,meps,'r-')
legend(strcat(sbrtn,' residuals'), 'm*eps', 'Location','NorthEast')
xlabel('test case index')
ylabel(['$\|A-U\Sigma V^*\|_{F} / \|A\|_{F}$'],'Interpreter','Latex')
title('Computed residuals')
grid
end
%..........................................................................
% relative errors in the computed singular values and error bounds 
% |d\sigma_i| / \sigma_i (reference=ZGESVJ)
scond = A(:,1); 
ifig = ifig + 1 ;
figure(ifig), semilogy(1:mexp,A(:,4),'b.-',1:mexp,scond.*meps,'r-', ...
    1:mexp,scond.*neps/n,'m-.',1:mexp,meps,'-.')
legend(strcat(sbrtn,' rel. errors'), 'm*eps*scond', 'eps*scond','m*eps',...
                    'Location','NorthWest')
xlabel('test case index')
ylabel(['$\max_i |\delta \sigma_i|/\sigma_i$'],'Interpreter','Latex')
title(strcat('Relative errors in the singular values (reference= ',ref,')'))
grid
%..........................................................................
% relative errors in the computed singular values and error bounds
% |d\sigma_i| / \sigma_1  (reference=ZGESVJ)
ifig = ifig + 1 ;
figure(ifig), semilogy(1:mexp,A(:,13),'b.-', 1:mexp,meps,'r-')
legend(sbrtn, 'm*eps','Location','SouthEast')
xlabel('test case index')
ylabel(['$\max_i |\delta \sigma_i|/\sigma_1$'],'Interpreter','Latex')
title(strcat('Relative errors in the singular values (reference= ',ref,')'))
grid
%..........................................................................
% accuracy of the right singular vectors
% check compliance with perturbation theory
if ( A(1,11) > -1 )
ifig = ifig + 1 ;
figure(ifig), semilogy(1:mexp,A(:,11),'b.-',1:mexp,scond.*meps/m,'m-', ...
                    1:mexp,scond.*neps,'r-')
legend('max_i ||dv_i||*gap_i',  'eps*scond', 'n*eps*scond', ...
       'Location','NorthWest')
xlabel('test case index')
ylabel(['$\max_{i}\|\delta v_i\|_2\cdot gap_i$'],'Interpreter','Latex')
title(strcat('Accuracy of the right singular vectors',' (',sbrtn,')'))
grid
end
%..........................................................................
% accuracy of the left singular vectors
% check compliance with perturbation theory
if ( A(1,12) > -1 )
ifig = ifig + 1 ;    
figure(ifig), semilogy(1:mexp,A(:,12),'b.-',1:mexp,scond.*meps/m,'m-', ...
                    1:mexp,scond.*neps,'r-')
legend('max_i ||du_i||*gap_i',  'eps*scond', 'n*eps*scond', ...
       'Location','NorthWest')
xlabel('test case index')
ylabel(['$\max_{i}\|\delta u_i\|_2\cdot gap_i$'],'Interpreter','Latex')
title(strcat('Accuracy of the left singular vectors',' (',sbrtn,')'))
grid
end
...........................................................................
% orthogonality of the left singular vectors 
if ( A(1,5) > -1 )
ifig = ifig + 1 ;
figure(ifig), semilogy(1:mexp,A(:,5),'b.-',1:mexp,A(:,6),'g.-', ...
                    1:mexp,A(:,7),'m.-',1:mexp,meps,'r-')
legend('u_i^* u_j','min_i||u_i||', 'max_i||u_i||', 'm*eps', ...
       'Location','NorthEast')
xlabel('test case index')
ylabel(['$\max_{i\neq j} |u_i^* u_j|, \max_{i} |u_i^* u_i|,\min_{i}|u_i^* u_i|$'],...
       'Interpreter','Latex')
title(strcat('Orthogonality of the left singular vectors',' (',sbrtn,')'))
grid
end
%..........................................................................
% orthogonality of the right singular vectors
if ( A(1,8) > -1 )
ifig = ifig + 1 ; 
figure(ifig), semilogy(1:mexp, A(:,8), 'b.-', 1:mexp,A(:,9),'g.-', ...
                    1:mexp,A(:,10), 'm.-', 1:mexp,meps,'r-')
legend('v_i^* v_j', 'min_i||v_i||', 'max_i||v_i||', 'm*eps','Location','NorthEast')
xlabel('test case index')
ylabel(['$\max_{i\neq j} |v_i^* v_j|, \max_{i} |v_i^* v_i|,\min_{i}|v_i^* v_i|$'],...
       'Interpreter','Latex')
title(strcat('Orthogonality of the right singular vectors',' (',sbrtn,')'))
grid
end
%..........................................................................
%
if ( A(1,2) > -1 )
ifig = ifig + 1 ;
sf=n^(3/4) ;
%figure(6), semilogy(1:mexp,B(:,1),'b.-',1:mexp,sf*B(:,2),'go-')
figure(ifig), semilogy(1:mexp,(A(:,1))./(sf*A(:,2)),'b.-')
%legend('cond', 'estimate','Location','NorthEast')
xlabel('test case index')
ylabel(['$\kappa_{scaled}(A)/estimate$'],'Interpreter','Latex')
title('Scaled condition estimates')
grid
end

