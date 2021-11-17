% IIPBF_plots runs the 18 cases used to evaluate the performance and
% accuracy of IIPBF. The script's output includes plots of Actual error vs
% function evaluations, and Estimated error vs. function evaluations. For
% details see reference paper.
% 
%   Author: Sirong Zhang & Jung Hun Kim
% 
%   Reference: Ratnanather, J. T., Kim, J. H., Zhang, S., Davis, A. M. J., 
%              and Lucas, S. K. 2012. IIPBF: a MATLAB toolbox for infinite 
%              integrals of product of Bessel functions. Pending review
%              from ACM Transactions on Mathematical Software
% 
%   Revision Date: Aug  31, 2012 
%   Author: 10/08/2010 Jung Hun Kim

path = path(path,'..');

clear all
clc
close all
[res01,abs01,rel01,n01,~]=testcases(1);
[res02,abs02,rel02,n02,~]=testcases(2);
[res03,abs03,rel03,n03,~]=testcases(3);
[res04,abs04,rel04,n04,~]=testcases(4);
[res05,abs05,rel05,n05,~]=testcases(5);
[res06,abs06,rel06,n06,~]=testcases(6,1,2,1);
[res07,abs07,rel07,n07,~]=testcases(7,1,2,1);
[res08,abs08,rel08,n08,~]=testcases(8,1.5);
[res09,abs09,rel09,n09,~]=testcases(9,.2);
[res10,abs10,rel10,n10,~]=testcases(10,2);

[res11,abs11,rel11,n11,~]=testcases(11,1,1,3);
[res12,abs12,rel12,n12,~]=testcases(12,1,1,3);
[res13,abs13,rel13,n13,~]=testcases(13,1,1,3);

[res14,abs14,rel14,n14,~]=testcases(14,3,1);
[res15,abs15,rel15,n15,~]=testcases(15,3,1);

[res16,abs16,rel16,n16,~]=testcases(16,3,1);
[res17,abs17,rel17,n17,~]=testcases(17,1,1,1);
[res18,abs18,rel18,n18,~]=testcases(18,1,3,1);

%Optional
%[res19,abs19,rel19,n19,~]=testcases(19,1,2);
%[res20,abs20,rel20,n20,~]=testcases(20,1,3);
%[res21,abs21,rel21,n21,~]=testcases(21,1,2);
%[res22,abs22,rel22,n22,~]=testcases(22,2);
%[res23,abs23,rel23,n23,~]=testcases(23,1,2,1,2);

figure

semilogy (n01, abs01, 'ko--','Linewidth',3)
hold on
semilogy (n02, abs02, 'mo:','Linewidth',2)
semilogy (n03, abs03, 'm*-')
semilogy (n04, abs04, 'r.- ')
semilogy (n05, abs05, 'gd- ','Linewidth',3)
semilogy (n06, abs06, 'ks-','Markersize',10)
semilogy (n07, abs07, 'b-','Linewidth',3)
semilogy (n08, abs08, 'bx--','Markersize',12)
semilogy (n09, abs09, 'rs-','Linewidth',3)
semilogy (n10, abs10,'kd-.','Markersize',12)
semilogy (n11, abs11, 'k-.','Linewidth',3)
semilogy (n12, abs12,'bs-')
semilogy (n13, abs13,'b+-','Linewidth',3,'Markersize',12);
semilogy (n14, abs14,'rx-','Markersize',10);
semilogy (n15, abs15,'gx-','Linewidth',2);
semilogy (n16, abs16,'r--','Linewidth',2);
semilogy (n17, abs17, 'yd-','Linewidth',3)
semilogy (n18, abs18, 'yd-')
set(gca,'Fontsize',14,'Fontname','Helvetica','Linewidth',3,'Xlim',[0 120],'Xtick',[0 20 40 60 80 100],'Ylim',[0 1e-8],'Ytick',[1e-22 1e-20 1e-18 1e-16 1e-14 1e-12 1e-10 1e-8])
xlabel('Number of Evaluations','Fontsize',14,'Fontname','Helvetica')
ylabel('Actual Error','Fontsize',14,'Fontname','Helvetica')
legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','Location','NorthEast')
print -dpng Fig1.png

hold off

figure

semilogy (n01, rel01, 'ko--','Linewidth',3)
hold on
semilogy (n02, rel02, 'mo:','Linewidth',2)
semilogy (n03, rel03, 'm*-')
semilogy (n04, rel04, 'r.- ')
semilogy (n05, rel05, 'gd- ','Linewidth',3)
semilogy (n06, rel06, 'ks-','Markersize',10)
semilogy (n07, rel07, 'b-','Linewidth',3)
semilogy (n08, rel08, 'bx--','Markersize',12)
semilogy (n09, rel09, 'rs-','Linewidth',3)
semilogy (n10, rel10,'kd-.','Markersize',12)
semilogy (n11, rel11, 'k-.','Linewidth',3)
semilogy (n12, rel12,'bs-')
semilogy (n13, rel13,'b+-','Linewidth',3,'Markersize',12);
semilogy (n14, rel14,'rx-','Markersize',10);
semilogy (n15, rel15,'gx-','Linewidth',2);
semilogy (n16, rel16,'r--','Linewidth',2);
semilogy (n17, rel17, 'yd-','Linewidth',3)
semilogy (n18, rel18, 'yd-')
set(gca,'Fontsize',14,'Fontname','Helvetica','Linewidth',3,'Xlim',[0 120],...
    'Xtick',[0 20 40 60 80 100],'Ylim',[0 1e-8],'Ytick',[1e-16 1e-14 1e-12 1e-10 1e-8 1e-4])
xlabel('Number of Evaluations','Fontsize',14,'Fontname','Helvetica')
ylabel('Estimated Error','Fontsize',14,'Fontname','Helvetica')
legend('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','Location','NorthEast')
print -dpng Fig2.png

hold off
