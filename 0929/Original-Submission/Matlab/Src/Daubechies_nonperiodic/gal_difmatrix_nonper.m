
function [Diff]=gal_difmatrix_nonper(D,N)
% [Diff]=gal_difmatrix_nonper(D,N)
% Differentiation projection matrix for wavelets defined on the interval
% Taken by the paper of ''The differentiation matrix for Daubechies-based
%wavelets on an interval'' by Leland Jameson.
% This function generates the Differentiation pojection matrix in the 
% non-periodic case of Galerkin approach.
% Dependencies
% cascade.m, conn.m, L_ro.m, R_ro.m, L_phi_origin.m, R_phi_origin.m,
% L_alpha.m, R_alpha.m

L=N;

%....................Grid x_1 and y_1.......................
for l=0:N-1
    x_1(l+1)=l+.5;
    y_1(l+1)=l*(N/(N-1));
end
M=D/2;
%...............M is the number of vanishing moment of the wavelet (D-1) is support.....
q=14;
[x1,phi] = cascade(D,q); % Get phi:		
h=x1(2)-x1(1);
con_1=conn(1,D);
supp=(D-1)*(1/h)+1;
L_sol=L_ro(D,L,h,supp,phi);
R_sol= R_ro(D,L,h,supp,phi);
 for l=0:N-1
    for k=0:N-1
        %--------------------First part of the row
        if(l>=0 & l<=M-1)
            %................First part of column.....
            if(k>=0 & k<=M-1)
                if(k==l)
                  [y,ro]=L_phi_origin(k,D,h,supp,phi) ;
                  D_j1(l+1,k+1)=ro;
                 else
                 ind=(l)*M+k;
                 D_j1(l+1,k+1)=L_sol(ind);
                 end
                  %................second part of column.....
              elseif(k>=M & k<=N-2*M+1)
               D_j1(l+1,k+1)=L_alpha(k,l,D,L);
            end
            
            %-------------------Second part of the row
    elseif(l>=M & l<=N-2*M+1)
        %................First part of column.....
            if(k>=0 & k<=M-1)
                D_j1(l+1,k+1)=-L_alpha(l,k,D,L); % l,k is due to definition
                 %................second part of column.....
            elseif(k>=M & k<=(N-2*M+1))
                ind=l-k+(D-1);
                   if(ind>=1 & ind<=2*(D-2)+1)
                   D_j1(l+1,k+1)=con_1(ind);
                   end
                    %................Third part of column.....
             elseif(k>=(N-2*M+2) & k<=N-1)
                 if(-(l-(N-2*M+1) )>=0 & -(l-(N-2*M+1))<=M-1 & -(k-(N-1))>=0 & -(k-(N-1))<=M-1)
                  D_j1(l+1,k+1)=-R_alpha(-(l-(N-1)),-(k-(N-1)),D,L);
                  else
                  D_j1(l+1,k+1)=0;
                  end
              end            
             %-------------------Third part of the row
        elseif(l>=(N-2*M+2) & l<=N-1)
            %................second part of column.....
            if(k>=M & k<=N-2*M+1)
              
                if(-(l-(N-1) )>=0 & -(l-(N-1))<=M-1 & -(k-(N-2*M+1))>=0 & -(k-(N-2*M+1))<=M-1)
                 D_j1(l+1,k+1)=R_alpha(-(k-(N-1)),-(l-(N-1)),D,L);
                 else
                  D_j1(l+1,k+1)=0;
                  end
                  %................First part of column.....
             elseif(k>=(N-2*M+2) & k<=N-1)
                if(k==l)
                     [y,ro]=R_phi_origin(-(k-(N-1)),D,h,supp,phi) ;
                      D_j1(l+1,k+1)=ro;
                else
                ind=(l-(N-2*M+2))*(2*M-2)+(k-(N-2*M+2));
                if(ind<=length(R_sol)) %dbm
                D_j1(l+1,k+1)=R_sol(ind);
                end
                end
           end
       end
    end
end
Diff=D_j1;
