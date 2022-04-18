function testeigentest()
%
%  This program tests the package Eigentest.  It runs 64 testcases
%  probing various aspects of the package, as described below.  The
%  numbers in the output should be within about two orders of
%  magnitude of the rounding unit.

   % These values may be adjusted at compile time to change the
   % tests, which involve computing (A - shift*I) op B,
   % where op is *, '*,  \, or '\ (' denoting transpose).

   n = 10;       % Order of A.
   ncols = 3;    % Number of columns in B.
   shift = -1;   % A shift.
   
   % Loop on test case number. 

   for cs = 0:63
       disp(strcat('test case:',num2str(cs)))
       
       % Choose the type of the outer hsvd Y.  It can be an identity
       % (atype == 0) or a random hsvd (ytype == 1).       

       ytype = mod(cs,2);

       % Choose the typep of the inner hsvd Z.
       %
       %   If ztype == 0, Z is an identity.
       %   If ztype == 1, Z has two blocks: (1:6)(7:10).
       %   The first block is identity.
       %   If ztype == 2, Z has three blocks: (1:2)(3:8)(9:10),
       %     and the second block is an indentity.
       %   If ztype == 3, Z has two blocks: (1:4)(5::10),
       %     and the second block is identity.

       ztype = mod(floor(cs/2), 4);
       bs = [1 2 3 2];
       
       % Chose the types of eigenvalues.
       %
       %   If atype == 0, all eigenvalues are real.
       %   If atype == 1, all eigenvalues are complex.
       %   If atype == 2, eigenvalues 1,2,5,6,9,10 are complex.
       %   If atype == 3, eigenvalues 3,4,7,8 are complex.
       %   If atype == 4, all eigenvalues are in a Jordan block.
       %   If atype == 5, eigenvalues (3,4),(7,8,9) are in Jordan blocks.
       %   If atype == 6, eigenvalues (1-3),(8-10) are in a Jordan block.
       %   If atype == 7, eigenvalues (2-4) are in a Jordan block.
       %                  eigenvalues 6,7,8,9 are complex
       atype = mod(floor(cs/8), 8);       
       
       % Initialize A.

       A = EigenmatGen(n, bs(ztype+1), ytype == 0, ztype == 0);

       % Generate eigenvalues.

       A.eig = rand(n,1) + 0.5;
       switch(atype)
           case 0
               A.type = ones(1,n);

           case 1
               A.type(1:2:n) = 2;
               A.type(2:2:n) = 3;

           case 2
               A.type([3,4,7,8]) = 1;
               A.type([1,5,9]) = 2;
               A.type([2,6,10]) = 3;

           case 3
               A.type([1,2,5,6,9,10]) = 1;
               A.type([3,7]) = 2;
               A.type([4,8]) = 3;
               
           case 4
               A.type(1)   = -10;
               A.type(2:10) = -1;

           case 5
               A.type(3)   = -2;
               A.type(7)   = -3;
               A.type([4,8,9]) = -1;
               A.type([1,2,5,6,10]) = 1;

           case 6
               A.type([1,8]) = -3;
               A.type(2:3)   = -1;
               A.type(9:10)  = -1;
               A.type(4:7)   =  1;
               
           case 7
               A.type(2)     = -3;
               A.type(3:4)    = -1;
               A.type([1,5,10])  =  1;
               A.type([6,8]) = 2;
               A.type([7,9]) = 3;
               
       end  

       disp(strcat('A type = ',num2str(A.type)))

       % Generate Y

       if(ytype == 1)
           A.Y.sig = rand(n,1) + 1.0;
           A.Y.u = hscal(rand(n,1) - 0.5);
           A.Y.v = hscal(rand(n,1) - 0.5);
       end

       disp(strcat('Y bs = ',num2str(A.Y.bs)))

       % Generate Z.

       if(ztype >0)
           A.Z.sig = rand(n,1) + 1.0;
           switch(ztype)
               case 1
                   % If ztype == 1, Z has two blocks: (1:6)(7:10).
                   % The first block is identity.

                   A.Z.bs = [1,-7,n+1];

               case 2
                   % If ztype == 2, Z has three blocks: (1:2)(3:8)(9:10).
                   % The second block is indentity.

                   A.Z.bs = [1 3 -9 n+1];

               case 3
                   % If ztype == 3, Z has two blocks: (1:4)(5::10).
                   % The second block is identity.

                   A.Z.bs = [1 5 -n-1];
           end
           for i = 1:A.Z.nblocks
               i1 = abs(A.Z.bs(i));
               i2 = A.Z.bs(i+1)-1;
               if(i2>0)
                   A.Z.u(i1:i2,1) = hscal(rand(i2-i1+1,1) - 0.5);
                   A.Z.v(i1:i2,1) = hscal(rand(i2-i1+1,1) - 0.5);
               end
           end
       end

       disp(strcat('Z bs = ',num2str(A.Z.bs)))

       % Generate B for A, Y, and Z to operate on.

       b = rand(n,ncols) - 0.5;
       
       % Test Z.

       c = HsvdProd(A.Z, b, 'ab');
       c = HsvdProd(A.Z, c, 'aib');
       disp(strcat('|Z\Zb-b|   =',num2str(norm(b-c,1))))

       c = HsvdProd(A.Z, b, 'aib');
       c = HsvdProd(A.Z, c, 'ab');
       disp(strcat('|ZZ\b-b|   =',num2str(norm(b-c,1))))

       c = HsvdProd(A.Z, b, 'atb');
       c = HsvdProd(A.Z, c, 'aitb');
       disp(strcat('|Zt\Ztb-b| =',num2str(norm(b-c,1))))

       c = HsvdProd(A.Z, b, 'aitb');
       c = HsvdProd(A.Z, c, 'atb');
       disp(strcat('|ZtZt\b-b| =',num2str(norm(b-c,1))))

       % Test Y.

       c = HsvdProd(A.Y, b, 'ab');
       c = HsvdProd(A.Y, c, 'aib');
       disp(strcat('|Y\Yb-b|   =',num2str(norm(b-c,1))))

       c = HsvdProd(A.Y, b, 'aib');
       c = HsvdProd(A.Y, c, 'ab');
       disp(strcat('|YY\b-b|   =',num2str(norm(b-c,1))))

       c = HsvdProd(A.Y, b, 'atb');
       c = HsvdProd(A.Y, c, 'aitb');
       disp(strcat('|Yt\Ytb-b| =',num2str(norm(b-c,1))))

       c = HsvdProd(A.Y, b, 'aitb');
       c = HsvdProd(A.Y, c, 'atb');
       disp(strcat('|YtYt\b-b| =',num2str(norm(b-c,1))))

       % Test A.

       c = EigenmatProd(A, b, shift, 'ab');
       c = EigenmatProd(A, c, shift, 'aib');
       disp(strcat('|A\Ab-b|   =',num2str(norm(b-c,1))))

       c = EigenmatProd(A, b, shift, 'aib');
       c = EigenmatProd(A, c, shift, 'ab');
       disp(strcat('|AA\b-b|   =',num2str(norm(b-c,1))))

       c = EigenmatProd(A, b, shift, 'atb');
       d = EigenmatProd(A, c, shift, 'aitb');
       disp(strcat('|At\Atb-b| =',num2str(norm(b-d,1))))

       c = EigenmatProd(A, b, shift, 'aitb');
       d = EigenmatProd(A, c, shift, 'atb');
       disp(strcat('|AtAt\b-b| =',num2str(norm(b-d,1))))
       
       % Test the eigenvector calculations by computing all
       % left and right eigenvectors and their residuals.

       i = 1;
       ev1 = zeros(n-1,1);
       ev  = zeros(n,1);
       while (i<=n)
           [ev(i) x(1:n,i) y(1:n,i)] = EigenmatVecs(A, i, 'b');
           if(A.type(i) == 1 || A.type(i) < -1)
               i = i + 1;
           elseif(A.type(i) == -1)
               ev1(i-1) = A.eig(i); 
               i = i + 1; 
           elseif(A.type(i) == 2)
               x(1:n,i+1) = conj(x(1:n,i));
               y(1:n,i+1) = conj(y(1:n,i));
               ev(i+1) = conj(ev(i));
               i = i + 2;
           end  
       end    
  
       % Right eigensystem. 
       ax = EigenmatProd(A, x, 0, 'ab');
       E  = diag(ev);
       E1 = diag(ev1, 1);
       disp(strcat('|AX-E|     =',num2str(norm(ax-x*E-x*E1))))
       
       % Left eigensystem.

       ay = EigenmatProd(A, y, 0, 'atb');
       disp(strcat('|AtY-YEt|  =',num2str(norm(ay-y*conj(E)-y*E1'))))
   end

