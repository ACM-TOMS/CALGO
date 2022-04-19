function A = EigenmatInit

    % Order of A.
    n = 100;       
    
    % initialize A
    A = EigenmatGen(n, false, false);
    A.eig = 0.95.^(0:n-1)';
    A.type = ones(n,1);

    % Initialize Y
    A.Y.sig = rand(n,1);
    A.Y.u = rand(n,1) - 0.5;
    A.Y.v = rand(n,1) - 0.5;
    nu = sqrt(2)/norm(A.Y.u);
    A.Y.u = nu*A.Y.u;
    nu = sqrt(2)/norm(A.Y.v);
    A.Y.v = nu*A.Y.v;

    % Initialize Z
    A.Z.sig = rand(n,1);
    A.Z.u = rand(n,1) - 0.5;
    A.Z.v = rand(n,1) - 0.5;
    A.Z.nblocks = 2;
    A.Z.bs = [1;n/2+1;n+1];
    for i = 1:A.Z.nblocks
        i1 = abs(A.Z.bs(i));
        i2 = A.Z.bs(i+1)-1;
        if(i2>0)
            nu = sqrt(2)/norm(A.Z.u(i1:i2,1));
            A.Z.u(i1:i2,1) = nu*A.Z.u(i1:i2,1);
            nu = sqrt(2)/norm(A.Z.v(i1:i2,1));
            A.Z.v(i1:i2,1) = nu*A.Z.v(i1:i2,1);
        end
    end

