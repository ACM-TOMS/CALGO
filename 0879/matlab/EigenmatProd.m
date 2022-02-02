function C = EigenmatProd(A, B, shift, job)
%
%  C = EigenmatProd(A, B, shift, job)
%
%  computes product, perhaps shifted and inverted, of the
%  eigenmat A and a matrix B, placing the result in C.
%
%  A       The eigenmat
%  B       The array containing the matrix B.
%  shift   A shift.
%  job     A string specifying the operation to be performed.
%
%          'ab'    C = (A - shift*I)*B
%          'atb'   C = (A - shift*I)'*B
%          'aib'   C = (A - shift*I)\B
%          'aitb'  C = (A - shift*I)'\B
%
%  EigenmatProd useds HsvdProd to compute hsvdmat products.


    n = A.n;          % Get the order of the matrix.
    
    C = B;
    if (strcmp(job,'ab'))

        % Compute C = (A-shift*I)*B  = Y*(Z*((L-shift*I)*(Z\(Y\C)))).

        % Compute C = Z\(Y\C).

        C = HsvdProd(A.Y, C, 'aib');
        C = HsvdProd(A.Z, C, 'aib');

        % Compute (L-shift*I)*C.  The index i points to successive blocks
        % of L.

        i = 1;
        while (i<=n)

            if(A.type(i) == 1)

                % Real eigenvalue.

                C(i,:) = (A.eig(i)-shift)*C(i,:);
                i = i+1;

            elseif(A.type(i)==2 && i==n)
                error('Error in EigenmatProd: 2x2 block starts at eig(n).');

            elseif(A.type(i) == 2 && A.type(i+1) == 3)

                % Complex eigenvalue.

                mu = A.eig(i)-shift;
                nu = A.eig(i+1);
                temp = mu*C(i,:) + nu*C(i+1,:);
                C(i+1,:) = -nu*C(i,:) + mu*C(i+1,:);
                C(i,:) = temp;
                i = i+2;

            elseif(A.type(i) < -1 )

                % Jordan block.
                % from i to j
                j = i+abs(A.type(i))-1;
                if (j>A.n)
                   error(strcat('Error in EigenmatProd:', ... 
                            'Jordan block too large.'));
                end 
                if (any(A.type(i+1:j)~=-1))
                    error(strcat('Error in EigenmatProd:', ... 
                            'Illegal type for Jordan block.'));
                end
                mu = A.eig(i) - shift;
                for k=1:size(C,2)
                    C(i:j-1,k) = mu*C(i:j-1,k) + A.eig(i+1:j).*C(i+1:j,k);
                end
                C(j,:) = mu*C(j,:);
                
                i = j+1;
                
            else
                error('Error in EigenmatProd: Illegal type.');
            end
        end

        % Compute C = Y*(Z*C).

        C = HsvdProd(A.Z, C, 'ab');
        C = HsvdProd(A.Y, C, 'ab');

    elseif(strcmp(job,'atb'))

        % Compute C = (A-shift*I)'*C = Y'\(Z'\(((L-shift*I)'*(Z'*(Y'*C)))).

        % Compute C = Z'*(Y'*C).

        C = HsvdProd(A.Y, C, 'atb');
        C = HsvdProd(A.Z, C, 'atb');

        % Compute (L-shift*I)*C.  The index i points to successive blocks
        % of L.

        i = 1;
        while (i<=n)

            if(A.type(i) == 1)

                % Real eigenvalue.

                C(i,:) = (A.eig(i)-shift)*C(i,:);
                i = i+1;

            elseif(i == n)
                error('Error in EigenmatProd: 2x2 block starts at eig(n).');

            elseif(A.type(i) == 2 && A.type(i+1) == 3)

                % Complex eigenvalue

                mu = A.eig(i) - shift;
                nu = A.eig(i+1);
                temp = mu*C(i,:) - nu*C(i+1,:);
                C(i+1,:) = nu*C(i,:) + mu*C(i+1,:);
                C(i,:) = temp;
                i = i+2;

            elseif(A.type(i) < -1 )

                % Jordan block.
                % from i to j
                j = i+abs(A.type(i))-1;
                if (j>A.n)
                   error(strcat('Error in EigenmatProd:', ... 
                            'Jordan block too large.'));
                end 
                if (any(A.type(i+1:j)~=-1))
                    error(strcat('Error in EigenmatProd:', ... 
                            'Illegal type for Jordan block.'));
                end
                mu = A.eig(i) - shift;
                for k=1:size(C,2)
                    C(i+1:j,k) = mu*C(i+1:j,k)+A.eig(i+1:j).*C(i:j-1,k);
                end
                C(i,:) = mu*C(i,:);
                
                i = j+1;
                
            else
                error('Error in EigenmatProd: Illegal type.');
            end
        end

        % Compute C = Y'\(Z'\C).

        C = HsvdProd(A.Z, C, 'aitb');
        C = HsvdProd(A.Y, C, 'aitb');

    elseif (strcmp(job,'aib'))

        % Compute C = (A-shift*I)*C = Y*(Z*((L-shift*I)\(Z\(Y\C)))).

        % Compute C = Z\*(Y\C).

        C = HsvdProd(A.Y, C, 'aib');
        C = HsvdProd(A.Z, C, 'aib');

        % Compute (L-shift*I)\C.  The index i points to successive blocks
        % of L.

        i = 1;
        while (i<=n)
            if(A.type(i) == 1)

                % Real eigenvalue.

                C(i,:) = C(i,:)./(A.eig(i)-shift);
                i = i+1;

            elseif(i == n)
                error('Error in EigenmatProd: 2x2 block starts at eig(n).');

            elseif(A.type(i) == 2 && A.type(i+1) == 3)

                % Complex eigenvalue.

                mu = A.eig(i) - shift;
                nu = A.eig(i+1);
                Li = [mu, nu; -nu, mu];
                C(i:i+1,:) = Li\C(i:i+1,:);

                i = i+2;

            elseif(A.type(i) < -1 )

                % Jordan block.
                % from i to j
                j = i+abs(A.type(i))-1;
                if (j>A.n)
                   error(strcat('Error in EigenmatProd:', ... 
                            'Jordan block too large.'));
                end 
                if (any(A.type(i+1:j)~=-1))
                    error(strcat('Error in EigenmatProd:', ... 
                            'Illegal type for Jordan block.'));
                end

                mu = A.eig(i) - shift;
                
                C(j,:) = C(j,:)/mu;
                for k=j-1:-1:i
                    C(k,:) = (C(k,:)-C(k+1,:)*A.eig(k+1))/mu;
                end
                
                i = j+1;
                
            else
                error('Error in EigenmatProd: Illegal type.');
            end
        end

        % Compute C = Y*(Z*C).

        C = HsvdProd(A.Z, C, 'ab');
        C = HsvdProd(A.Y, C, 'ab');

    elseif(strcmp(job,'aitb'))

        % Compute C = A^-T*C = Y'\(Z'\((L-shift*I)'\(Z'*(Y'*C)))).

        % Compute C = Z'*(Y'*C).

        C = HsvdProd(A.Y, C, 'atb');
        C = HsvdProd(A.Z, C, 'atb');

        % Compute (L-shift*I)'\C.  The index i points to successive blocks
        % of L.

        i = 1;
        while (i<=n)
            if(A.type(i) == 1)

                % Real eigenvalue.

                C(i,:) = C(i,:)./(A.eig(i)-shift);
                i = i+1;

            elseif(i == n)
                error('Error in EigenmatProd: 2x2 block starts at eig(n).')

            elseif(A.type(i) == 2 && A.type(i+1) == 3)

                % Complex eigenvalue

                mu = A.eig(i) - shift;
                nu = A.eig(i+1);
                Li = [mu, -nu; nu, mu];
                C(i:i+1,:) = Li\C(i:i+1,:);
                i = i+2;
            elseif(A.type(i) < -1 )

                % Jordan block.
                % from i to j
                j = i+abs(A.type(i))-1;
                if (j>A.n)
                   error(strcat('Error in EigenmatProd:', ... 
                            'Jordan block too large.'));
                end 
                if (any(A.type(i+1:j)~=-1))
                    error(strcat('Error in EigenmatProd:', ... 
                            'Illegal type for Jordan block.'));
                end

                mu = A.eig(i) - shift;
                
                C(i,:) = C(i,:)/mu;
                for k=i+1:j
                    C(k,:) = (C(k,:)-C(k-1,:)*A.eig(k))/mu;
                end
                
                i = j+1;
            else
                error('Error in EigenmatProd: Illegal type.');
            end
        end

        % Compute C = Y'\(Z'\C).

        C = HsvdProd(A.Z, C, 'aitb');
        C = HsvdProd(A.Y, C, 'aitb');

    else
        error('Error in EigenmatProd: Illegal operation.')
    end

    return
