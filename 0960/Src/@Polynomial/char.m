      % CHAR    CHAR(OBJ) creates a formated display of the polynom as
      % linear combination of Bernstein polynomials
      function str = char(obj) 
          if all(obj.BernsCoeff == 0)
              str = ['0'];
          else
              grado = obj.getDegree();
              n = grado;
              s = cell(1,n);
              ind = 1;
              for a = obj.BernsCoeff;
                  if a ~= 0;
                      if ind ~= 1
                          if a > 0
                              s(ind) = {' + '};
                              ind = ind + 1;
                          else
                              s(ind) = {' - '};
                              a = -a; %#ok<FXSET>
                              ind = ind + 1;
                          end
                      end
                      if a ~= 1 || grado == 0
                          if a == -1
                              s(ind) = {'-'};
                              ind = ind + 1;                            
                          else                              
                              s(ind) = {num2str(a)};
                              ind = ind + 1;
                              if grado~=0
                                  s(ind) = {'*'};
                                  ind = ind + 1;
                              end
                          end
                      end
                      switch grado
                          case 0
                              % Not to do anything
                          case 1
                              if n==1
                                  s(ind)={['bin(' int2str(grado) ',' int2str(grado-n) ')*(1-x)']};
                              else
                                  s(ind)={['bin(' int2str(grado) ',' int2str(grado-n) ')*x']};
                              end
                              ind = ind + 1;
                          case 2
                              if n==2
                                  s(ind)={['bin(' int2str(grado) ',' int2str(grado-n) ')*(1-x)^2']};
                              elseif n==1
                                  s(ind)={['bin(' int2str(grado) ',' int2str(grado-n) ')*x*(1-x)']};
                              else
                                  s(ind)={['bin(' int2str(grado) ',' int2str(grado-n) ')*x^2']};
                              end
                              ind = ind + 1;
                          otherwise
                              if n==grado
                                  s(ind)={['bin(' int2str(grado) ',' int2str(grado-n) ')*(1-x)^' int2str(n)]};
                              elseif n==grado-1
                                  s(ind) = {['bin(' int2str(grado) ',' int2str(grado-n) ')*x(1-x)^' int2str(n)]};
                              elseif n==1
                                  s(ind) = {['bin(' int2str(grado) ',' int2str(grado-n) ')*x^' int2str(grado-n) '(1-x)']};
                              elseif n==0
                                  s(ind) = {['bin(' int2str(grado) ',' int2str(grado-n) ')*x^' int2str(grado-n)]};
                              else
                                  s(ind) = {['bin(' int2str(grado) ',' int2str(grado-n) ')*x^' int2str(grado-n) '(1-x)^' int2str(n)]};
                              end
                              ind = ind + 1;
                      end
                  end
                  n = n - 1;
              end
              str = [s{:}];
          end
      end 

