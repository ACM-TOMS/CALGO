function [varargout] = subsasgn(poly, S, value)


if (strcmp(S.type,'()'))
    if (iscell(S.subs{1}))
        S.subs = S.subs{1};
    end
    if (length(S.subs) ~= length(poly.vertices))
        error('Index exceeds polynomial dimensions.');
        return
    end
    
    jump = 1;
    numelem = 1;
    colon = [];
    expon = []; %Which exponents to search
    totalsimpl = [];
    for cont=1:length(poly.vertices)
        if (cont > 1)
            jump(cont) = numelem;
        end
        deg = sum(poly.data(1).exponent{cont});
        N = poly.vertices(cont);
        if (N > 0)
            totalsimpl(cont) = (factorial(N+deg-1)/(factorial(deg)*factorial(N-1)));
            numelem = numelem*totalsimpl(cont);
        end
        if ((isstr(S.subs{cont})) && (strcmp(S.subs{cont},':')))
            colon = [colon cont];
        end
    end
    
    ind = search(poly,S.subs,jump,colon);
    for cont=1:length(ind)
        poly.data(ind(cont)).value = value;
    end
    varargout{1} = poly;
end

return


function index = search(poly,exponent,jump,colon)
    
  noncolon = setxor(colon,1:length(poly.vertices));
  index = 1:length(poly.data);
  for cont = 1:length(noncolon)
      indexini = 1;
      while (sum(abs(exponent{noncolon(cont)} - poly.data(indexini).exponent{noncolon(cont)})) > 0)
          indexini = indexini + jump(noncolon(cont));
          if (indexini > length(poly.data))
              error('rolmip:idn','Index does not exist');
          end
      end
      aux = [];
      while (indexini <= length(poly.data))
          aux = [aux indexini:indexini + jump(noncolon(cont)) - 1];
          if (noncolon(cont) < length(jump))
              indexini = indexini + jump(noncolon(cont)+1);
          else
              indexini = length(poly.data) + 1;
          end
      end
          
      index = intersect(index,aux);
  end

return