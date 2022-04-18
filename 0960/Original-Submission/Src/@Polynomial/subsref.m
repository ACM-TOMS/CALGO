% SUBSREF Implementing the following syntax:
% obj([1 ...])
% obj.BernsCoeff
% obj.plot
% out = obj.method(args)
% out = obj.method
function b = subsref(a,s)  
    if isempty(s)
        error('Polynomial.subsref: missing index');
    end
    
    switch s(1).type
        case '()'
            if length(s.subs)>1
                x = s.subs{1};
                prec = s.subs{2};
                b = a.Eval(x,prec);
            else
                x = s.subs{:};
                b = a.Eval(x);
            end
        case '{}'
            error('Polynomial.subsref: cell arrays not implemented for polynomials');
        case '.'             
             switch s(1).subs
                 case 'BernsCoeff'
                     error('Polynomial.BernsCoeff: it is a private member, uses getCoeff or double methods instead');
                 case 'disp'
                     a.disp;
                 case 'plot'
                     if length(s)>1
                         a.plot(s(2).subs{:});
                     else
                         a.plot;
                     end
                 case 'Eval'
                     if length(s)>1 & ~isempty(s(2).subs)
                         if length(s(2).subs)==1
                             b = Eval(a,s(2).subs{1});
                         elseif length(s(2).subs)==2
                             b = a.Eval(s(2).subs{1},s(2).subs{2});
                         else
                             error('Polynomial.Eval: Too many input argumets');
                         end
                     else
                         error('Polynomial.Eval: Input argumets are missed');
                     end
                 case 'char'
                     b = a.char;
                 otherwise
                     if length(s)>1
                         b = a.(s(1).subs)(s(2).subs{:});
                     else
                         b = a.(s.subs);
                     end
             end
    end
end


