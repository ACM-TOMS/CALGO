with Complex_Numbers;                use Complex_Numbers;

package body Permute_Operations is

  function "*" ( p : Permutation; v : Natural_Vectors.Vector )
	       return Natural_Vectors.Vector is

    r : Natural_Vectors.Vector(v'range);

  begin
    for i in p'range loop
      if p(i) >= 0
       then r(i) := v(p(i));
       else r(i) := -v(-p(i));
      end if;
    end loop;
    return r;
  end "*";

  function "*" ( p : Permutation; v : Integer_Vectors.Vector )
	       return Integer_Vectors.Vector is

    r : Integer_Vectors.Vector(v'range);

  begin
    for i in p'range loop
      if p(i) >= 0
       then r(i) := v(p(i));
       else r(i) := -v(-p(i));
      end if;
    end loop;
    return r;
  end "*";

  function "*" ( p : Permutation; v : Float_Vectors.Vector )
               return Float_Vectors.Vector is

    r : Float_Vectors.Vector(v'range);

  begin
    for i in p'range loop
      if p(i) >= 0
       then r(i) := v(p(i));
       else r(i) := -v(-p(i));
      end if;
    end loop;
    return r;
  end "*";

  function "*" ( p : Permutation; v : Complex_Vectors.Vector )
               return Complex_Vectors.Vector is

    r : Complex_Vectors.Vector(v'range);

  begin
    for i in p'range loop
      if p(i) >= 0
       then r(i) := v(p(i));
       else r(i) := -v(-p(i));
      end if;
    end loop;
    return r;
  end "*";

  function Permutable ( v1,v2 : Natural_Vectors.Vector ) return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then return false;  -- the dimensions must correspond !
     else declare
            p : Permutation(v1'first..v1'last);
          begin
            for k in p'range loop
              p(k) := 0;
              for l in v2'range loop
                if v2(l) = v1(k)
                 then p(k) := l;
                      for j in 1..(k-1) loop
                        if p(j) = l
                         then p(k) := 0;
                        end if;
                      end loop;
                end if;
                exit when p(k) /= 0;
              end loop;
              if p(k) = 0
               then return false;
              end if;
            end loop;
          end;
          return true;
    end if;
  end Permutable;

  function Permutable ( v1,v2 : Integer_Vectors.Vector ) return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then return false;  -- the dimensions must correspond !
     else declare
            p : Permutation(v1'first..v1'last);
          begin
            for k in p'range loop
              p(k) := 0;
              for l in v2'range loop
                if v2(l) = v1(k)
                 then p(k) := l;
                      for j in 1..(k-1) loop
                        if p(j) = l
                         then p(k) := 0;
                        end if;
                      end loop;
                end if;
                exit when p(k) /= 0;
              end loop;
              if p(k) = 0
               then return false;
              end if;
            end loop;
          end;
          return true;
    end if;
  end Permutable;

  function Permutable ( v1,v2 :   Float_Vectors.Vector ) return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then return false;  -- the dimensions must correspond !
     else declare
            p : Permutation(v1'first..v1'last);
          begin
            for k in p'range loop
              p(k) := 0;
              for l in v2'range loop
                if v2(l) = v1(k)
                 then p(k) := l;
                      for j in 1..(k-1) loop
                        if p(j) = l
                         then p(k) := 0;
                        end if;
                      end loop;
                end if;
                exit when p(k) /= 0;
              end loop;
              if p(k) = 0
               then return false;
              end if;
            end loop;
          end;
          return true;
    end if;
  end Permutable;

  function Permutable ( v1,v2 : Complex_Vectors.Vector ) return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then return false;  -- the dimensions must correspond !
     else declare
            p : Permutation(v1'first..v1'last);
          begin
            for k in p'range loop
              p(k) := 0;
              for l in v2'range loop
                if v2(l) = v1(k)
                 then p(k) := l;
                      for j in 1..(k-1) loop
                        if p(j) = l
                         then p(k) := 0;
                        end if;
                      end loop;
                end if;
                exit when p(k) /= 0;
              end loop;
              if p(k) = 0
               then return false;
              end if;
            end loop;
          end;
          return true;
    end if;
  end Permutable;

  function Permutable ( v1,v2 :   Float_Vectors.Vector; tol : double_float )
                      return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then return false;  -- the dimensions must correspond !
     else declare
            p : Permutation(v1'first..v1'last);
          begin
            for k in p'range loop
              p(k) := 0;
              for l in v2'range loop
                if ABS(v2(l) - v1(k)) <= tol
                 then p(k) := l;
                      for j in 1..(k-1) loop
                        if p(j) = l
                         then p(k) := 0;
                        end if;
                      end loop;
                end if;
                exit when p(k) /= 0;
              end loop;
              if p(k) = 0
               then return false;
              end if;
            end loop;
          end;
          return true;
    end if;
  end Permutable;

  function Permutable ( v1,v2 : Complex_Vectors.Vector; tol : double_float )
                      return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then return false;  -- the dimensions must correspond !
     else declare
            p : Permutation(v1'first..v1'last);
          begin
            for k in p'range loop
              p(k) := 0;
              for l in v2'range loop
                if (ABS(REAL_PART(v2(l)) - REAL_PART(v1(k))) <= tol)
                  and then (ABS(IMAG_PART(v2(l)) - IMAG_PART(v1(k))) <= tol)
                 then p(k) := l;
                      for j in 1..(k-1) loop
                        if p(j) = l
                         then p(k) := 0;
                        end if;
                      end loop;
                end if;
                exit when p(k) /= 0;
              end loop;
              if p(k) = 0
               then return false;
              end if;
            end loop;
          end;
          return true;
    end if;
  end Permutable;

  function Sign_Permutable ( v1,v2 : Natural_Vectors.Vector ) return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then return false;  -- the dimensions must correspond !
     else declare
            p : Permutation(v1'first..v1'last);
          begin
            for k in p'range loop
              p(k) := 0;
              for l in v2'range loop
                if v2(l) = v1(k) or else v2(l) = -v1(k)
                 then p(k) := l;
                      for j in 1..(k-1) loop
                        if p(j) = l
                         then p(k) := 0;
                        end if;
                      end loop;
                end if;
                exit when p(k) /= 0;
              end loop;
              if p(k) = 0
               then return false;
              end if;
            end loop;
          end;
          return true;
    end if;
  end Sign_Permutable;

  function Sign_Permutable ( v1,v2 : Integer_Vectors.Vector ) return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then return false;  -- the dimensions must correspond !
     else declare
            p : Permutation(v1'first..v1'last);
          begin
            for k in p'range loop
              p(k) := 0;
              for l in v2'range loop
                if v2(l) = v1(k) or else v2(l) = -v1(k)
                 then p(k) := l;
                      for j in 1..(k-1) loop
                        if p(j) = l
                         then p(k) := 0;
                        end if;
                      end loop;
                end if;
                exit when p(k) /= 0;
              end loop;
              if p(k) = 0
               then return false;
              end if;
            end loop;
          end;
          return true;
    end if;
  end Sign_Permutable;

  function Sign_Permutable ( v1,v2 :   Float_Vectors.Vector ) return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then return false;  -- the dimensions must correspond !
     else declare
            p : Permutation(v1'first..v1'last);
          begin
            for k in p'range loop
              p(k) := 0;
              for l in v2'range loop
                if v2(l) = v1(k) or else v2(l) = -v1(k)
                 then p(k) := l;
                      for j in 1..(k-1) loop
                        if p(j) = l
                         then p(k) := 0;
                        end if;
                      end loop;
                end if;
                exit when p(k) /= 0;
              end loop;
              if p(k) = 0
               then return false;
              end if;
            end loop;
          end;
          return true;
    end if;
  end Sign_Permutable;

  function Sign_Permutable ( v1,v2 : Complex_Vectors.Vector ) return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then return false;  -- the dimensions must correspond !
     else declare
            p : Permutation(v1'first..v1'last);
          begin
            for k in p'range loop
              p(k) := 0;
              for l in v2'range loop
                if v2(l) = v1(k) or else v2(l) = -v1(k)
                 then p(k) := l;
                      for j in 1..(k-1) loop
                        if p(j) = l
                         then p(k) := 0;
                        end if;
                      end loop;
                end if;
                exit when p(k) /= 0;
              end loop;
              if p(k) = 0
               then return false;
              end if;
            end loop;
          end;
          return true;
    end if;
  end Sign_Permutable;

  function Sign_Permutable ( v1,v2 :   Float_Vectors.Vector;
                             tol : double_float ) return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then return false;  -- the dimensions must correspond !
     else declare
            p : Permutation(v1'first..v1'last);
          begin
            for k in p'range loop
              p(k) := 0;
              for l in v2'range loop
                if (ABS(v2(l) - v1(k)) <= tol)
                   or else  (ABS(v2(l) + v1(k)) <= tol)
                 then p(k) := l;
                      for j in 1..(k-1) loop
                        if p(j) = l
                         then p(k) := 0;
                        end if;
                      end loop;
                end if;
                exit when p(k) /= 0;
              end loop;
              if p(k) = 0
               then return false;
              end if;
            end loop;
          end;
          return true;
    end if;
  end Sign_Permutable;

  function Sign_Permutable ( v1,v2 : Complex_Vectors.Vector;
                             tol : double_float ) return boolean is
  begin
    if v1'first /= v2'first or else v1'last /= v2'last
     then return false;  -- the dimensions must correspond !
     else declare
            p : Permutation(v1'first..v1'last);
          begin
            for k in p'range loop
              p(k) := 0;
              for l in v2'range loop
                if ((ABS(REAL_PART(v2(l)) - REAL_PART(v1(k))) <= tol)
                    and then (ABS(IMAG_PART(v2(l)) - IMAG_PART(v1(k))) <= tol))
                  or else ((ABS(REAL_PART(v2(l)) + REAL_PART(v1(k))) <= tol)
                    and then (ABS(IMAG_PART(v2(l)) + IMAG_PART(v1(k))) <= tol))
                 then p(k) := l;
                      for j in 1..(k-1) loop
                        if p(j) = l
                         then p(k) := 0;
                        end if;
                      end loop;
                end if;
                exit when p(k) /= 0;
              end loop;
              if p(k) = 0
               then return false;
              end if;
            end loop;
          end;
          return true;
    end if;
  end Sign_Permutable;

  function "*" ( p : Permutation; t : Complex_Multivariate_Polynomials.Term )
               return Complex_Multivariate_Polynomials.Term is

    res : Complex_Multivariate_Polynomials.Term;

  begin
    res.cf := t.cf;
    res.dg := new Natural_Vectors.Vector(t.dg'range);
    for i in p'range loop
      if p(i) >= 0
       then res.dg(i) := t.dg(p(i));
       else res.dg(i) := t.dg(-p(i));
            res.cf := -res.cf;
      end if;
    end loop;
    return res;
  end "*";

  function "*" ( p : Permutation; s : Complex_Multivariate_Polynomials.Poly )
               return Complex_Multivariate_Polynomials.Poly is

    use Complex_Multivariate_Polynomials;
    res : Poly := Null_Poly;

    procedure Permute_Term ( t : in Term; continue : out boolean ) is
      tt : Term := p*t;
    begin
      Plus_Term(res,tt);
      Clear(tt);
      continue := true;
    end Permute_Term;
    procedure Permute_Terms is new Visiting_Iterator(Permute_Term);

  begin
    Permute_Terms(s);
    return res;
  end "*";

  function "*" ( p : Permutation;
                 t : Complex_Multivariate_Laurent_Polynomials.Term )
               return Complex_Multivariate_Laurent_Polynomials.Term is

    res : Complex_Multivariate_Laurent_Polynomials.Term;

  begin
    res.cf := t.cf;
    res.dg := new Integer_Vectors.Vector(t.dg'range);
    for i in p'range loop
      if p(i) >= 0
       then res.dg(i) := t.dg(p(i));
       else res.dg(i) := t.dg(-p(i));
            res.cf := -res.cf;
      end if;
    end loop;
    return res;
  end "*";

  function "*" ( p : Permutation;
                 s : Complex_Multivariate_Laurent_Polynomials.Poly )
               return Complex_Multivariate_Laurent_Polynomials.Poly is

    use Complex_Multivariate_Laurent_Polynomials;
    res : Poly := Null_Poly;

    procedure Permute_Term ( t : in Term; continue : out boolean ) is
      tt : Term := p*t;
    begin
      Plus_Term(res,tt);
      Clear(tt);
      continue := true;
    end Permute_Term;
    procedure Permute_Terms is new Visiting_Iterator(Permute_Term);

  begin
    Permute_Terms(s);
    return res;
  end "*";

  function "*" ( s : Poly_Sys; p : Permutation ) return Poly_Sys is

    res : Poly_Sys(s'range);

  begin
    for k in res'range loop
      res(k) := p*s(k);
    end loop;
    return res;
  end "*";

  function "*" ( s : Laur_Sys; p : Permutation ) return Laur_Sys is

    res : Laur_Sys(s'range);

  begin
    for k in res'range loop
      res(k) := p*s(k);
    end loop;
    return res;
  end "*";

  function "*" ( p : Permutation; s : Poly_Sys ) return Poly_Sys is

    r : Poly_Sys(s'range);
    use Complex_Multivariate_Polynomials;

  begin
    for i in p'range loop
      if p(i) >= 0
       then Copy(s(p(i)),r(i));
       else r(i) := -s(-p(i));
      end if;
    end loop;
    return r;
  end "*";

  function "*" ( p : Permutation; s : Laur_Sys ) return Laur_Sys is

    r : Laur_Sys(s'range);
    use Complex_Multivariate_Laurent_Polynomials;

  begin
    for i in p'range loop
      if p(i) >= 0
       then Copy(s(p(i)),r(i));
       else r(i) := -s(-p(i));
      end if;
    end loop;
    return r;
  end "*";

end Permute_Operations;
