with unchecked_deallocation;
with Floating_Point_Numbers;              use Floating_Point_Numbers;
with Natural_Vectors,Jacobi_Matrices;     use Jacobi_Matrices;
with Complex_Vectors,Complex_Matrices;    use Complex_Vectors,Complex_Matrices;
with Complex_Multivariate_Polynomials;    use Complex_Multivariate_Polynomials;

package body Homotopy is

-- There is an amount of internal data in order to
-- perform the numeric routines as fast as possible.

  type homtype is (nat,art);  -- natural or artificial parameter homotopy

  type homdata ( ht : homtype; n,n1 : natural ) is record

    p : Poly_Sys(1..n);
    pe : Eval_Poly_Sys(1..n);
    dh : Jacobi(1..n,1..n1);
    dhe : Eval_Jacobi(1..n,1..n1);

    case ht is
      when nat =>
        i : integer;
      when art =>
        q,h : Poly_Sys(1..n);
        qe,he : Eval_Poly_Sys(1..n);
        dpe,dqe : Eval_Jacobi(1..n,1..n);
        k : positive;
        a : double_complex;
      end case;

  end record;
 
  type homtp is access homdata;
  hom : homtp;

-- CONSTRUCTORS :

  function Create ( p,q : in Poly_Sys; k : in positive; a : in double_complex )
                  return Poly_Sys is

    h : Poly_Sys(p'range);
    tempp,tempq,q_fac,p_fac : Poly;
    n : natural := p'length;

    function Compute_q_fac ( n,k : in natural; a : in double_complex )
                           return Poly is

    -- DESCRIPTION : returns a*(1-t)^k, where t = x_n+1.

      res,temp : Poly;
      t : Term;

    begin
      t.cf := a;
      t.dg := new Natural_Vectors.Vector'(1..n+1 => 0);
      res := Create(t);
      t.cf := CMPLX(1.0);
      temp := Create(t);
      t.dg(n+1) := 1;
      Min_Term(temp,t);
      Natural_Vectors.Clear(Natural_Vectors.Link_to_Vector(t.dg));
      for i in 1..k loop
        Mult_Poly(res,temp);
      end loop;
      Clear(temp);
      return res;
    end Compute_q_fac;

    function Compute_p_fac ( n,k : in natural; a : in double_complex ) 
                           return Poly is

    -- DESCRIPTION : returns t^k, where t = x_n+1.

      res : Poly;
      t : Term;

    begin
      t.cf := CMPLX(1.0);
      t.dg := new Natural_Vectors.Vector'(1..n+1 => 0);
      t.dg(n+1) := k;
      res := Create(t);
      Natural_Vectors.Clear(Natural_Vectors.Link_to_Vector(t.dg));
      return res;
    end Compute_p_fac;

    function Plus_Unknown ( p : in Poly ) return Poly is

    -- DESCRIPTION :
    --   The returning polynomial has an additional unknown.

      res : Poly;

      procedure Plus_Unknown_In_Term (t : in out Term; c : out boolean) is

        temp : Natural_Vectors.Vector(1..(n+1));

      begin
        temp(1..n) := t.dg.all;
        temp(n+1) := 0;
        Natural_Vectors.Clear(Natural_Vectors.Link_to_Vector(t.dg));
        t.dg := new Natural_Vectors.Vector'(temp);
        c := true;
      end Plus_Unknown_In_Term;
      procedure Plus_Unknown_In_Terms is 
         new Changing_Iterator (process => Plus_Unknown_In_Term);

    begin
      Copy(p,res);
      Plus_Unknown_In_Terms(res);
      return res;
    end Plus_Unknown;

  begin
    q_fac := Compute_q_fac(n,k,a);
    p_fac := Compute_p_fac(n,k,a);
    for i in h'range loop
      tempq := Plus_Unknown(q(i));
      tempp := Plus_Unknown(p(i));
      Mult_Poly(tempq,q_fac);
      Mult_Poly(tempp,p_fac);
      h(i) := tempq + tempp;
      Clear(tempq); Clear(tempp);
    end loop;
    Clear(q_fac); Clear(p_fac);
    return h;
  end Create;

  procedure Create ( p,q : in Poly_Sys; k : in positive;
                     a : in double_complex ) is

    n : constant natural := p'last-p'first+1;
    dp, dq : Jacobi(1..n,1..n);
    ho : homdata(art,n,n+1);

  begin
    Copy(p,ho.p); Copy(q,ho.q);
    ho.h := Create(p,q,k,a);
    ho.pe := Create(ho.p);
    ho.qe := Create(ho.q);
    ho.he := Create(ho.h);
    dp := Create(ho.p);
    dq := Create(ho.q);
    ho.dh := Create(ho.h);
    ho.dpe := Create(dp);
    ho.dqe := Create(dq);
    ho.dhe := Create(ho.dh);
    Clear(dp); Clear(dq);
    ho.k := k;
    ho.a := a;
    hom := new homdata'(ho);
  end Create; 

  procedure Create ( p : in Poly_Sys; k : in integer ) is

    n : constant natural := p'last-p'first+1; 
    ho : homdata(nat,n,n+1);

  begin
    Copy(p,ho.p);
    ho.pe := Create(ho.p);
    ho.dh := Create(ho.p);
    ho.dhe := Create(ho.dh);
    ho.i := k;
    hom := new homdata'(ho);
  end Create;

-- SELECTOR :

  function Homotopy_System return Poly_Sys is

    ho : homdata renames hom.all;

  begin
    case ho.ht is
      when nat => return ho.p;
      when art => return ho.h;
    end case;
  end Homotopy_System;

-- SYMBOLIC ROUTINES :

  function Eval ( t : double_complex ) return Poly_Sys is

    ho : homdata renames hom.all;

  begin
    case ho.ht is
      when nat =>  -- t = x(ho.i)
        return Eval(ho.p,t,ho.i);
      when art =>  -- compute : a * ((1 - t)^k) * q + (t^k) * p
        declare
          p_factor,q_factor : double_complex;
          res,tmp : Poly_Sys(1..ho.n);
        begin
          if modulus(t) = 0.0
           then return ho.a * ho.q;
           elsif abs(REAL_PART(t) - 1.0 ) + 1.0 = 1.0 
                and then abs(IMAG_PART(t)) + 1.0 = 1.0
               then copy(ho.p,res);
                    return res;
               else q_factor := ho.a;
                    p_factor := CMPLX(1.0);
                    for i in 1..ho.k loop
                      q_factor := (CMPLX(1.0)-t) * q_factor;
                      p_factor := t * p_factor;
                    end loop;
                    res := p_factor * ho.p; 
                    tmp := q_factor * ho.q;
                    Plus_Vector(res,tmp);
                    Clear(tmp);
                    return res;
          end if;
        end;
    end case;
  end Eval;

  function Diff ( t : double_complex ) return Poly_Sys is

    ho : homdata renames hom.all;

  begin
    case ho.ht is
      when nat =>  -- t = x(ho.i)
        return Diff(ho.p,ho.i);
      when art =>  -- compute  - a * k * (1 - t)^(k-1) * q + k * t^(k-1) *p
        declare
          q_factor,p_factor : double_complex;
          tmp : Poly_Sys(1..ho.n);
          res : Poly_Sys(1..ho.n);
        begin
          q_factor := (-ho.a) * CMPLX(double_float(ho.k));
          p_factor := CMPLX(double_float(ho.k));
          if modulus(t) = 0.0
           then if ho.k = 1
                 then res := p_factor * ho.p;
                      tmp := q_factor * ho.q;
                      Plus_Vector(res,tmp);
                      Clear(tmp);
                      return res;
                 else return q_factor * ho.q;
                end if;
           elsif abs( REAL_PART(t) - 1.0 ) + 1.0 = 1.0
                and then abs(IMAG_PART(t)) + 1.0 = 1.0
               then return CMPLX(double_float(ho.k)) * ho.p;
               else for i in 1..(ho.k-1) loop
                      q_factor := (CMPLX(1.0)-t) * q_factor;
                      p_factor := t * p_factor;
                    end loop;
                    res := p_factor * ho.p;
                    tmp := q_factor * ho.q;
                    Plus_Vector(res,tmp);
                    Clear(tmp);
                    return res;
          end if;
        end;
    end case;
  end Diff;

-- NUMERIC ROUTINES :

  function Eval ( x : Vector; t : double_complex ) return Vector is

    ho : homdata renames hom.all;

  begin
    case ho.ht is
      when nat =>
        declare
          y : Vector(x'first..x'last+1);
        begin
          y(1..ho.i-1) := x(1..ho.i-1);
          y(ho.i) := t;
          y(ho.i+1..y'last) := x(ho.i..x'last);
          return Eval(ho.pe,y);
        end;
      when art =>
        if modulus(t) + 1.0 = 1.0
         then return ho.a * Eval(ho.qe,x);
         elsif abs( REAL_PART(t) - 1.0 ) + 1.0 = 1.0
              and then abs(IMAG_PART(t)) + 1.0 = 1.0
             then return Eval(ho.pe,x);
             else declare
                    y : Vector(x'first..x'last+1);
                  begin
                    y(x'range) := x;
                    y(y'last) := t;
                    return Eval(ho.he,y);
                  end;
       end if;
    end case;
  end Eval;

  function Diff ( x : Vector; t : double_complex ) return Vector is

    n : natural := x'length;

  begin
    case hom.ht is
      when nat => return Diff(x,t,hom.i);
      when art => return Diff(x,t,n+1);
    end case;
  end Diff;
 
  function Diff ( x : Vector; t : double_complex ) return matrix is

    ho : homdata renames hom.all;
    n : natural renames ho.n;

  begin
    case ho.ht is
      when nat =>
        declare
          m : matrix(1..n,1..n);
          y : Vector(1..n+1);
        begin
          y(1..ho.i-1) := x(1..ho.i-1);
          y(ho.i) := t;
          y(ho.i+1..n+1) := x(ho.i..n);
          for i in 1..n loop
            for j in 1..n loop
              m(i,j) := Eval(ho.dhe(i,j),y);
            end loop;
          end loop;
          return m;
        end;
      when art =>
        if modulus(t) + 1.0 = 1.0
         then declare
                m : Matrix(1..n,1..n) := Eval(ho.dqe,x);
              begin
                for i in 1..n loop
                  for j in 1..n loop
                    m(i,j) := m(i,j) * ho.a;
                  end loop;
                end loop;
                return m;
              end;
         elsif abs( REAL_PART(t) - 1.0 ) + 1.0 = 1.0
              and then abs(IMAG_PART(t)) + 1.0 = 1.0
             then return Eval(ho.dpe,x);
             else declare
                    m : matrix(1..n,1..n);
                    y : Vector(1..n+1);
                  begin
                    y(1..n) := x;
                    y(n+1) := t;
                    for i in 1..n loop
                      for j in 1..n loop
                        m(i,j) := Eval(ho.dhe(i,j),y);
                      end loop;
                    end loop;
                    return m;
                  end;
        end if;
    end case;
  end Diff;

  function Diff ( x : Vector; t : double_complex; k : natural )
                return Vector is

    ho : homdata renames hom.all;
    n : natural renames ho.n;
    y : Vector(1..n+1);
    res : Vector(1..n);

  begin
    case ho.ht is
      when nat => y(1..ho.k-1) := x(1..ho.i-1);
                  y(ho.i) := t;
                  y(ho.i+1..n+1) := x(ho.i..n);
      when art => y(1..n) := x;
                  y(n+1) := t;
    end case;
    for i in 1..n loop
      res(i) := Eval(ho.dhe(i,k),y);
    end loop;
    return res;
  end Diff;

  function Diff ( x : Vector; t : double_complex; k : natural )
                return matrix is

    ho : homdata renames hom.all;
    n : natural renames ho.n;
    y : Vector(1..n+1);
    res : matrix(1..n,1..n);

  begin
    case ho.ht is
      when nat => y(1..ho.i-1) := x(1..ho.i-1);
                  y(ho.i) := t;
                  y(ho.i+1..n+1) := x(ho.i..n); 
      when art => y(1..n) := x;
                  y(n+1) := t;
    end case;
    for j in 1..(k-1) loop
      for i in 1..n loop
        res(i,j) := Eval(ho.dhe(i,j),y);
      end loop;
    end loop;
    for j in (k+1)..(n+1) loop
      for i in 1..n loop
        res(i,j-1) := Eval(ho.dhe(i,j),y);
      end loop;
    end loop;
    return res;
  end Diff;

-- DESTRUCTOR :

  procedure free is new unchecked_deallocation (homdata,homtp);

  procedure Clear is
  begin
    if hom /= null
     then declare
            ho : homdata renames hom.all;
          begin
            Clear(ho.p);  Clear(ho.pe);
            Clear(ho.dh); Clear(ho.dhe);
            case ho.ht is
              when nat => null;
              when art =>
                Clear(ho.q);   Clear(ho.qe);
                Clear(ho.h);   Clear(ho.he);
                Clear(ho.dpe); Clear(ho.dqe);
            end case;
          end;
          free(hom);
    end if;
  end Clear;

end Homotopy;
