with Floating_Point_Numbers;    use Floating_Point_Numbers;
with Complex_Numbers;           use Complex_Numbers;
with Natural_Vectors;           use Natural_Vectors;

package body Reduction_of_Polynomials is

-- AUXILIARIES :

  function Leading_Term ( p : Poly ) return Term is

  -- DESCRIPTION :
  --   Returns the leading term of the polynomial p.

    res : Term;

    procedure First_Term ( t : in Term; continue : out boolean ) is
    begin
      Natural_Vectors.Copy(Link_to_Vector(t.dg),Link_to_Vector(res.dg));
      res.cf := t.cf;
      continue := false;
    end First_Term;
    procedure Get_First_Term is new Visiting_Iterator(First_Term);

  begin
    Get_First_Term(p);
    return res;
  end Leading_Term;

  function LEQ ( d1,d2 : Degrees ) return boolean is

  -- DESCRIPTION :
  --   Returns true if all degrees of d1 are lower than
  --   or equal to the degrees of d2.

  begin
    for i in d1'range loop
      if d1(i) > d2(i)
       then return false;
      end if;
    end loop;
    return true;
  end LEQ;

  function Search (p : Poly; pt : Term) return Term is

  -- DESCRIPTION :
  --   Returns the first term of p for which the following holds :
  --   LEQ(result,pt).

    result : Term;

    procedure Search_Term (t : in Term; continue : out boolean) is
    begin
      if LEQ(t.dg,pt.dg)
       then result.cf := t.cf;
            result.dg := new Natural_Vectors.Vector(t.dg'range);
            for i in t.dg'range loop
              result.dg(i) := t.dg(i);
            end loop;
            continue := false;
       else continue := true;
      end if;
    end Search_Term;
    procedure Search_Terms is new Visiting_Iterator(Search_Term);

  begin
    result.cf := CMPLX(0.0);
    Search_Terms(p);
    return result;
  end Search;

  procedure Purify ( p : in out Poly; tol : in double_float ) is

  -- DESCRIPTION :
  --   All terms of which the coefficient are in modulus smaller
  --   than tol are deleted.

    procedure Purify_Term ( t : in out Term; continue : out boolean ) is
    begin
      if modulus(t.cf) < tol
       then t.cf := CMPLX(0.0);
      end if;
      continue := true;
    end Purify_Term;
    procedure Purify_Terms is new Changing_Iterator(Purify_Term);

  begin
    Purify_Terms(p);
    if Number_Of_Unknowns(p) = 0
     then Clear(p);
          p := Null_Poly;
    end if;
  end Purify;

-- TARGET ROUTINES :

  function Spoly ( p,q : poly ) return Poly is

    S,pp,qq : Poly;
    tfp,tfq,facq,facp : Term;
    tol : constant double_float := 10.0**(-13);

  begin
    if (p = Null_Poly) or else (Number_Of_Unknowns(p) = 0)
      or else (q = Null_Poly) or else (Number_Of_Unknowns(q) = 0)
     then return Null_Poly;
    end if;

    tfp := Leading_Term(p);
    tfq := Leading_Term(q);

    facp.dg := new Natural_Vectors.Vector'(tfp.dg'range => 0);
    facq.dg := new Natural_Vectors.Vector'(tfq.dg'range => 0);
    for i in tfp.dg'range loop
      if tfp.dg(i) > tfq.dg(i)
       then facq.dg(i) := tfp.dg(i) - tfq.dg(i);
       elsif tfp.dg(i) < tfq.dg(i)
           then facp.dg(i) := tfq.dg(i) - tfp.dg(i);
      end if;
    end loop;
    if modulus(tfp.cf) > modulus(tfq.cf)
     then facp.cf := CMPLX(1.0);
          facq.cf := - tfp.cf / tfq.cf;
     else facp.cf := tfq.cf / tfp.cf;
          facq.cf := CMPLX(-1.0);
    end if;
    pp := facp * p;
    qq := facq * q;
    S := pp + qq;
    Clear(pp); Clear(qq);
    Natural_Vectors.Clear(Link_to_Vector(facp.dg));
    Natural_Vectors.Clear(Link_to_Vector(tfp.dg));
    Natural_Vectors.Clear(Link_to_Vector(facq.dg));
    Natural_Vectors.Clear(Link_to_Vector(tfq.dg));
    Purify(S,tol);
    return S;
  end Spoly;

  function Rpoly ( p,q : Poly ) return Poly is

    tol : constant double_float := 10.0**(-13);

  begin
    if (p = Null_Poly) or else (Number_Of_Unknowns(p) = 0)
     then return Null_Poly;
     elsif (q = Null_Poly) or else (Number_Of_Unknowns(q) = 0)
         then null;
         else declare
                ltp,tq : Term;
              begin
                ltp := Leading_Term(p);
                tq := Search(q,ltp);
                if tq.cf /= CMPLX(0.0)
                 then
                   declare
                     R,qq : Poly;
                     fac : Term;
                   begin
                     fac.dg := new Natural_Vectors.Vector'(ltp.dg'range => 0);
                     for i in ltp.dg'range loop
                       if ltp.dg(i) > tq.dg(i)
                        then fac.dg(i) := ltp.dg(i) - tq.dg(i);
                       end if;
                     end loop;
                     fac.cf := -ltp.cf / tq.cf;
                     qq := fac * q;
                     R := p + qq;
                     Clear(qq);
                     Natural_Vectors.Clear(Link_to_Vector(fac.dg));
                     Purify(R,tol);
                     return R;
                   end;
                end if;
                Natural_Vectors.Clear(Link_to_Vector(ltp.dg));
                Natural_Vectors.Clear(Link_to_Vector(tq.dg));
              end;
    end if;
    declare
      temp : Poly;
    begin
      Copy(p,temp);
      return temp;
    end;
  end Rpoly;

end Reduction_of_Polynomials;
