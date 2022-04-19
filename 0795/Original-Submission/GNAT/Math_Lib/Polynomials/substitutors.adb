with Floating_Point_Numbers;     use Floating_Point_Numbers;
with Natural_Vectors;

package body Substitutors is

  function Substitute ( k : integer; c : double_complex; t : Term )
                      return Term is

    res : Term;

  begin
    res.cf := t.cf;
    for l in 1..t.dg(k) loop
      res.cf := res.cf*c;
    end loop;
    res.dg := new Natural_Vectors.Vector(t.dg'first..t.dg'last-1);
    for l in t.dg'range loop
      if l < k
       then res.dg(l) := t.dg(l);
       elsif l > k
           then res.dg(l-1) := t.dg(l);
      end if;
    end loop;
    return res;
  end Substitute;

  function Substitute ( k : integer; c : double_complex; p : Poly )
                      return Poly is
    res : Poly;

    procedure Substitute_Term ( t : in Term; cont : out boolean ) is

      st : Term := Substitute(k,c,t);

    begin
      Plus_Term(res,st);
      Clear(st);
      cont := true;
    end Substitute_Term;
    procedure Substitute_Terms is new Visiting_Iterator (Substitute_Term);

  begin
    Substitute_Terms(p);
    return res;
  end Substitute; 

  function Substitute ( k : integer; c : double_complex; p : Poly_Sys )
                      return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for l in res'range loop
      res(l) := Substitute(k,c,p(l));
    end loop;
    return res;
  end Substitute;

  procedure Substitute ( k : in integer; c : in double_complex; 
                         t : in out Term ) is

    tmp : Natural_Vectors.Vector(t.dg'first..t.dg'last-1);

  begin
    for l in 1..t.dg(k) loop
      t.cf := t.cf*c;
    end loop;
    for l in t.dg'range loop
      if l < k
       then tmp(l) := t.dg(l);
       elsif l > k
           then tmp(l-1) := t.dg(l);
      end if;
    end loop;
    Clear(t);
    t.dg := new Natural_Vectors.Vector'(tmp);
  end Substitute;

--  procedure Shuffle ( p : in out Poly ) is
--
--   -- DESCRIPTION :
--   --   Changes the position of the terms in p back to the normal order.
--
--    res : Poly := Null_Poly;
--
--    procedure Shuffle_Term ( t : in Term; cont : out boolean ) is
--    begin
--      Plus_Term(res,t);
--      cont := true;
--    end Shuffle_Term;
--    procedure Shuffle_Terms is new Visiting_Iterator(Shuffle_Term);
--
--  begin
--    Shuffle_Terms(p);
--    Clear(p); Copy(res,p); Clear(res);
--  end Shuffle;

  procedure Substitute ( k : in integer; c : in double_complex;
                         p : in out Poly ) is

   -- NOTE :
   --   An obvious thing to do would be to visit and change all terms,
   --   and leaving the term order unchanged.
   --   However, the terms in the resulting polynomial must be ordered
   --   according to the term order with which the package has been
   --   instantiated.  Therefore the procedure Shuffle is called.

    procedure Substitute_Term ( t : in out Term; cont : out boolean ) is
    begin
      Substitute(k,c,t);
      cont := true;
    end Substitute_Term;
    procedure Substitute_Terms is new Changing_Iterator ( Substitute_Term );

  begin
    Substitute_Terms(p);
--    Shuffle(p);
  end Substitute;

  procedure Substitute ( k : in integer; c : in double_complex;
                         p : in out Poly_Sys ) is
  begin
    for l in p'range loop
      Substitute(k,c,p(l));
    end loop;
  end Substitute;

  procedure Purify ( p : in out Poly; tol : in double_float ) is

  -- DESCRIPTION :
  --   All terms of which the coefficient are in modulus smaller
  --   than tol are deleted.

    procedure Purify_Term (t : in out Term; continue : out boolean) is
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

  function Substitute_Factor ( k : integer; h : Vector ) return Poly is
   
  -- DESCRIPTION :
  --   returns the polynomial which will replace the kth unknown.

    res : Poly;
    rt : Term;

  begin
    rt.dg := new Natural_Vectors.Vector'((h'first+1)..(h'last-1) => 0);
    rt.cf := -h(0)/h(k);
    res := Create(rt);
    for i in rt.dg'range loop
      rt.dg(i) := 1;
      if i < k
       then rt.cf := -h(i)/h(k);
       else rt.cf := -h(i+1)/h(k);
      end if;
      if modulus(rt.cf) > 10.0**(-10)
       then Plus_Term(res,rt);
      end if;
      rt.dg(i) := 0;
    end loop;
    Clear(rt);
    return res;
  end Substitute_Factor;

  function Substitute ( k : integer; h : Vector; p : Poly ) return Poly is

    res,sub : Poly;

  begin
    sub := Substitute_Factor(k,h);
    res := Substitute(k,sub,p);
    Clear(sub);
    return res;
  end Substitute;

  function Substitute ( k : integer; s,p : Poly ) return Poly is

    res : Poly := Null_Poly;

    procedure Substitute_Term ( t : in Term; continue : out boolean ) is
      rt : Term;
      fac : Poly;
    begin
      rt.cf := t.cf;
      rt.dg := new Natural_Vectors.Vector(t.dg'first..(t.dg'last-1));
      for i in rt.dg'range loop
	if i < k
	 then rt.dg(i) := t.dg(i);
	 else rt.dg(i) := t.dg(i+1);
        end if;
      end loop;
      if t.dg(k) = 0
       then Plus_Term(res,rt);
       else fac := Create(rt);
	    for i in 1..t.dg(k) loop
	      Mult_Poly(fac,s);
            end loop;
	    Purify(fac,10.0**(-10));
	    Plus_Poly(res,fac);
            Clear(fac);
      end if;
      Clear(rt);
      continue := true;
    end Substitute_Term;
    procedure Substitute_Terms is new Visiting_Iterator (Substitute_Term);

  begin
    Substitute_Terms(p);
    return res;
  end Substitute;

  procedure Substitute ( k : in integer; h : in Vector; p : in out Poly ) is

    res : Poly;

  begin
    res := Substitute(k,h,p);
    Clear(p); Copy(res,p); Clear(res);
  end Substitute;

  procedure Substitute ( k : in integer; s : in Poly; p : in out Poly ) is

    res : Poly;

  begin
    res := Substitute(k,s,p);
    Clear(p); Copy(res,p); Clear(res);
  end Substitute;

  function Substitute ( k : integer; h : Vector; p : Poly_Sys )
                      return Poly_Sys is

    res : Poly_Sys(p'range);
    s : Poly := Substitute_Factor(k,h);

  begin
    res := Substitute(k,s,p);
    Clear(s);
    return res;
  end Substitute;

  procedure Substitute ( k : in integer; h : in Vector;
                         p : in out Poly_Sys ) is

    s : Poly := Substitute_Factor(k,h);

  begin
    Substitute(k,s,p);
    Clear(s);
  end Substitute;

  function Substitute ( k : integer; s : Poly; p : Poly_Sys ) return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Substitute(k,s,p(i));
    end loop;
    return res;
  end Substitute;

  procedure Substitute ( k : in integer; s : in Poly; p : in out Poly_Sys ) is
  begin
    for i in p'range loop
      Substitute(k,s,p(i));
    end loop;
  end Substitute;

end Substitutors;
