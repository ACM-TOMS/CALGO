with Complex_Numbers;               use Complex_Numbers;
with Integer_Vectors_Utilities;     use Integer_Vectors_Utilities;

package body Transforming_Laurent_Systems is

  function Initial_Link_to_Vector ( p : Poly ) return Link_to_Vector is

  -- DESCRIPTION :
  --   Returns the initial degrees of the polynomial p.

    init : Link_to_Vector;

    procedure Init_Term ( t : in Term; cont : out boolean ) is
    begin
      init := new Integer_Vectors.Vector'(t.dg.all);
      cont := false;
    end Init_Term;
    procedure Initial_Term is new Visiting_Iterator (Init_Term);

  begin
    Initial_Term(p);
    return init;
  end Initial_Link_to_Vector;

  procedure Shift ( p : in out Poly ) is

    init : Link_to_Vector := Initial_Link_to_Vector(p);

    procedure Shift_Term ( t : in out Term; cont : out boolean ) is
    begin
      Min_Vector(Link_to_Vector(t.dg),init);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is new Changing_Iterator (Shift_Term);

  begin
    if p /= Null_Poly
     then Shift_Terms(p);
    end if;
    Clear(init);
  end Shift;

  function Shift ( p : Poly ) return Poly is

    res : Poly := Null_Poly;
    init : Link_to_Vector := Initial_Link_to_Vector(p);

    procedure Shift_Term ( t : in Term; cont : out boolean ) is
      rt : Term;
    begin
      rt.cf := t.cf;
      rt.dg := t.dg - Degrees(init);
      Plus_Term(res,rt);
      Clear(rt);
      cont := true;
    end Shift_Term;
    procedure Shift_Terms is new Visiting_Iterator (Shift_Term);

  begin
    if p /= Null_Poly
     then Shift_Terms(p);
    end if;
    Clear(init);
    return res;
  end Shift;

  procedure Shift ( l : in out Laur_Sys ) is
  begin
    for k in l'range loop
      Shift(l(k));
    end loop;
  end Shift;

  function Shift ( l : Laur_Sys ) return Laur_Sys is
    res : Laur_Sys (l'range);
  begin
    for k in l'range loop
      res(k) := Shift(l(k));
    end loop;
    return res;
  end Shift;

  procedure Transform ( t : in Transfo; p : in out Poly ) is

    procedure Transform_Term ( tt : in out Term; cont : out boolean ) is
    begin
      Apply(t,Link_to_Vector(tt.dg));
      cont := true;
    end Transform_Term;
    procedure Transform_Terms is new Changing_Iterator (Transform_Term);

  begin
    Transform_Terms(p);
  end Transform;

  function Transform ( t : Transfo; p : Poly ) return Poly is
    res : Poly;
  begin
    Copy(p,res);
    Transform(t,res);
    return res;
  end Transform;

  function  Transform2 ( t : Transfo; p : Poly ) return Poly is

  -- IMPORTANT : This function might change the term order !

    res : Poly := Null_Poly;

    procedure Transform_Term ( tt : in Term; cont : out boolean ) is
      rt : Term;
    begin
      rt.cf := tt.cf;
      rt.dg := Degrees(t*Link_to_Vector(tt.dg));
      Plus_Term(res,rt);
      Clear(rt);
      cont := true;
    end Transform_Term;
    procedure Transform_Terms is new Visiting_Iterator (Transform_Term);

  begin
    Transform_Terms(p);
    return res;
  end Transform2;

  procedure Transform ( t : in Transfo; l : in out Laur_Sys ) is
  begin
    for i in l'range loop
      Transform(t,l(i));
    end loop;
  end Transform;

  function  Transform ( t : Transfo; l : Laur_Sys ) return Laur_Sys is
    res : Laur_Sys(l'range);
  begin
    for i in l'range loop
      res(i) := Transform(t,l(i));
    end loop;
    return res;
  end Transform;

  function Maximal_Support ( p : Poly; v : Vector ) return integer is

    res : integer;
    first : boolean := true;

    procedure Scan_Term ( t : in Term; cont : out boolean ) is

      sp : integer := t.dg.all*v;

    begin
      if first 
       then res := sp; first := false;
       elsif sp > res
           then res := sp;
      end if;
      cont := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator (Scan_Term);

  begin
    Scan_Terms(p);
    return res;
  end Maximal_Support;

  function Maximal_Support ( p : Poly; v : Link_to_Vector ) return integer is
  begin
    return Maximal_Support(p,v.all);
  end Maximal_Support;

  procedure Face ( i,m : in integer; p : in out Poly ) is

    procedure Face_Term ( t : in out Term; cont : out boolean ) is
    begin
      if t.dg(i) /= m
       then t.cf := CMPLX(0.0);
      end if;
      cont := true;
    end Face_Term;
    procedure Face_Terms is new Changing_Iterator(Face_Term);

  begin
    Face_Terms(p);
  end Face;

  function Face ( i,m : integer; p : Poly ) return Poly is

    res : Poly;

  begin
    Copy(p,res);
    Face(i,m,res);
    return res;
  end Face;

  function  Face2 ( i,m : integer; p : Poly ) return Poly is

  -- IMPORTANT : This function might change the term order !

    res : Poly := Null_Poly;

    procedure Face_Term ( t : in Term; cont : out boolean ) is
    begin
      if t.dg(i) = m
       then Plus_Term(res,t);
      end if;
      cont := true;
    end Face_Term;
    procedure Face_Terms is new Visiting_Iterator(Face_Term);

  begin
    Face_Terms(p);
    return res;
  end Face2;

  procedure Face ( i,m : in integer; l : in out Laur_Sys ) is
  begin
    for j in l'range loop
      Face(i,m,l(j));
    end loop;
  end Face;

  function Face ( i,m : integer; l : Laur_Sys ) return Laur_Sys is

    res : Laur_Sys(l'range);

  begin
    for j in l'range loop
      res(j) := Face(i,m,l(j));
    end loop;
    return res;
  end Face;

  procedure Face ( v : in Vector; m : in integer; p : in out Poly ) is

    procedure Face_Term ( t : in out Term; cont : out boolean ) is
    begin
      if t.dg.all*v /= m
       then t.cf := CMPLX(0.0);
      end if;
      cont := true;
    end Face_Term;
    procedure Face_Terms is new Changing_Iterator(Face_Term);

  begin
    Face_Terms(p);
  end Face;

  function Face ( v : Vector; m : integer; p : Poly ) return Poly is

    res : Poly;

  begin
    Copy(p,res);
    Face(v,m,res);
    return res;
  end Face;

  function  Face2 ( v : Vector; m : integer; p : Poly ) return Poly is

  -- IMPORTANT : This procedure might change the term order !

    res : Poly := Null_Poly;

    procedure Face_Term ( t : in Term; cont : out boolean ) is
    begin
      if t.dg.all*v = m
       then Plus_Term(res,t);
      end if;
      cont := true;
    end Face_Term;
    procedure Face_Terms is new Visiting_Iterator(Face_Term);

  begin
    Face_Terms(p);
    return res;
  end Face2;

  procedure Face ( v,m : in Vector; l : in out Laur_Sys ) is
  begin
    for i in l'range loop
      Face(v,m(i),l(i));
    end loop;
  end Face;

  function Face ( v,m : Vector; l : Laur_Sys ) return Laur_Sys is

    res : Laur_Sys(l'range);

  begin
    for i in l'range loop
      res(i) := Face(v,m(i),l(i));
    end loop;
    return res;
  end Face;

  procedure Reduce ( i : in integer; p : in out Poly ) is

    procedure Reduce_Term ( t : in out Term; cont : out boolean ) is
    begin
      Reduce(Link_to_Vector(t.dg),i);
      cont := true;
    end Reduce_Term;
    procedure Reduce_Terms is new Changing_Iterator (Reduce_Term);

  begin
    Reduce_Terms(p);
  end Reduce;

  function Reduce ( i : integer; p : Poly ) return Poly is
    res : Poly;
  begin
    Copy(p,res);
    Reduce(i,res);
    return res;
  end Reduce;

  function  Reduce2 ( i : integer; p : Poly ) return Poly is

  -- IMPORTANT : This function might change the term order !

    res : Poly := Null_Poly;

    procedure Reduce_Term ( t : in Term; cont : out boolean ) is
      rt : Term;
    begin
      rt.cf := t.cf;
      rt.dg := Degrees(Reduce(Link_to_Vector(t.dg),i));
      Plus_Term(res,rt);
      Clear(rt);
      cont := true;
    end Reduce_Term;
    procedure Reduce_Terms is new Visiting_Iterator (Reduce_Term);

  begin
    Reduce_Terms(p);
    return res;
  end Reduce2;

  procedure Reduce ( i : in integer; l : in out Laur_Sys ) is
  begin
    for j in l'range loop
      Reduce(i,l(j));
    end loop;
  end Reduce;

  function  Reduce ( i : integer; l : Laur_Sys ) return Laur_Sys is
    res : Laur_Sys(l'range);
  begin
    for j in l'range loop
      res(j) := Reduce(i,l(j));
    end loop;
    return res;
  end Reduce;

  procedure Insert ( i,d : in integer; p : in out Poly ) is

    procedure Insert_Term ( t : in out Term; cont : out boolean ) is
    begin
      Insert(Link_to_Vector(t.dg),i,d);
      cont := true;
    end Insert_Term;
    procedure Insert_Terms is new Changing_Iterator (Insert_Term);

  begin
    Insert_Terms(p);
  end Insert;

  function Insert ( i,d : integer; p : Poly ) return Poly is
    res : Poly;
  begin
    Copy(p,res);
    Insert(i,d,res);
    return res;
  end Insert;

  function  Insert2 ( i,d : integer; p : Poly ) return Poly is

  -- IMPORTANT : This function might change the term order !

    res : Poly := Null_Poly;

    procedure Insert_Term ( t : in Term; cont : out boolean ) is
      rt : Term;
    begin
      rt.cf := t.cf;
      rt.dg := Degrees(Insert(Link_to_Vector(t.dg),i,d));
      Plus_Term(res,rt);
      Clear(rt);
      cont := true;
    end Insert_Term;
    procedure Insert_Terms is new Visiting_Iterator (Insert_Term);

  begin
    Insert_Terms(p);
    return res;
  end Insert2;

  procedure Insert ( i,d : in integer; l : in out Laur_Sys ) is
  begin
    for j in l'range loop
      Insert(i,d,l(j));
    end loop;
  end Insert;

  function  Insert ( i,d : integer; l : Laur_Sys ) return Laur_Sys is
    res : Laur_Sys(l'range);
  begin
    for j in l'range loop
      res(j) := Insert(i,d,l(j));
    end loop;
    return res;
  end Insert;

end Transforming_Laurent_Systems;
