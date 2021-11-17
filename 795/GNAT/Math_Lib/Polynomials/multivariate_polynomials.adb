with unchecked_deallocation;
with Lists;

package body Multivariate_Polynomials is

-- REPRESENTATION INVARIANT :
--   1. Only terms with a coefficient different from zero are stored.
--   2. The terms in any polynomial are ordered from high to low degree
--      according to the instantiating order on the degrees.

-- DATA STRUCTURES : 

  package Term_List is new Lists(Term);
  type Poly_Rep is new Term_List.List;

  type kind is (polynomial,number);
  type Poly_Rec ( k : kind := number ) is record
    case k is
      when number     => c : coefftp;
      when polynomial => p : Eval_Poly;
    end case;
  end record;
  type Coeff_Poly_Rec ( k : kind := number ) is record
    case k is
      when number     => c : integer;
      when polynomial => p : Eval_Coeff_Poly;
    end case;
  end record;

  Null_Poly_Rec : constant Poly_Rec(number) := (number,zero);
  type Eval_Poly_Rep is array(integer range <>) of Poly_Rec;
  type Eval_Coeff_Poly_Rep is array(integer range <>) of Coeff_Poly_Rec;

  procedure free is new unchecked_deallocation(Poly_Rep,Poly);
  procedure free is new unchecked_deallocation(Eval_Poly_Rep,Eval_Poly);
  procedure free is 
    new unchecked_deallocation(Eval_Coeff_Poly_Rep,Eval_Coeff_Poly);
 
-- AUXILIARY OPERATIONS :

  function Convert ( c : coefftp; n : natural ) return natural is

  -- DESCRIPTION :
  --   Returns the corresponding value for c, when it lies in 1..n,
  --   otherwise 0 is returned.

  begin
    for i in 1..n loop
      if c = Convert(i)
       then return i;
      end if;
    end loop;
    return 0;
  end Convert;

  procedure Shuffle ( p : in out Poly ) is

  -- DESCRIPTION :
  --   Changes the position of the terms in p back to the normal order.
  --   Needed to guarantee the second representation invariant.

    res : Poly := Null_Poly;

    procedure Shuffle_Term ( t : in Term; cont : out boolean ) is
    begin
      Plus_Term(res,t);
      cont := true;
    end Shuffle_Term;
    procedure Shuffle_Terms is new Visiting_Iterator(Shuffle_Term);

  begin
    Shuffle_Terms(p);
    Clear(p); Copy(res,p); Clear(res);
  end Shuffle;

  procedure Append_Copy ( first,last : in out Poly_Rep; t : in Term ) is

  -- DESCRIPTION :
  --   Appends a copy of the term to the list.

    tt : Term;

  begin
    Copy(t,tt);
    Append(first,last,tt);
  end Append_Copy;

-- CONSTRUCTORS :
 
  function Create ( t : Term ) return Poly is

    p : Poly;

  begin
    if equal(t.cf,zero)
     then p := Null_Poly;
     else declare
            tt : Term;
          begin
            Copy(t,tt);
            p := new Poly_Rep;
            Construct(tt,p.all);
          end;
    end if;
    return p;
  end Create;

  function Create ( p : Poly; n : natural; d : integer ) return Eval_Poly is

  -- DESCRIPTION :
  --   An evaluable polynomial is returned for p, with d = Degree(p,x1)
  --   and n = Number_of_Unknowns(p) >= 1.

    res : Eval_Poly;
    evpr : Eval_Poly_Rep(0..d);
    terms : array(0..d) of Poly := (0..d => Null_Poly);

    procedure Add_Term1 ( t : in Term; cont : out boolean ) is

      pr : Poly_Rec(number);

    begin
      copy(t.cf,pr.c);
      evpr(t.dg(t.dg'first)) := pr;
      cont := true;
    end Add_Term1;
    procedure Add_Terms1 is new Visiting_Iterator(Add_Term1);

    procedure Add_Term ( t : in Term; cont : out boolean ) is

      nt : Term;

    begin
      nt.cf := t.cf;
      nt.dg := new Natural_Vectors.Vector(t.dg'first..t.dg'last-1);
      for i in nt.dg'range loop
        nt.dg(i) := t.dg(i+1);
      end loop;
      Plus_Term(terms(t.dg(t.dg'first)),nt);
      cont := true;
    end Add_Term;
    procedure Add_Terms is new Visiting_Iterator(Add_Term);

  begin
    if n = 1
     then Add_Terms1(p);
     else Add_Terms(p);
          for i in terms'range loop
            declare
              pr : Poly_Rec(polynomial);
            begin
              pr.p := Create(terms(i));
              evpr(i) := pr;
              Clear(terms(i));
            end;
          end loop;
    end if;
    res := new Eval_Poly_Rep'(evpr);
    return res;
  end Create;

  function Create ( p : Poly; n,nb : natural; d : integer )
                  return Eval_Coeff_Poly is

  -- DESCRIPTION :
  --   An evaluable polynomial is returned for p, with d = Degree(p,x1),
  --   n = Number_of_Unknowns(p) >= 1 and nb = Number_of_Terms(p).
  --   The coefficients of p are converted natural numbers.

    res : Eval_Coeff_Poly;
    evpr : Eval_Coeff_Poly_Rep(0..d);
    terms : array(0..d) of Poly := (0..d => Null_Poly);

    procedure Add_Term1 ( t : in Term; cont : out boolean ) is

      pr : Coeff_Poly_Rec(number);

    begin
      pr.c := Convert(t.cf,nb);
      evpr(t.dg(t.dg'first)) := pr;
      cont := true;
    end Add_Term1;
    procedure Add_Terms1 is new Visiting_Iterator(Add_Term1);

    procedure Add_Term ( t : in Term; cont : out boolean ) is

      nt : Term;

    begin
      nt.cf := t.cf;
      nt.dg := new Natural_Vectors.Vector(t.dg'first..t.dg'last-1);
      for i in nt.dg'range loop
        nt.dg(i) := t.dg(i+1);
      end loop;
      Plus_Term(terms(t.dg(t.dg'first)),nt);
      cont := true;
    end Add_Term;
    procedure Add_Terms is new Visiting_Iterator(Add_Term);

  begin
    if n = 1
     then for i in evpr'range loop                   -- initialization
            declare
              nullpr : Coeff_Poly_Rec(polynomial);
            begin
              nullpr.p := null;
              evpr(i) := nullpr;
            end;
          end loop;
          Add_Terms1(p);
     else Add_Terms(p);
          for i in terms'range loop
            declare
              pr : Coeff_Poly_Rec(polynomial);
              ind : integer;
            begin
              if terms(i) = Null_Poly
               then pr.p := null;
               else ind := Head_Of(terms(i).all).dg'first;
                    pr.p := Create(terms(i),n-1,nb,Degree(terms(i),ind));
              end if;
              evpr(i) := pr;
              Clear(terms(i));
            end;
          end loop;
    end if;
    res := new Eval_Coeff_Poly_Rep'(evpr);
    return res;
  end Create;

  function Create ( p : Poly ) return Eval_Poly is

    n : constant natural := Number_of_Unknowns(p);
    ind : integer;

  begin
    if (p = Null_Poly) or else (n = 0)
     then return null;
     else ind := Head_Of(p.all).dg'first;
          return Create(p,n,Degree(p,ind));
    end if;
  end Create;

  function Create ( p : Poly ) return Eval_Coeff_Poly is

    res : Eval_Coeff_Poly;
    lp : Poly := Null_Poly;
    n : constant natural := Number_of_Unknowns(p);
    nb : constant natural := Number_of_Terms(p);
    cnt : natural := 0;
    ind : integer;

    procedure Label_Term ( t : in Term; cont : out boolean ) is

      lt : Term;

    begin
      cnt := cnt + 1;
      lt.cf := Convert(cnt);
      lt.dg := new Natural_Vectors.Vector'(t.dg.all);
      Plus_Term(lp,lt);
      cont := true;
    end Label_Term;
    procedure Label_Terms is new Visiting_Iterator(Label_Term);

  begin
    if (p = Null_Poly) or else (nb = 0)
     then return null;
     else Label_Terms(p);
          ind := Head_Of(p.all).dg'first;
          res := Create(lp,n,nb,Degree(p,ind));
          Clear(lp);
    end if;
    return res;
  end Create;

  procedure Copy ( t1 : in Term; t2 : in out Term ) is
  begin
    Clear(t2);
    Natural_Vectors.Copy(Link_to_Vector(t1.dg),Link_to_Vector(t2.dg));
    copy(t1.cf,t2.cf);
  end Copy;

  procedure Copy ( p: in Poly_Rep; q : in out Poly_Rep ) is
 
    lp,lq : Poly_Rep;
    t : Term;
 
  begin
    Clear(q);
    if not Is_Null(p)
     then lp := p;
          while not Is_Null(lp) loop
            t := Head_Of(lp);
            Append_Copy(q,lq,t);
            lp := Tail_Of(lp);
          end loop;
    end if;
  end Copy;

  procedure Copy ( p : in Poly; q : in out Poly ) is

    l : Poly_Rep;

  begin
    Clear(q);
    if p /= Null_Poly
     then Copy(p.all,l);
          q := new Poly_Rep'(l);
     else q := Null_Poly;
    end if;
  end Copy;

-- SELECTORS :

  function Equal ( t1,t2 : Term ) return boolean is
  begin
    if not Equal(t1.dg,t2.dg)
     then return false;
     else return equal(t1.cf,t2.cf);
    end if;
  end Equal;

  function Equal ( p,q : Poly_Rep ) return boolean is

    lp, lq : Poly_Rep;

  begin
    lp := p;
    lq := q;
    while not Is_Null(lp) and not Is_Null(lq) loop
      if not Equal(Head_Of(lp),Head_Of(lq))
       then return false;
       else lp := Tail_Of(lp);
            lq := Tail_Of(lq);
      end if;
    end loop;
    if Is_Null(lp) and Is_Null(lq)
     then return true;
     else return false;
    end if;
  end Equal;

  function Equal ( p,q : Poly ) return boolean is
  begin
    if (p = Null_Poly) and (q = Null_Poly)
     then return true;
      elsif (p = Null_Poly) or (q = Null_Poly)
         then return false;
         else return Equal(p.all,q.all);
    end if;
  end Equal;
 
  function Number_Of_Unknowns ( p : Poly ) return natural is

    temp : Term;

  begin
    if (p = Null_Poly) or else Is_Null(p.all)
     then return 0;
     else temp := Head_Of(p.all);
          if temp.dg = null
           then return 0;
           else return temp.dg'length;
          end if;
    end if;
  end Number_Of_Unknowns;

  function Number_Of_Terms ( p : Poly ) return natural is
  begin
     if (p = Null_Poly) or else Is_Null(p.all)
      then return 0;
      else return Length_Of(p.all);
     end if;
  end Number_Of_Terms;
 
  function Degree ( p : Poly ) return integer is

    temp : Term;

  begin
     if (p = Null_Poly) or else Is_Null(p.all)
      then return -1;
      else temp := Head_Of(p.all);
           if temp.dg = null
            then return 0;
            else return integer(Sum(temp.dg));
           end if;
     end if;
  end Degree;
 
  function Degree ( p : Poly; i : integer ) return integer is

    res : integer := 0;
   
    procedure Degree_Term (t : in Term; continue : out boolean) is
      index,temp : integer;
    begin
      if t.dg /= null
       then index := t.dg'first + i - 1;
            temp := t.dg(index);
            if (temp > 0) and then (temp > res)
             then res := temp;
            end if;
            continue := true;
      end if;
    end Degree_Term;
    procedure Degree_Terms is new Visiting_Iterator(process => Degree_Term);

  begin
    if p = Null_Poly or else Is_Null(p.all)
     then return -1;
     else Degree_Terms(p);
          return res;
    end if;
  end Degree;

  function "<" ( d1,d2 : Degrees ) return boolean is
  begin
    return Link_to_Vector(d1) < Link_to_Vector(d2);
  end "<";

  function ">" ( d1,d2 : Degrees ) return boolean is
  begin
    return Link_to_Vector(d1) > Link_to_Vector(d2);
  end ">";

  function Coeff ( p : Poly; d : degrees ) return coefftp is

    l : Poly_Rep;
    t : term;
    res : coefftp;

  begin 
    if p = Null_Poly
     then return zero; 
     else l := p.all;
          while not Is_Null(l) loop
            t := Head_Of(l);
            if t.dg < d
             then return zero;
             elsif Equal(t.dg,d)
                 then copy(t.cf,res);
                      return res;
                 else l := Tail_Of(l);
            end if;
          end loop;
          return zero;
    end if;
  end Coeff;

  function Coeff ( p : Poly ) return Vector_of_coefftp is

    res : Vector_of_coefftp(1..Number_of_Terms(p));
    cnt : natural := 0;

    procedure Visit_Term ( t : in Term; cont : out boolean ) is
    begin
      cnt := cnt + 1;
      copy(t.cf,res(cnt));
      cont := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Coeff;

-- ARITHMETICAL OPERATIONS :
 
  function "+" ( p : Poly; t : Term ) return Poly is

    temp : Poly;

  begin
    Copy(p,temp);
    Plus_Term(temp,t);
    return temp;
  end "+";

  function "+" ( t : Term; p : Poly ) return Poly is
  begin
    return p+t;
  end "+";
 
  function "+" ( p,q : Poly ) return Poly is

    temp : Poly;

  begin
    Copy(p,temp);
    Plus_Poly(temp,q);
    return temp;
  end "+";
 
  function "-" ( p : Poly; t : Term ) return Poly is

    temp : Poly;

  begin
    Copy(p,temp);
    Min_Term(temp,t);
    return temp;
  end "-";
 
  function "-" ( t : Term; p : Poly ) return Poly is

    temp : Poly;

  begin
    temp := Create(t);
    Min_Poly(temp,p);
    return temp;
  end "-";
   
  function "-" ( p : Poly ) return Poly is

    temp : Poly;

  begin
    Copy(p,temp);
    Min_Poly(temp);
    return temp;
  end "-";
 
  function "-" ( p,q : Poly ) return Poly is
 
    temp : Poly;

  begin
    Copy(p,temp);
    Min_Poly(temp,q);
    return temp;
  end "-";
 
  function "*" ( p : Poly; a : coefftp ) return Poly is

    temp : Poly;

  begin
    Copy(p,temp);
    Mult_Coeff(temp,a);
    return temp;
  end "*";
 
  function "*" ( a : coefftp; p : Poly ) return Poly is
  begin
    return p*a;
  end "*";
 
  function "*" ( p : Poly; t : Term ) return Poly is

    temp : Poly;

  begin
    Copy(p,temp);
    Mult_Term(temp,t);
    return temp;
  end "*";
 
  function "*" ( t : Term; p : Poly ) return Poly is
  begin
    return p*t;
  end "*";
 
  function "*" ( p,q : Poly ) return Poly is

    temp : Poly;

  begin
    Copy(p,temp);
    Mult_Poly(temp,q);
    return temp;
  end "*";
   
  procedure Plus_Term ( p : in out Poly; t : in Term ) is

    l1,l2,temp : Poly_Rep;
    tt,tp : Term;

  begin
    if t.cf /= zero
     then Copy(t,tt);
          if p = Null_Poly
           then p := new Poly_Rep;
                Construct(tt,p.all);
           else tp := Head_Of(p.all);
                if tt.dg > tp.dg
                 then Construct(tt,p.all);
                 elsif Equal(tt.dg,tp.dg)
                     then Plus_Coeff(tp.cf,tt.cf);
                          if tp.cf /= zero
                           then Set_Head(p.all,tp);
                           else Clear(tp);
                                if Is_Null(Tail_Of(p.all))
                                 then Term_List.Clear(Term_List.List(p.all));
                                      free(p);
                                 else Swap_Tail(p.all,l1);
                                      Term_List.Clear(Term_List.List(p.all));
                                      p.all := l1;
                                end if;
                          end if;
                          Clear(tt);
                     else l1 := p.all;
                          l2 := Tail_Of(l1);
                          while not Is_Null(l2) loop
                            tp := Head_Of(l2);
                            if tt.dg > tp.dg
                             then Construct(tt,temp);
                                  Swap_Tail(l1,temp);
                                  l1 := Tail_Of(l1);
                                  Swap_Tail(l1,temp);
                                  return;
                             elsif Equal(tt.dg,tp.dg)
                                 then Plus_Coeff(tp.cf,tt.cf);
                                      if tp.cf /= zero
                                       then Set_Head(l2,tp);
                                       else Clear(tp);
                                            temp := Tail_Of(l2);
                                            Swap_Tail(l1,temp);
                                      end if;
                                      Clear(tt);
                                      return;
                                 else l1 := l2;
                                      l2 := Tail_Of(l1);
                            end if;
                          end loop;
                          Construct(tt,temp);
                          Swap_Tail(l1,temp);
                end if;
          end if;
    end if;
  end Plus_Term;

  procedure Plus_Poly ( p : in out Poly; q : in Poly ) is

    procedure Plus_Term ( t : in Term; continue : out boolean ) is
    begin
      Plus_Term(p,t);
      continue := true;
    end Plus_Term;
    procedure Plus_Terms is new Visiting_Iterator(Plus_Term);
 
  begin
    Plus_Terms(q);
  end Plus_Poly;
 
  procedure Min_Term ( p : in out Poly; t : in Term ) is

    temp : Term;

  begin
    Natural_Vectors.Copy(Link_to_Vector(t.dg),Link_to_Vector(temp.dg));
    temp.cf := -t.cf;
    Plus_Term(p,temp);
    Natural_Vectors.Clear(Link_to_Vector(temp.dg));
    clear(temp.cf);
  end Min_Term;

  procedure Min_Poly ( p : in out Poly; q : in Poly ) is

    temp : Poly := -q;

  begin
    Plus_Poly(p,temp);
    Clear(temp);
  end Min_Poly;
 
  procedure Min_Poly ( p : in out Poly ) is
 
    procedure Min_Term ( t : in out Term; continue : out boolean ) is
    begin
      Min_Coeff(t.cf);
      continue := true;
    end Min_Term;
    procedure Min_Poly1 is new Changing_Iterator (process => Min_Term);
 
  begin
    Min_Poly1(p);
  end Min_Poly;
 
  procedure Mult_Coeff ( p : in out Poly; a : in coefftp ) is
 
    procedure Mult_Term ( t : in out Term; continue : out boolean ) is
    begin
      Mult_Coeff(t.cf,a);
      continue := true;
    end Mult_Term;
    procedure Mult_Coeff1 is new Changing_Iterator (process => Mult_Term);
 
  begin
    Mult_Coeff1(p);
  end Mult_Coeff;

  procedure Mult_Term ( p : in out Poly; t : in Term ) is

    procedure Mult_Term ( tp : in out Term; continue : out boolean ) is
    begin
      Natural_Vectors.Plus_Vector(Link_to_Vector(tp.dg),Link_to_Vector(t.dg));
      Mult_Coeff(tp.cf,t.cf);
      continue := true;
    end Mult_Term;
    procedure Mult_Terms is new Changing_Iterator (process => Mult_Term);

  begin
    Mult_Terms(p);
  end Mult_Term;
 
  procedure Mult_Poly ( p : in out Poly; q : in Poly ) is

    res : Poly;
    l : Poly_Rep;
    t : Term;

  begin
    if (q = Null_Poly) or else Is_Null(q.all)
     then Clear(p);
     else l := q.all;
          res := Null_Poly;
          while not Is_Null(l) loop
            t := Head_Of(l);
            declare
              temp : Poly;
            begin
              temp := p * t;
              Plus_Poly(res,temp);
              Clear(temp);
            end;
            l := Tail_Of(l);
          end loop;
          Copy(res,p); Clear(res);
    end if;
  end Mult_Poly;
 
  function Eval ( p : Poly; x : coefftp; i : integer ) return Poly is

    res : Poly := Null_Poly;

    procedure Eval_Term ( t : in Term; cont : out boolean ) is

      nt : Term;

    begin
      copy(t.cf,nt.cf);
      nt.dg := new Natural_Vectors.Vector(t.dg'first..t.dg'last-1);
      for j in t.dg'range loop
        if j < i
         then nt.dg(j) := t.dg(j);
         elsif j > i
             then nt.dg(j-1) := t.dg(j);
             else for k in 1..t.dg(i) loop
                    mult_coeff(nt.cf,x);
                  end loop;
        end if;
      end loop;
      Plus_Term(res,nt);
      cont := true;
    end Eval_Term;
    procedure Eval_Terms is new Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Eval;

  function Eval ( d : Degrees; c : coefftp; x : Vector_of_coefftp )
                return coefftp is

    res : coefftp;

  begin
    copy(c,res);
    for i in d'range loop
      for j in 1..d(i) loop
        mult_coeff(res,x(i));
      end loop;
    end loop;
    return res;
  end Eval;

  function Eval ( t : Term; x : Vector_of_coefftp ) return coefftp is
  begin
    return Eval(t.dg,t.cf,x);
  end Eval;

  function Eval ( t : Term; c : coefftp; x : Vector_of_coefftp )
                return coefftp is
  begin
    return Eval(t.dg,c,x);
  end Eval;
 
  function Eval ( p : Poly; x : Vector_of_coefftp ) return coefftp is

    res : coefftp := zero;

    procedure Eval_Term ( t : in Term; cont : out boolean ) is
    begin
      plus_coeff(res,Eval(t,x));
      cont := true;
    end Eval_Term;
    procedure Eval_Terms is new Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Eval;

  function Eval ( p : Poly; c,x : Vector_of_coefftp ) return coefftp is

    res : coefftp := zero;
    cnt : natural := 0;

    procedure Eval_Term ( t : in Term; cont : out boolean ) is
    begin
      cnt := cnt + 1;
      plus_coeff(res,Eval(t,c(cnt),x));
      cont := true;
    end Eval_Term;
    procedure Eval_Terms is new Visiting_Iterator(Eval_Term);

  begin
    Eval_Terms(p);
    return res;
  end Eval;
 
  function Eval ( vp : Eval_Poly_Rep; x : Vector_of_coefftp; i : integer )
                return coefftp;
 
  function Eval ( vprec : Poly_Rec; x : Vector_of_coefftp; i : integer )
                return coefftp is

    res : coefftp;

  begin
    if vprec.k = number
     then copy(vprec.c,res);
          return res;
     else return Eval(vprec.p.all,x,i);
    end if;
  end Eval;
 
  function Eval ( vp : Eval_Poly_Rep; x : Vector_of_coefftp; i : integer )
                return coefftp is

    deg : natural := vp'length-1;
    res : coefftp;

  begin
    if deg = 0
     then return Eval(vp(0),x,i+1);
     else res := Eval(vp(deg),x,i+1);
          for j in reverse 0..(deg-1) loop
            Mult_Coeff(res,x(i));
            if (vp(j).k = number) or (vp(j).p /= null)
             then declare
                    temp : coefftp := Eval(vp(j),x,i+1);
                  begin
                    Plus_Coeff(res,temp);
                    clear(temp);
                  end;
            end if;
          end loop;
          return res;
    end if;
  end Eval;
 
  function Eval ( p : Eval_Poly; x : Vector_of_coefftp ) return coefftp is
  begin
    if p = null
     then return zero;
     else return Eval(p.all,x,x'first);
    end if;
  end Eval;

  function Eval ( vp : Eval_Coeff_Poly_Rep; c,x : Vector_of_coefftp;
                  i : integer ) return coefftp;

  function Eval ( vprec : Coeff_Poly_Rec; c,x : Vector_of_coefftp; i : integer )
                return coefftp is

    res : coefftp;

  begin
    if vprec.k = number
     then copy(c(vprec.c),res);
          return res;
     else return Eval(vprec.p.all,c,x,i);
    end if;
  end Eval;

  function Eval ( vp : Eval_Coeff_Poly_Rep; c,x : Vector_of_coefftp;
                  i : integer ) return coefftp is

    deg : natural := vp'length-1;
    res : coefftp;

  begin
    if deg = 0
     then return Eval(vp(0),c,x,i+1);
     else res := Eval(vp(deg),c,x,i+1);
          for j in reverse 0..(deg-1) loop
            Mult_Coeff(res,x(i));
            if (vp(j).k = number) or (vp(j).p /= null)
             then declare
                    temp : coefftp := Eval(vp(j),c,x,i+1);
                  begin
                    Plus_Coeff(res,temp);
                    clear(temp);
                  end;
            end if;
          end loop;
          return res;
    end if;
  end Eval;

  function Eval ( p : Eval_Coeff_Poly; c,x : Vector_of_coefftp )
                return coefftp is
  begin
    if p = null
     then return zero;
     else return Eval(p.all,c,x,x'first);
    end if;
  end Eval;
 
  procedure Diff ( p : in out Poly; i : in integer ) is

    procedure Diff_Term ( t : in out Term; continue : out boolean ) is
      temp : coefftp;
      index : integer := t.dg'first + i - 1;
    begin
      if t.dg(index) = 0
       then Clear(t);
            t.cf := zero;
       else temp := convert(t.dg(index));
            Mult_Coeff(t.cf,temp);
            t.dg(index) := t.dg(index) - 1;
      end if;
      continue := true;
    end Diff_Term;
    procedure Diff_Terms is new Changing_Iterator( process => Diff_Term );

  begin
    Diff_Terms(p);
  end Diff;

  function Diff ( p : Poly; i : integer ) return Poly is

    temp : Poly;

  begin
    Copy(p,temp);
    Diff(temp,i);
    return temp;
  end Diff;

  procedure Diff ( p : in Poly; i : in integer;
                   cp : out Eval_Coeff_Poly; m : out Vector_of_coefftp ) is

    nb : constant natural := Number_of_Terms(p);
    n : constant natural := Number_of_Unknowns(p);
    ind,cnt : integer;
    dp : Poly := Null_Poly;

    procedure Diff_Term ( t : in Term; cont : out boolean ) is

      dt : Term;

    begin
      cnt := cnt + 1;
      if t.dg(i) > 0
       then dt.cf := convert(cnt);
            dt.dg := new Natural_Vectors.Vector'(t.dg.all);
            m(cnt) := convert(t.dg(i));
            dt.dg(i) := dt.dg(i) - 1;
            Plus_Term(dp,dt);
            Clear(dt);
       else m(cnt) := convert(0);
      end if;
      cont := true;
    end Diff_Term;
    procedure Diff_Terms is new Visiting_Iterator(Diff_Term);

  begin
    cnt := 0;
    Diff_Terms(p);
    if dp /= null
     then ind := Head_Of(dp.all).dg'first;
          cp := Create(dp,n,nb,Degree(dp,ind));
    end if;
  end Diff;

-- ITERATORS :

  procedure Visiting_Iterator ( p : in Poly ) is

    l : Poly_Rep;
    temp : Term;
    continue : boolean;

  begin
    if p /= Null_Poly
     then l := p.all;
          while not Is_Null(l) loop 
            temp := Head_Of(l);
            process(temp,continue);
            exit when not continue;
            l := Tail_Of(l);
          end loop;
    end if;
  end Visiting_Iterator;

  procedure Changing_Iterator ( p : in out Poly ) is

    q,lq,lp : Poly_Rep;
    t : Term;
    continue : boolean := true;

    procedure copy_append is
      temp : Term;
    begin
      Copy(t,temp);
      if continue
       then process(temp,continue);
      end if;
      if temp.cf /= zero
       then Append(q,lq,temp);
       else Clear(temp);
      end if;
      Clear(t);
    end copy_append;

  begin
    if p = Null_Poly
     then return;
     else lp := p.all;
          while not Is_Null(lp) loop
            t := Head_Of(lp);
            copy_append;
            lp := Tail_Of(lp);
          end loop;
          Term_List.Clear(Term_List.List(p.all));  free(p);
          if Is_Null(q)
           then p := Null_Poly;
           else p := new Poly_Rep'(q);
          end if;
          Shuffle(p);  -- ensure the second representation invariant
    end if;
  end Changing_Iterator;

-- DESTRUCTORS :

  procedure Clear ( t : in out Term ) is
  begin
    Natural_Vectors.Clear(Link_to_Vector(t.dg));
    clear(t.cf);
  end Clear;
 
  procedure Clear ( p : in out Poly_Rep ) is

    l : Poly_Rep := p;
    t : Term;

  begin
    if Is_Null(p)
     then return;
     else while not Is_Null(l) loop
            t := Head_Of(l);
            Clear(t);
            l := Tail_Of(l);
          end loop;
          Term_List.Clear(Term_List.List(p));
    end if;
  end Clear;
 
  procedure Clear ( p : in out Poly ) is
  begin
    if p /= Null_Poly
     then Clear(p.all);
          free(p);
          p := Null_Poly;
    end if;
  end Clear;

  procedure Clear ( p : in out Eval_Poly ) is
  begin
    if p /= null
     then declare
            vp : Eval_Poly_Rep renames p.all;
          begin
            for i in vp'range loop
              if vp(i).k = number
               then clear(vp(i).c);
               else Clear(vp(i).p);
              end if;
            end loop;
          end;
          free(p);
    end if;
  end Clear;

  procedure Clear ( p : in out Eval_Coeff_Poly ) is
  begin
    if p /= null
     then declare
            vp : Eval_Coeff_Poly_Rep renames p.all;
          begin
            for i in vp'range loop
              if vp(i).k /= number
               then Clear(vp(i).p);
              end if;
            end loop;
          end;
          free(p);
    end if;
  end Clear;

end Multivariate_Polynomials;
