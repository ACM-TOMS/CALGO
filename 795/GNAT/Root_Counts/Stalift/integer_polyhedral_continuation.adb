with integer_io;                         use integer_io;
with Floating_Point_Numbers;             use Floating_Point_Numbers;
with Mathematical_Functions;             use Mathematical_Functions;
with Float_Vectors;
with Float_Vectors_of_Vectors;
with Integer_Vectors_of_Vectors;
with Complex_Multivariate_Polynomials;
with Complex_Multivariate_Laurent_Polynomials;
with Complex_Polynomial_Systems;         use Complex_Polynomial_Systems;
with Laurent_to_Polynomial_Converters;   use Laurent_to_Polynomial_Converters;

with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Power_Lists;                        use Power_Lists;
with Integer_Lifting_Utilities;          use Integer_Lifting_Utilities;

with Substitutors;                       use Substitutors;
with Transforming_Integer_Vector_Lists;  use Transforming_Integer_Vector_Lists;
with Transforming_Laurent_Systems;       use Transforming_Laurent_Systems;

with Complex_Matrices;                   use Complex_Matrices;
with Complex_Numbers,Complex_Norms;      use Complex_Numbers,Complex_Norms;
with Complex_Vectors;
with Continuation_Parameters;
with Increment_and_Fix_Continuation;     use Increment_and_Fix_Continuation;

with Fewnomials;
with BKK_Bound_Computations;             use BKK_Bound_Computations;

with Integer_Vectors_of_Vectors_io;      use Integer_Vectors_of_Vectors_io;
with Integer_Vectors_io;                 use Integer_Vectors_io;
with Complex_Numbers_io;                 use Complex_Numbers_io;

package body Integer_Polyhedral_Continuation is

  procedure Write ( file : in file_type;
                    c : in Complex_Vectors_of_Vectors.Vector ) is
  begin
    for i in c'range loop
      put(file,i,1); put(file," : ");
      for j in c(i)'range loop
        if REAL_PART(c(i)(j)) = 0.0
         then put(file," 0");
         else put(file,c(i)(j),2,3,0);
        end if;
      end loop;
      new_line(file);
    end loop;
  end Write;

-- UTILITIES FOR POLYHEDRAL COEFFICIENT HOMOTOPY :

  function Power ( x,m : double_float ) return double_float is

  -- DESCRIPTION :
  --   Computes x**m for high powers of m to avoid RANGE_ERROR.

   -- intm : integer := integer(m);
   -- fltm : double_float;

  begin
   -- if m < 10.0
   --  then
    return (x**m);  -- FOR GNAT 3.07
   --  else if double_float(intm) > m
   --        then intm := intm-1;
   --       end if;
   --       fltm := m - double_float(intm);
   --       return ((x**intm)*(x**fltm));
   -- end if;
  end Power;

  function Minimum ( v : Integer_Vectors.Vector ) return integer is

  -- DESCRIPTION :
  --   Returns the minimal element in the vector v, different from zero.

    tmp,min : integer := 0;

  begin
    for i in v'range loop
      if v(i) /= 0 
       then if v(i) < 0
             then tmp := -v(i);
            end if;
            if min = 0
             then min := tmp;
             elsif tmp < min
                 then min := tmp;
            end if;
      end if;
    end loop;
    return min;
  end Minimum;

  function Scale ( v : Integer_Vectors.Vector ) return Float_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the vector v divided by its minimal element different from zero,
  --   such that the smallest positive element in the scaled vector equals one.

    res : Float_Vectors.Vector(v'range);
    min : constant integer := Minimum(v);

  begin
    if (min = 0) or (min = 1)
     then for i in res'range loop
            res(i) := double_float(v(i));
          end loop;
     else for i in res'range loop
            res(i) := double_float(v(i))/double_float(min);
          end loop;
    end if;
    return res;
  end Scale;

  function Scale ( v : Integer_Vectors_of_Vectors.Vector )
                 return Float_Vectors_of_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns an array of scaled vectors.

    res : Float_Vectors_of_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      declare
        sv : constant Float_Vectors.Vector := Scale(v(i).all);
      begin
        res(i) := new Float_Vectors.Vector(sv'range);
        for j in sv'range loop
          res(i)(j) := sv(j);
        end loop;
      end;  -- detour set up for GNAT 3.07
     -- res(i) := new Float_Vectors.Vector'(Scale(v(i).all));
    end loop;
    return res;
  end Scale;

  procedure Shift ( v : in out Integer_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Shifts the elements in v, such that the minimal element equals zero.

    min : integer := v(v'first);

  begin
    for i in v'first+1..v'last loop
      if v(i) < min
       then min := v(i);
      end if;
    end loop;
    if min /= 0
     then for i in v'range loop
            v(i) := v(i) - min;
          end loop;
    end if;
  end Shift;

  function Create ( e : Integer_Vectors_of_Vectors.Vector;
                    l : List; normal : Vector ) 
                  return Integer_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector with all inner products of the normal with
  --   the exponents in the list, such that minimal value equals zero.

    res : Integer_Vectors.Vector(e'range);
    lei : Integer_Vectors.Vector(normal'first..normal'last-1);
    found : boolean;
    lif : integer;

  begin
    for i in e'range loop
      lei := e(i).all;
      Search_Lifting(l,lei,found,lif);
      if not found
       then res(i) := 1000;
       else res(i) := lei*normal(lei'range) + normal(normal'last)*lif;
      end if;
    end loop;
    Shift(res);
    return res;
  end Create;

  function Create ( e : Exponent_Vectors_Array;
                    l : Array_of_Lists; mix,normal : Vector )
                  return Integer_Vectors_of_Vectors.Vector is

    res : Integer_Vectors_of_Vectors.Vector(e'range);
    cnt : natural := res'first;
 
  begin
    for i in mix'range loop
      declare
        rescnt : constant Integer_Vectors.Vector
               := Create(e(cnt).all,l(i),normal);
      begin
        res(cnt) := new Integer_Vectors.Vector(rescnt'range);
        for j in rescnt'range loop
          res(cnt)(j) := rescnt(j);
        end loop;      
      end;
      Shift(res(cnt).all);
      for k in 1..(mix(i)-1) loop
        res(cnt+k) := new Integer_Vectors.Vector(res(cnt)'range);
        for j in res(cnt)'range loop
          res(cnt+k)(j) := res(cnt)(j);
        end loop;
      end loop;
      cnt := cnt + mix(i);
    end loop;
    return res;
  end Create;

  procedure Eval ( c : in Complex_Vectors.Vector;
                   t : in double_float; m : in Integer_Vectors.Vector;
                   ctm : in out Complex_Vectors.Vector ) is

  -- DESCRIPTION : ctm = c*t**m.

  begin
    for i in ctm'range loop
      ctm(i) := c(i)*CMPLX((t**m(i)));
    end loop;
  end Eval;

  procedure Eval ( c : in Complex_Vectors.Vector;
                   t : in double_float; m : in Float_Vectors.Vector;
                   ctm : in out Complex_Vectors.Vector ) is

  -- DESCRIPTION : ctm = c*t**m.

  begin
    for i in ctm'range loop
     -- ctm(i) := c(i)*CMPLX((t**m(i)));
      if (REAL_PART(c(i)) = 0.0) and (IMAG_PART(c(i)) = 0.0)
       then ctm(i) := CMPLX(0.0);
       else ctm(i) := c(i)*CMPLX(Power(t,m(i)));
      end if;
    end loop;
  end Eval;

  procedure Eval ( c : in Complex_Vectors_of_Vectors.Vector;
                   t : in double_float; 
                   m : in Integer_Vectors_of_Vectors.Vector;
                   ctm : in out Complex_Vectors_of_Vectors.Vector ) is

  -- DESCRIPTION : ctm = c*t**m.

  begin
    for i in ctm'range loop
      Eval(c(i).all,t,m(i).all,ctm(i).all);
    end loop;
  end Eval;

  procedure Eval ( c : in Complex_Vectors_of_Vectors.Vector;
                   t : in double_float; m : in Float_Vectors_of_Vectors.Vector;
                   ctm : in out Complex_Vectors_of_Vectors.Vector ) is

  -- DESCRIPTION : ctm = c*t**m.

  begin
    for i in ctm'range loop
      Eval(c(i).all,t,m(i).all,ctm(i).all);
    end loop;
  end Eval;

 -- pragma inline(Eval);

-- HOMOTOPY CONSTRUCTOR :

  function Construct_Homotopy ( p : Laur_Sys; normal : Vector )
                              return Laur_Sys is

  -- DESCRIPTION :
  --   Given a Laurent polynomial system of dimension n*(n+1) and a
  --   normal, a homotopy will be constructed, with t = x(n+1)
  --   and so that the support of the start system corresponds with
  --   all points which give the minimal product with the normal.

    res : Laur_Sys(p'range);
    n : constant natural := p'length;
    use Complex_Multivariate_Laurent_Polynomials;
  
    function Construct_Polynomial ( p : Poly; v : Vector ) return Poly is

      res : Poly := Null_Poly;

      procedure Construct_Term ( t : in Term; cont : out boolean ) is
        rt : term;
      begin
        rt.cf := t.cf;
        rt.dg := new Integer_Vectors.Vector'(t.dg.all);
        rt.dg(n+1) := t.dg.all*v;
        Plus_Term(res,rt);
        Clear(rt);
        cont := true;
      end Construct_Term;
      procedure Construct_Terms is 
        new Complex_Multivariate_Laurent_Polynomials.Visiting_Iterator
              (Process => Construct_Term);

    begin
      Construct_Terms(p);
      return res;
    end Construct_Polynomial;

  begin
   -- SUBSTITUTIONS :
    for k in p'range loop
      res(k) := Construct_Polynomial(p(k),normal);
    end loop;
   -- SHIFT :
    for k in res'range loop
      declare
        d : integer := Minimal_Degree(res(k),n+1);
        t : Term;
      begin
        t.cf := CMPLX(1.0);
        t.dg := new Integer_Vectors.Vector'(1..n+1 => 0);
        t.dg(n+1) := -d;
        Mult_Term(res(k),t);
        Clear(t);
      end;
    end loop;
    return res;
  end Construct_Homotopy;

  function Determine_Power ( n : natural; h : Laur_Sys ) return positive is

  -- DESCRIPTION :
  --   Returns the smallest power of the last unknown,
  --   over all polynomials in h.

    res : positive := 1;
    first : boolean := true;
    d : integer;
    use Complex_Multivariate_Laurent_Polynomials;

    procedure Scan ( t : in Term; cont : out boolean ) is
    begin
      if (t.dg(n+1) > 0) and then (t.dg(n+1) < d)
       then d := t.dg(n+1);
      end if;
      cont := true;
    end Scan;
    procedure Search_Positive_Minimum is new Visiting_Iterator(Scan);

  begin
    for i in h'range loop
      d := Maximal_Degree(h(i),n+1);
      if d > 0
       then Search_Positive_Minimum(h(i));
            if first
             then res := d;
                  first := false;
             elsif d < res
                 then res := d;
            end if;
      end if;
      exit when (d=1);
    end loop;
    if res = 1
     then return res;
     else return 2;
    end if;
  end Determine_Power;

  procedure Extract_Regular ( sols : in out Solution_List ) is

    function To_Be_Removed ( flag : in natural ) return boolean is
    begin
      return ( flag /= 1 );
    end To_Be_Removed;
    procedure Extract_Regular_Solutions is new Solutions.Delete(To_Be_Removed);

  begin
    Extract_Regular_Solutions(sols);
  end Extract_Regular;

  procedure Refine ( file : in file_type; p : in Laur_Sys;
                     sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Given a polyhedral homotopy p and a list of solution for t=1,
  --   this list of solutions will be refined.

    pp : Poly_Sys(p'range) := Laurent_to_Polynomial_System(p);
    n : constant natural := p'length;
    eps : constant double_float := 10.0**(-12);
    tolsing : constant double_float := 10.0**(-8);
    max : constant natural := 3;
    numb : natural := 0;

  begin
    pp := Laurent_to_Polynomial_System(p);
    Substitute(n+1,CMPLX(1.0),pp);
   -- Reporting_Root_Refiner(file,pp,sols,eps,eps,tolsing,numb,max,false);
    Clear(pp); Extract_Regular(sols);
  end Refine;

-- FIRST LAYER OF CONTINUATION ROUTINES :

  procedure Mixed_Continuation
                 ( p : in Laur_Sys; normal : in Vector;
                   sols : in out Solution_List ) is

    n : constant natural := p'length;

    h   : Laur_Sys(p'range) := Construct_Homotopy(p,normal);
    hpe : Eval_Laur_Sys(h'range) := Create(h);
    j  : Jacobi(h'range,h'first..h'last+1) := Create(h);
    je : Eval_Jacobi(j'range(1),j'range(2)) := Create(j);

    use Complex_Multivariate_Laurent_Polynomials;

    function Eval ( x : Complex_Vectors.Vector; t : double_complex )
                  return Complex_Vectors.Vector is

      xt : Complex_Vectors.Vector(x'first..x'last+1);

    begin
      xt(x'range) := x;
      xt(xt'last) := t;
      return Eval(hpe,xt);
    end Eval;

    function dHt ( x : Complex_Vectors.Vector; t : double_complex )
                 return Complex_Vectors.Vector is

      res : Complex_Vectors.Vector(p'range);
      xt : Complex_Vectors.Vector(x'first..x'last+1);

    begin
      xt(x'range) := x;
      xt(xt'last) := t;
      for i in res'range loop
        res(i) := Eval(je(i,xt'last),xt);
      end loop;
      return res;
    end dHt;

    function dHx ( x : Complex_Vectors.Vector; t : double_complex )
                 return matrix is

      m : matrix(x'range,x'range);
      xt : Complex_Vectors.Vector(x'first..x'last+1);

    begin
      xt(x'range) := x;
      xt(xt'last) := t;
      for i in m'range(1) loop
        for j in m'range(1) loop
          m(i,j) := Eval(je(i,j),xt);
        end loop;
      end loop;
      return m;
    end dHx;

   -- pragma inline(Eval,dHt,dHx);

    procedure Laur_Cont is new Silent_Continue(Norm1,Eval,dHt,dHx);

  begin
   -- Continuation_Parameters.power_of_t := Determine_Power(h'length,h);
    Laur_Cont(sols,false);
    Clear(h); Clear(hpe); Clear(j); Clear(je);
    Extract_Regular(sols);
  end Mixed_Continuation;

  procedure Mixed_Continuation
                 ( file : in file_type; p : in Laur_Sys;
                   normal : in Vector; sols : in out Solution_List ) is

    h   : Laur_Sys(p'range) := Construct_Homotopy(p,normal);
    hpe : Eval_Laur_Sys(h'range) := Create(h);
    j  : Jacobi(h'range,h'first..h'last+1) := Create(h);
    je : Eval_Jacobi(j'range(1),j'range(2)) := Create(j);

    use Complex_Multivariate_Laurent_Polynomials;

    function Eval ( x : Complex_Vectors.Vector; t : double_complex )
                  return Complex_Vectors.Vector is

      xt : Complex_Vectors.Vector(x'first..x'last+1);

    begin
      xt(x'range) := x;
      xt(xt'last) := t;
      return Eval(hpe,xt);
    end Eval;

    function dHt ( x : Complex_Vectors.Vector; t : double_complex )
                 return Complex_Vectors.Vector is

      res : Complex_Vectors.Vector(p'range);
      xt : Complex_Vectors.Vector(x'first..x'last+1);

    begin
      xt(x'range) := x;
      xt(xt'last) := t;
      for i in res'range loop
        res(i) := Eval(je(i,xt'last),xt);
      end loop;
      return res;
    end dHt;

    function dHx ( x : Complex_Vectors.Vector; t : double_complex )
                 return matrix is

      m : matrix(x'range,x'range);
      xt : Complex_Vectors.Vector(x'first..x'last+1);

    begin
      xt(x'range) := x;
      xt(xt'last) := t;
      for i in m'range(1) loop
        for j in m'range(1) loop
           m(i,j) := Eval(je(i,j),xt);
        end loop;
      end loop;
      return m;
    end dHx;

   -- pragma inline(Eval,dHt,dHx);

    procedure Laur_Cont is new Reporting_Continue(Norm1,Eval,dHt,dHx);

  begin
   -- Continuation_Parameters.power_of_t := Determine_Power(h'length,h);
    Laur_Cont(file,sols,false);
    Clear(h); Clear(hpe); Clear(j); Clear(je);
    Extract_Regular(sols);
  end Mixed_Continuation;

  procedure Mixed_Continuation
                ( mix : in Integer_Vectors.Vector;
                  lifted : in Array_of_Lists; h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  normal : in Vector; sols : in out Solution_List ) is

    pow : Integer_Vectors_of_Vectors.Vector(c'range)
        := Create(e,lifted,mix,normal);
   -- scapow : Float_Vectors_of_Vectors.Vector(c'range) := Scale(pow);
    ctm : Complex_Vectors_of_Vectors.Vector(c'range);

    use Complex_Multivariate_Laurent_Polynomials;

    function Eval ( x : Complex_Vectors.Vector; t : double_complex )
                  return Complex_Vectors.Vector is
    begin
     -- Eval(c,REAL_PART(t),scapow,ctm);
      Eval(c,REAL_PART(t),pow,ctm);
      return Eval(h,ctm,x);
    end Eval;

    function dHt ( x : Complex_Vectors.Vector; t : double_complex )
                 return Complex_Vectors.Vector is

      res : Complex_Vectors.Vector(h'range);
      xtl : constant integer := x'last+1;

    begin
     -- Eval(c,REAL_PART(t),scapow,ctm);
      Eval(c,REAL_PART(t),pow,ctm);
      for i in res'range loop
        res(i) := Eval(j(i,xtl),m(i,xtl).all,ctm(i).all,x);
      end loop;
      return res;
    end dHt;

    function dHx ( x : Complex_Vectors.Vector; t : double_complex )
                 return matrix is

      mt : Matrix(x'range,x'range);

    begin
     -- Eval(c,REAL_PART(t),scapow,ctm);
      Eval(c,REAL_PART(t),pow,ctm);
      for k in mt'range(1) loop
        for l in mt'range(2) loop
          mt(k,l) := Eval(j(k,l),m(k,l).all,ctm(k).all,x);
        end loop;
      end loop;
      return mt;
    end dHx;

   -- pragma inline(Eval,dHt,dHx);

    procedure Laur_Cont is new Silent_Continue(Norm1,Eval,dHt,dHx);

  begin
    for i in c'range loop
      ctm(i) := new Complex_Vectors.Vector(c(i)'range);
    end loop;
    Laur_Cont(sols,false);
    Integer_Vectors_of_Vectors.Clear(pow);
   -- Float_Vectors_of_Vectors.Clear(scapow);
    Complex_Vectors_of_Vectors.Clear(ctm);
    Extract_Regular(sols);
  end Mixed_Continuation;

  procedure Mixed_Continuation
                ( file : in file_type; mix : in Integer_Vectors.Vector;
                  lifted : in Array_of_Lists; h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  normal : in Vector; sols : in out Solution_List ) is

    pow : Integer_Vectors_of_Vectors.Vector(c'range) 
        := Create(e,lifted,mix,normal);
   -- scapow : Float_Vectors_of_Vectors.Vector(c'range) := Scale(pow);
    ctm : Complex_Vectors_of_Vectors.Vector(c'range);

    use Complex_Multivariate_Laurent_Polynomials;

    function Eval ( x : Complex_Vectors.Vector; t : double_complex )
                  return Complex_Vectors.Vector is
    begin
     -- Eval(c,REAL_PART(t),scapow,ctm);
      Eval(c,REAL_PART(t),pow,ctm);
      return Eval(h,ctm,x);
    end Eval;

    function dHt ( x : Complex_Vectors.Vector; t : double_complex )
                 return Complex_Vectors.Vector is

      res : Complex_Vectors.Vector(h'range);
      xtl : constant integer := x'last+1;

    begin
     -- Eval(c,REAL_PART(t),scapow,ctm);
      Eval(c,REAL_PART(t),pow,ctm);
      for i in res'range loop
        res(i) := Eval(j(i,xtl),m(i,xtl).all,ctm(i).all,x);
      end loop;
      return res;
    end dHt;

    function dHx ( x : Complex_Vectors.Vector; t : double_complex )
                 return matrix is

      mt : Matrix(x'range,x'range);

    begin
     -- Eval(c,REAL_PART(t),scapow,ctm);
      Eval(c,REAL_PART(t),pow,ctm);
      for k in m'range(1) loop
        for l in m'range(1) loop
          mt(k,l) := Eval(j(k,l),m(k,l).all,ctm(k).all,x);
        end loop;
      end loop;
      return mt;
    end dHx;

   -- pragma inline(Eval,dHt,dHx);

    procedure Laur_Cont is new Reporting_Continue(Norm1,Eval,dHt,dHx);

  begin
   -- put(file,"The normal : "); put(file,normal); new_line(file);
   -- put_line(file,"The exponent vector : ");
   -- for i in pow'range loop
   --   put(file,pow(i)); new_line(file);
   -- end loop;
   -- put_line(file,"The coefficient vector : "); Write(file,c);
    for i in c'range loop
      ctm(i) := new Complex_Vectors.Vector(c(i)'range);
    end loop;
    Laur_Cont(file,sols,false);
    Integer_Vectors_of_Vectors.Clear(pow);
   -- Float_Vectors_of_Vectors.Clear(scapow);
    Complex_Vectors_of_Vectors.Clear(ctm);
    Extract_Regular(sols);
  end Mixed_Continuation;

-- UTILITIES FOR SECOND LAYER :

  function Sub_Lifting ( q : Laur_Sys; mix : Vector; mic : Mixed_Cell )
                       return Array_of_Lists is

  -- DESCRIPTION :
  --   Returns the lifting used to subdivide the cell.

    res : Array_of_Lists(mix'range);
    sup : Array_of_Lists(q'range);
    n : constant natural := q'last;
    cnt : natural := sup'first;

  begin
    for i in mic.pts'range loop
      sup(cnt) := Reduce(mic.pts(i),q'last+1);
      for j in 1..(mix(i)-1) loop
        Copy(sup(cnt),sup(cnt+j));
      end loop;
      cnt := cnt + mix(i);
    end loop;
    res := Induced_Lifting(n,mix,sup,mic.sub.all);
    Deep_Clear(sup);
    return res;
  end Sub_Lifting;

  procedure Refined_Mixed_Solve
               ( q : in Laur_Sys; mix : in Vector; mic : in Mixed_Cell;
                 qsols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Solves a subsystem q using the refinement of the cell mic.

  -- REQUIRED : mic.sub /= null.

    lif : Array_of_Lists(mix'range) := Sub_Lifting(q,mix,mic);
    lifq : Laur_Sys(q'range) := Perform_Lifting(q'last,mix,lif,q);

  begin
    Mixed_Solve(lifq,mix,mic.sub.all,qsols);
    Deep_Clear(lif); Clear(lifq);
  end Refined_Mixed_Solve;

  procedure Refined_Mixed_Solve
               ( file : in file_type;
                 q : in Laur_Sys; mix : in Vector; mic : in Mixed_Cell;
                 qsols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Solves a subsystem q using the refinement of the cell mic.

  -- REQUIRED : mic.sub /= null.

    lif : Array_of_Lists(mix'range) := Sub_Lifting(q,mix,mic);
    lifq : Laur_Sys(q'range) := Perform_Lifting(q'last,mix,lif,q);

  begin
    Mixed_Solve(file,lifq,mix,mic.sub.all,qsols);
    Deep_Clear(lif); Clear(lifq);
  end Refined_Mixed_Solve;

  function Sub_Polyhedral_Homotopy
               ( l : List; e : Integer_Vectors_of_Vectors.Vector;
                 c : Complex_Vectors.Vector ) return Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   For every vector in e that does not belong to l, the corresponding
  --   index in c will be set to zero, otherwise it is copied to the result.

    res : Complex_Vectors.Vector(c'range);
    found : boolean;
    lif : integer;

  begin
    for i in e'range loop
      Search_Lifting(l,e(i).all,found,lif);
      if not found
       then res(i) := CMPLX(0.0);
       else res(i) := c(i);
      end if;
    end loop;
    return res;
  end Sub_Polyhedral_Homotopy;

  function Sub_Polyhedral_Homotopy
               ( mix : Vector; mic : Mixed_Cell;
                 e : Exponent_Vectors_Array;
                 c : Complex_Vectors_of_Vectors.Vector )
               return Complex_Vectors_of_Vectors.Vector is

  -- DESCRIPTION :
  --   Given a subsystem q of p and the coefficient vector of p, the
  --   vector on return will have only nonzero entries for coefficients
  --   that belong to q.

    res : Complex_Vectors_of_Vectors.Vector(c'range);

  begin
    for i in mix'range loop
      declare
        cri : constant Complex_Vectors.Vector
            := Sub_Polyhedral_Homotopy(mic.pts(i),e(i).all,c(i).all);
      begin
        res(i) := new Complex_Vectors.Vector'(cri);
        for j in 1..(mix(i)-1) loop
          declare
            crj : constant Complex_Vectors.Vector
                := Sub_Polyhedral_Homotopy(mic.pts(i),e(i+j).all,c(i+j).all);
          begin
            res(i+j) := new Complex_Vectors.Vector'(crj);
          end;
        end loop;
      end;
    end loop;
    return res;
  end Sub_Polyhedral_Homotopy;

  procedure Refined_Mixed_Solve
                ( q : in Laur_Sys; mix : in Vector; mic : in Mixed_Cell;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  qsols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Polyhedral coeffient-homotopy for subsystem q.

  -- REQUIRED : mic.sub /= null.

    lif : Array_of_Lists(mix'range) := Sub_Lifting(q,mix,mic);
    lq : Laur_Sys(q'range) := Perform_Lifting(q'last,mix,lif,q);
    cq : Complex_Vectors_of_Vectors.Vector(c'range)
       := Sub_Polyhedral_Homotopy(mix,mic,e,c);

  begin
    Mixed_Solve(lq,lif,h,cq,e,j,m,mix,mic.sub.all,qsols);
    Complex_Vectors_of_Vectors.Clear(cq); Deep_Clear(lif); Clear(lq);
  end Refined_Mixed_Solve;

  procedure Refined_Mixed_Solve
                ( file : in file_type; q : in Laur_Sys; mix : in Vector;
                  mic : in Mixed_Cell; h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  qsols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Polyhedral coeffient-homotopy for subsystem q.

  -- REQUIRED : mic.sub /= null.

    lif : Array_of_Lists(mix'range) := Sub_Lifting(q,mix,mic);
    lq : Laur_Sys(q'range) := Perform_Lifting(q'last,mix,lif,q);
    cq : Complex_Vectors_of_Vectors.Vector(c'range)
       := Sub_Polyhedral_Homotopy(mix,mic,e,c);

  begin
    Mixed_Solve(file,lq,lif,h,cq,e,j,m,mix,mic.sub.all,qsols);
    Complex_Vectors_of_Vectors.Clear(cq); Deep_Clear(lif); Clear(lq);
  end Refined_Mixed_Solve;

-- TARGET ROUTINES FOR SECOND LAYER :

  procedure Mixed_Solve 
               ( p : in Laur_Sys; mix : in Vector; mic : in Mixed_Cell;
                 sols,sols_last : in out Solution_List ) is

    q : Laur_Sys(p'range) := Select_Subsystem(p,mix,mic);
    sq : Laur_Sys(q'range);
    pq : Poly_Sys(q'range);
    qsols : Solution_List;
    len : natural := 0;
    fail : boolean;

  begin
    Reduce(q'last+1,q); sq := Shift(q);
    Fewnomials.Solve(sq,qsols,fail);
    if fail
     then if mic.sub = null
           then pq := Laurent_to_Polynomial_System(sq);
                qsols := Solve_by_Static_Lifting(pq);
                Clear(pq);
           else Refined_Mixed_Solve(q,mix,mic,qsols);
          end if;
          Set_Continuation_Parameter(qsols,CMPLX(0.0));
    end if;
    len := Length_Of(qsols);
    if len > 0
     then Mixed_Continuation(p,mic.nor.all,qsols);
          Concat(sols,sols_last,qsols);
    end if;
    Clear(sq); Clear(q); Shallow_Clear(qsols);
  end Mixed_Solve;

  procedure Mixed_Solve 
               ( file : in file_type; p : in Laur_Sys;
                 mix : in Vector; mic : in Mixed_Cell;
                 sols,sols_last : in out Solution_List ) is

    q : Laur_Sys(p'range) := Select_Subsystem(p,mix,mic);
    sq : Laur_Sys(q'range);
    pq : Poly_Sys(q'range);
    qsols : Solution_List;
    len : natural := 0;
    fail : boolean;

  begin
    Reduce(q'last+1,q); sq := Shift(q);
    Fewnomials.Solve(sq,qsols,fail);
    if not fail
     then put_line(file,"It is a fewnomial system.");
     else put_line(file,"No fewnomial system.");
          if mic.sub = null
           then put_line(file,"Calling the black box solver.");
                pq := Laurent_to_Polynomial_System(sq);
                qsols := Solve_by_Static_Lifting(file,pq);
                Clear(pq);
           else put_line(file,"Using the refinement of the cell.");
                Refined_Mixed_Solve(file,q,mix,mic,qsols);
          end if;
          Set_Continuation_Parameter(qsols,CMPLX(0.0));
    end if;
    len := Length_Of(qsols);
    put(file,len,1); put_line(file," solutions found.");
    if len > 0
     then Mixed_Continuation(file,p,mic.nor.all,qsols);
          Concat(sols,sols_last,qsols);
    end if;
    Clear(sq); Clear(q); Shallow_Clear(qsols);
  end Mixed_Solve;

  procedure Mixed_Solve
                ( p : in Laur_Sys; lifted : in Array_of_Lists;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  mix : in Vector; mic : in Mixed_Cell;
                  sols,sols_last : in out Solution_List ) is

    q : Laur_Sys(p'range) := Select_Subsystem(p,mix,mic);
    sq : Laur_Sys(q'range);
    pq : Poly_Sys(q'range);
    qsols : Solution_List;
    len : natural := 0;
    fail : boolean;

  begin
    Reduce(q'last+1,q); sq := Shift(q);
    Fewnomials.Solve(sq,qsols,fail);
    if fail
     then if mic.sub = null
           then pq := Laurent_to_Polynomial_System(sq);
                qsols := Solve_by_Static_Lifting(pq);
                Clear(pq);
           else -- Refined_Mixed_Solve(q,mix,mic,qsols);
                Refined_Mixed_Solve(q,mix,mic,h,c,e,j,m,qsols);
          end if;
          Set_Continuation_Parameter(qsols,CMPLX(0.0));
    end if;
    len := Length_Of(qsols);
    if len > 0
     then Mixed_Continuation(mix,lifted,h,c,e,j,m,mic.nor.all,qsols);
          Concat(sols,sols_last,qsols);
    end if;
    Clear(sq); Clear(q); Shallow_Clear(qsols);
  end Mixed_Solve;

  procedure Mixed_Solve
                ( file : in file_type;
                  p : in Laur_Sys; lifted : in Array_of_Lists;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  mix : in Vector; mic : in Mixed_Cell;
                  sols,sols_last : in out Solution_List ) is

    q : Laur_Sys(p'range) := Select_Subsystem(p,mix,mic);
    sq : Laur_Sys(q'range);
    pq : Poly_Sys(q'range);
    qsols : Solution_List;
    len : natural := 0;
    fail : boolean;

  begin
    Reduce(q'last+1,q); sq := Shift(q);
    Fewnomials.Solve(sq,qsols,fail);
    if not fail
     then put_line(file,"It is a fewnomial system.");
     else put_line(file,"No fewnomial system.");
          if mic.sub = null
           then put_line(file,"Calling the black box solver.");
                pq := Laurent_to_Polynomial_System(sq);
                qsols := Solve_by_Static_Lifting(file,pq);
                Clear(pq);
           else put_line(file,"Using the refinement of the cell.");
               -- Refined_Mixed_Solve(file,q,mix,mic,qsols);
                Refined_Mixed_Solve(file,q,mix,mic,h,c,e,j,m,qsols);
          end if;
          Set_Continuation_Parameter(qsols,CMPLX(0.0));
    end if;
    len := Length_Of(qsols);
    put(file,len,1); put_line(file," solutions found.");
    if len > 0
     then Mixed_Continuation(file,mix,lifted,h,c,e,j,m,mic.nor.all,qsols);
          Concat(sols,sols_last,qsols);
    end if;
    Clear(sq); Clear(q); Shallow_Clear(qsols);
  end Mixed_Solve;

-- THIRD LAYER :

  procedure Mixed_Solve
               ( p : in Laur_Sys;
                 mix : in Vector; mixsub : in Mixed_Subdivision;
                 sols : in out Solution_List ) is

    tmp : Mixed_Subdivision := mixsub;
    mic : Mixed_Cell;
    sols_last : Solution_List;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      Mixed_Solve(p,mix,mic,sols,sols_last);
      tmp := Tail_Of(tmp);
    end loop;
  end Mixed_Solve;

  procedure Mixed_Solve 
               ( file : in file_type; p : in Laur_Sys;
                 mix : in Vector; mixsub : in Mixed_Subdivision;
                 sols : in out Solution_List ) is

    tmp : Mixed_Subdivision := mixsub;
    mic : Mixed_Cell;
    sols_last : Solution_List;
    cnt : natural := 0;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      cnt := cnt + 1;
      new_line(file);
      put(file,"*** PROCESSING SUBSYSTEM "); put(file,cnt,1);
      put_line(file," ***");
      new_line(file);
      Mixed_Solve(file,p,mix,mic,sols,sols_last);
      tmp := Tail_Of(tmp);
    end loop;
  end Mixed_Solve;

  procedure Mixed_Solve
                ( p : in Laur_Sys; lifted : in Array_of_Lists;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  mix : in Vector; mixsub : in Mixed_Subdivision;
                  sols : in out Solution_List ) is

    tmp : Mixed_Subdivision := mixsub;
    mic : Mixed_Cell;
    sols_last : Solution_List;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      Mixed_Solve(p,lifted,h,c,e,j,m,mix,mic,sols,sols_last);
      tmp := Tail_Of(tmp);
    end loop;
  end Mixed_Solve;

  procedure Mixed_Solve
                ( file : in file_type;
                  p : in Laur_Sys; lifted : in Array_of_Lists;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  mix : in Vector; mixsub : in Mixed_Subdivision;
                  sols : in out Solution_List ) is

    tmp : Mixed_Subdivision := mixsub;
    mic : Mixed_Cell;
    sols_last : Solution_List;
    cnt : natural := 0;

  begin
    while not Is_Null(tmp) loop
      mic := Head_Of(tmp);
      cnt := cnt + 1;
      new_line(file);
      put(file,"*** PROCESSING SUBSYSTEM "); put(file,cnt,1);
      put_line(file," ***");
      new_line(file);
      Mixed_Solve(file,p,lifted,h,c,e,j,m,mix,mic,sols,sols_last);
      tmp := Tail_Of(tmp);
    end loop;
  end Mixed_Solve;

end Integer_Polyhedral_Continuation;
