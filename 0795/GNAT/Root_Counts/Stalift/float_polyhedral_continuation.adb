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

with Float_Integer_Convertors;           use Float_Integer_Convertors;
with Lists_of_Float_Vectors;             use Lists_of_Float_Vectors;
with Float_Lifting_Utilities;            use Float_Lifting_Utilities;

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

package body Float_Polyhedral_Continuation is

-- UTILITIES FOR POLYHEDRAL COEFFICIENT HOMOTOPY :

  function Power ( x,m : double_float ) return double_float is

  -- DESCRIPTION :
  --   Computes x**m for high powers of m to avoid RANGE_ERROR.

   -- intm : integer := integer(m);
   -- fltm : double_float;

  begin
   -- if m < 10.0
   --  then
    return (x**m);  -- GNAT is better at this
   --  else if double_float(intm) > m
   --        then intm := intm-1;
   --       end if;
   --       fltm := m - double_float(intm);
   --       return ((x**intm)*(x**fltm));
   -- end if;
  end Power;

  function Minimum ( v : Float_Vectors.Vector ) return double_float is

  -- DESCRIPTION :
  --   Returns the minimal element (/= 0)  in the vector v.

    tol : constant double_float := 10.0**(-7);
    min : double_float := 0.0;
    tmp : double_float;

  begin
    for i in v'range loop
      tmp := abs(v(i));
      if tmp > tol
       then if (min = 0.0) or else (tmp < min)
             then min := tmp;
            end if;
      end if;
    end loop;
    return min;
  end Minimum;

  function Scale ( v : Float_Vectors.Vector ) return Float_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the scaled vector such that the minimal element (/= 0)
  --   equals one.

    res : Float_Vectors.Vector(v'range);
    min : constant double_float := Minimum(v);

  begin
    if (min /= 0.0) and (min /= 1.0)
     then for i in res'range loop
            res(i) := double_float(v(i))/double_float(min);
          end loop;
     else res := v;
    end if;
    return res;
  end Scale;

  function Scale ( v : Float_Vectors_of_Vectors.Vector )
                 return Float_Vectors_of_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns an array of scaled vectors.

    res : Float_Vectors_of_Vectors.Vector(v'range);

  begin
    for i in v'range loop
      declare
        sv : constant Float_Vectors.Vector := Scale(v(i).all);
      begin
        res(i) := new Float_Vectors.Vector'(sv);
      end;
    end loop;
    return res;
  end Scale;

  procedure Shift ( v : in out Float_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Shifts the elements in v, such that the minimal element equals zero.

    min : double_float := v(v'first);

  begin
    for i in v'first+1..v'last loop
      if v(i) < min
       then min := v(i);
      end if;
    end loop;
    if min /= 0.0
     then for i in v'range loop
            v(i) := v(i) - min;
          end loop;
    end if;
  end Shift;

  function Allocate ( c : Complex_Vectors_of_Vectors.Vector )
                    return Complex_Vectors_of_Vectors.Vector is

    res : Complex_Vectors_of_Vectors.Vector(c'range);

  begin
    for i in res'range loop
      res(i) := new Complex_Vectors.Vector'(c(i)'range => CMPLX(0.0));
    end loop;
    return res;
  end Allocate;

  function Allocate ( n : natural; mix : Integer_Vectors.Vector;
                      l : Array_of_Lists )
                    return Float_Vectors_of_Vectors.Vector is

    res : Float_Vectors_of_Vectors.Vector(1..n);
    cnt : natural := 1;
 
  begin
    for i in mix'range loop
      res(cnt) := new Float_Vectors.Vector'(1..Length_Of(l(i)) => 0.0);
      for j in 1..(mix(i)-1) loop
        res(cnt+j) := new Float_Vectors.Vector'(res(cnt).all);
      end loop;
      cnt := cnt + mix(i);
    end loop;
    return res;
  end Allocate;

  function Create ( l : List; normal : Float_Vectors.Vector ) 
                  return Float_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns a vector with all inner products of the normal with
  --   the exponents in the list, such that minimal value equals zero.

    res : Float_Vectors.Vector(1..Length_Of(l));
    tmp : List := l;
    use Float_Vectors;

  begin
    for i in res'range loop
      res(i) := Head_Of(tmp).all*normal;
      tmp := Tail_Of(tmp);
    end loop;
    Shift(res);
    return res;
  end Create;

  function Create ( l : Array_of_Lists; mix : Integer_Vectors.Vector;
                    normal : Float_Vectors.Vector )
                  return Float_Vectors_of_Vectors.Vector is

    res : Float_Vectors_of_Vectors.Vector(normal'first..normal'last-1);
    cnt : natural := res'first; 

  begin
    for i in mix'range loop
      res(cnt) := new Float_Vectors.Vector'(Create(l(i),normal));
      for j in 1..(mix(i)-1) loop
        res(cnt+j) := new Float_Vectors.Vector'(res(cnt).all);
      end loop;
      cnt := cnt + mix(i);
    end loop;
    return res;
  end Create;

  procedure Eval ( c : in Complex_Vectors.Vector;
                   t : in double_float; m : in Float_Vectors.Vector;
                   ctm : in out Complex_Vectors.Vector ) is

  -- DESCRIPTION : ctm = c*t**m.

  begin
    for i in ctm'range loop
      -- ctm(i) := c(i)*CMPLX((t**m(i)));
      ctm(i) := c(i)*CMPLX(Power(t,m(i)));
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

-- USEFUL PRIMITIVES :

  function Is_In ( l : List;
                   d : Complex_Multivariate_Laurent_Polynomials.Degrees )
                 return boolean is

  -- DESCRIPTION :
  --   Returns true if the degrees belong to the list l.

    tmp : List := l;
    pt : Float_Vectors.Link_to_Vector;
    fld : Float_Vectors.Vector(d'range);
    equ : boolean;

  begin
    for i in fld'range loop
      fld(i) := double_float(d(i));
    end loop;
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      equ := true;
      for i in fld'range loop
        if pt(i) /= fld(i)
         then equ := false;
        end if;
        exit when not equ;
      end loop;
      if equ
       then return true;
       else tmp := Tail_Of(tmp);
      end if;
    end loop;
    return false;
  end Is_In;

  function Select_Terms ( p : Complex_Multivariate_Laurent_Polynomials.Poly;
                          l : List )
                        return Complex_Multivariate_Laurent_Polynomials.Poly is

    use Complex_Multivariate_Laurent_Polynomials;
    res : Poly := Null_Poly;

    procedure Visit_Term ( t : in Term; cont : out boolean ) is
    begin
      if Is_In(l,t.dg)
       then Plus_Term(res,t);
      end if;
      cont := true;
    end Visit_Term;
    procedure Visit_Terms is new Visiting_Iterator(Visit_Term);

  begin
    Visit_Terms(p);
    return res;
  end Select_Terms;

  function Select_Subsystem
               ( p : Laur_Sys; mix : Integer_Vectors.Vector; mic : Mixed_Cell )
               return Laur_Sys is

  -- DESCRIPTION :
  --   Given a Laurent polynomial system and a mixed cell,
  --   the corresponding subsystem will be returned.

  -- ON ENTRY :
  --   p          a Laurent polynomial system;
  --   mix        type of mixture: occurencies of the supports;
  --   mic        a mixed cell.

  -- REQUIRED :
  --   The polynomials in p must be ordered according to the type of mixture.

    res : Laur_Sys(p'range);
    cnt : natural := 0;

  begin
    for k in mix'range loop
      for l in 1..mix(k) loop
        cnt := cnt + 1;
        res(cnt) := Select_Terms(p(cnt),mic.pts(k));
      end loop;
    end loop;
    return res;
  end Select_Subsystem;

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
                ( mix : in Integer_Vectors.Vector;
                  lifted : in Array_of_Lists; h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  normal : in Float_Vectors.Vector;
                  sols : in out Solution_List ) is

    pow : Float_Vectors_of_Vectors.Vector(c'range)
        := Create(lifted,mix,normal);
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

    procedure Laur_Cont is new Silent_Continue(Norm1,Eval,dHt,dHx);

  begin
    for i in ctm'range loop
      ctm(i) := new Complex_Vectors.Vector'(c(i).all);
    end loop;
    Laur_Cont(sols,false);
    Float_Vectors_of_Vectors.Clear(pow);
   -- Float_Vectors_of_Vectors.Clear(scapow);
    Complex_Vectors_of_Vectors.Clear(ctm);
    Extract_Regular(sols);
  end Mixed_Continuation;

  procedure Mixed_Continuation
                ( file : in file_type; mix : in Integer_Vectors.Vector;
                  lifted : in Array_of_Lists; h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  normal : in Float_Vectors.Vector;
                  sols : in out Solution_List ) is

    pow : Float_Vectors_of_Vectors.Vector(c'range);
      --  := Allocate(h'length,mix,lifted);
      --  := Create(lifted,mix,normal);
   -- scapow : Float_Vectors_of_Vectors.Vector(c'range) := Scale(pow);
    ctm : Complex_Vectors_of_Vectors.Vector(c'range); -- := Allocate(c);

    cnt : natural;
    tmp : List;

    use Float_Vectors;

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

    procedure Laur_Cont is new Reporting_Continue(Norm1,Eval,dHt,dHx);

  begin
   -- put_line(file,"The coefficient vectors :" );
   -- for i in c'range loop
   --   put(file,c(i),3,3,3); new_line(file);
   -- end loop;
    cnt := pow'first;
    for i in mix'range loop
      tmp := lifted(i);
      pow(cnt) := new Float_Vectors.Vector(1..Length_Of(lifted(i)));
      for jj in pow(cnt)'range loop
        pow(cnt)(jj) := Head_Of(tmp).all*normal;
        tmp := Tail_Of(tmp);
      end loop;
      Shift(pow(cnt).all);
      for k in 1..(mix(i)-1) loop
        pow(cnt+k) := new Float_Vectors.Vector(pow(cnt)'range);
        for jj in pow(cnt)'range loop
          pow(cnt+k)(jj) := pow(cnt)(jj);
        end loop;
      end loop;
      cnt := cnt + mix(i);
    end loop;

    for i in c'range loop
      ctm(i) := new Complex_Vectors.Vector'(c(i).all'range => CMPLX(0.0));
    end loop;

   -- put(file,"The normal : "); put(file,normal,3,3,3); new_line(file);
   -- put_line(file,"The exponent vector : ");
   -- for i in pow'range loop
   --   put(file,pow(i),3,3,3); new_line(file);
   -- end loop;
   -- put_line(file,"The scaled exponent vector : ");
   -- for i in pow'range loop
   --   put(file,scapow(i),3,3,3); new_line(file);
   -- end loop;
    Laur_Cont(file,sols,false);
    Float_Vectors_of_Vectors.Clear(pow);
   -- Float_Vectors_of_Vectors.Clear(scapow);
    Complex_Vectors_of_Vectors.Clear(ctm);
    Extract_Regular(sols);
  end Mixed_Continuation;

-- UTILITIES FOR SECOND LAYER :

  function Remove_Lifting ( l : List ) return List is

  -- DESCRIPTION :
  --   Removes the lifting value from the vectors.

    tmp,res,res_last : List;

  begin
    tmp := l;
    while not Is_Null(tmp) loop
      declare
        d1 : constant Float_Vectors.Vector := Head_Of(tmp).all;
        d2 : constant Float_Vectors.Vector := d1(d1'first..d1'last-1);
      begin
        Append(res,res_last,d2);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    return res;
  end Remove_Lifting;

  function Sub_Lifting ( q : Laur_Sys; mix : Integer_Vectors.Vector;
                         mic : Mixed_Cell ) return Array_of_Lists is

  -- DESCRIPTION :
  --   Returns the lifting used to subdivide the cell.

    res : Array_of_Lists(mix'range);
    sup : Array_of_Lists(q'range);
    n : constant natural := q'last;
    cnt : natural := sup'first;

  begin
    for i in mic.pts'range loop
      sup(cnt) := Remove_Lifting(mic.pts(i));
      for j in 1..(mix(i)-1) loop
        Copy(sup(cnt),sup(cnt+j));
      end loop;
      cnt := cnt + mix(i);
    end loop;
    res := Induced_Lifting(n,mix,sup,mic.sub.all);
    Deep_Clear(sup);
    return res;
  end Sub_Lifting;

  function Sub_Polyhedral_Homotopy
               ( l : List; e : Integer_Vectors_of_Vectors.Vector;
                 c : Complex_Vectors.Vector ) return Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   For every vector in e that does not belong to l, the corresponding
  --   index in c will be set to zero, otherwise it is copied to the result.

    res : Complex_Vectors.Vector(c'range);
    found : boolean;
    lif : double_float;

  begin
    for i in e'range loop
      declare
        fei : constant Float_Vectors.Vector := Convert(e(i).all);
      begin
        Search_Lifting(l,fei,found,lif);
        if not found
         then res(i) := CMPLX(0.0);
         else res(i) := c(i);
        end if;
      end;
    end loop;
    return res;
  end Sub_Polyhedral_Homotopy;

  function Sub_Polyhedral_Homotopy
               ( mix : Integer_Vectors.Vector; mic : Mixed_Cell;
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
                ( q : in Laur_Sys; mix : in Integer_Vectors.Vector;
                  mic : in Mixed_Cell; h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  qsols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Polyhedral coeffient-homotopy for subsystem q.

  -- REQUIRED : mic.sub /= null.

    lif : Array_of_Lists(mix'range) := Sub_Lifting(q,mix,mic);
    cq : Complex_Vectors_of_Vectors.Vector(c'range)
       := Sub_Polyhedral_Homotopy(mix,mic,e,c);

  begin
    Mixed_Solve(q,lif,h,cq,e,j,m,mix,mic.sub.all,qsols);
    Complex_Vectors_of_Vectors.Clear(cq); Deep_Clear(lif);
  end Refined_Mixed_Solve;

  procedure Refined_Mixed_Solve
                ( file : in file_type; q : in Laur_Sys;
                  mix : in Integer_Vectors.Vector;
                  mic : in Mixed_Cell; h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  qsols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Polyhedral coeffient-homotopy for subsystem q.

  -- REQUIRED : mic.sub /= null.

    lif : Array_of_Lists(mix'range) := Sub_Lifting(q,mix,mic);
    cq : Complex_Vectors_of_Vectors.Vector(c'range)
       := Sub_Polyhedral_Homotopy(mix,mic,e,c);

  begin
    Mixed_Solve(file,q,lif,h,cq,e,j,m,mix,mic.sub.all,qsols);
    Complex_Vectors_of_Vectors.Clear(cq); Deep_Clear(lif);
  end Refined_Mixed_Solve;

-- SECOND LAYER :

  procedure Mixed_Solve
                ( p : in Laur_Sys; lifted : in Array_of_Lists;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  mix : in Integer_Vectors.Vector; mic : in Mixed_Cell;
                  sols,sols_last : in out Solution_List ) is

    q : Laur_Sys(p'range) := Select_Subsystem(p,mix,mic);
    sq : Laur_Sys(q'range) := Shift(q);
    pq : Poly_Sys(q'range);
    qsols : Solution_List;
    len : natural := 0;
    fail : boolean;

  begin
    Fewnomials.Solve(sq,qsols,fail);
    if fail
     then if mic.sub = null
           then pq := Laurent_to_Polynomial_System(sq);
                qsols := Solve_by_Static_Lifting(pq);
                Clear(pq);
           else Refined_Mixed_Solve(q,mix,mic,h,c,e,j,m,qsols);
          end if;
          Set_Continuation_Parameter(qsols,CMPLX(0.0));
    end if;
    len := Length_Of(qsols);
    if len > 0
     then Mixed_Continuation(mix,lifted,h,c,j,m,mic.nor.all,qsols);
          Concat(sols,sols_last,qsols);
    end if;
    Clear(sq); Clear(q); Clear(qsols);
  end Mixed_Solve;

  procedure Mixed_Solve
                ( file : in file_type;
                  p : in Laur_Sys; lifted : in Array_of_Lists;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  mix : in Integer_Vectors.Vector; mic : in Mixed_Cell;
                  sols,sols_last : in out Solution_List ) is

    q : Laur_Sys(p'range) := Select_Subsystem(p,mix,mic);
    sq : Laur_Sys(q'range) := Shift(q);
    pq : Poly_Sys(q'range);
    qsols : Solution_List;
    len : natural := 0;
    fail : boolean;

  begin
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
                Refined_Mixed_Solve(file,q,mix,mic,h,c,e,j,m,qsols);
          end if;
          Set_Continuation_Parameter(qsols,CMPLX(0.0));
    end if;
    len := Length_Of(qsols);
    put(file,len,1); put_line(file," solutions found.");
    if len > 0
     then 
       Mixed_Continuation(file,mix,lifted,h,c,j,m,mic.nor.all,qsols);
       Concat(sols,sols_last,qsols);
    end if;
    Clear(sq); Clear(q); Clear(qsols);
  end Mixed_Solve;

-- THIRD LAYER :

  procedure Mixed_Solve
                ( p : in Laur_Sys; lifted : in Array_of_Lists;
                  h : in Eval_Coeff_Laur_Sys;
                  c : in Complex_Vectors_of_Vectors.Vector;
                  e : in Exponent_Vectors_Array;
                  j : in Eval_Coeff_Jacobi; m : in Mult_Factors;
                  mix : in Integer_Vectors.Vector;
                  mixsub : in Mixed_Subdivision;
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
                  mix : in Integer_Vectors.Vector;
                  mixsub : in Mixed_Subdivision;
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

end Float_Polyhedral_Continuation;
