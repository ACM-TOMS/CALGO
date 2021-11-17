with Integer_Vectors;
with Integer_Vectors_io;                 use Integer_Vectors_io;
with Integer_Vectors_Utilities;          use Integer_Vectors_Utilities;
with Integer_Vectors_of_Vectors;
with Integer_Linear_System_Solvers;      use Integer_Linear_System_Solvers;
with Floating_Point_Numbers;             use Floating_Point_Numbers;
with Complex_Numbers,Complex_Vectors;    use Complex_Numbers;
with Complex_Norms,Complex_Matrices;     use Complex_Norms;
with Random_Number_Generators;

with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors;
with Lists_of_Vectors_Utilities;         use Lists_of_Vectors_Utilities;
with Integer_Support_Functions;          use Integer_Support_Functions;
with Arrays_of_Integer_Vector_Lists;     use Arrays_of_Integer_Vector_Lists;
with Arrays_of_Lists_Utilities;          use Arrays_of_Lists_Utilities;
with Volumes;

with Transformations;                    use Transformations;
with Transforming_Solutions;             use Transforming_Solutions;
with Transforming_Laurent_Systems;       use Transforming_Laurent_Systems;
with Transforming_Integer_Vector_Lists;  use Transforming_Integer_Vector_Lists;

with Power_Lists;                        use Power_Lists;
with Complex_Polynomial_Systems;         use Complex_Polynomial_Systems;
with Complex_Multivariate_Polynomials;
with Complex_Multivariate_Laurent_Polynomials;
with Polynomial_to_Laurent_Converters;   use Polynomial_to_Laurent_Converters;
with Laurent_to_Polynomial_Converters;   use Laurent_to_Polynomial_Converters;

with Durand_Kerner;
with Jacobi_Matrices;                   use Jacobi_Matrices;
with Increment_and_Fix_Continuation;    use Increment_and_Fix_Continuation;
with Root_Refiners;                     use Root_Refiners;

package body Mixed_Homotopy_Continuation is

-- INVARIANT CONDITION :
--   The procedures and functions in this package `mirror' corresponding
--   routines in the package Volumes.

  use Complex_Multivariate_Laurent_Polynomials;

-- AUXILIARIES :

  type bar is array ( integer range <> ) of boolean;

  procedure Pick_Degrees ( ar : in bar; dl : in List;
                           da : in out Integer_Vectors_of_Vectors.Vector ) is

    tmp : List;
    cnt : natural;

  begin
    tmp := Tail_Of(dl);
    cnt := da'first;
    for i in ar'range loop
      if ar(i)
       then da(cnt) := Head_Of(tmp);
            cnt := cnt + 1;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Pick_Degrees;

  function Interchange2 ( p : Laur_Sys; index : integer ) return Laur_Sys is

  -- DESCRIPTION :
  --   Returns a polynomial system where the first equation is interchanged
  --   with the equation given by the index.

    res : Laur_Sys(p'range);

  begin
    if index = p'first
     then res := p;
     else res(res'first) := p(index);
          res(index) := p(p'first);
          res(res'first+1..index-1) := p(p'first+1..index-1);
          res(index+1..res'last) := p(index+1..p'last);
    end if;
    return res;
  end Interchange2;

  procedure Interchange2 ( adl : in out Array_of_Lists; index : integer ) is

  -- DESCRIPTION :
  --   Interchanges the first list with the list given by the index.

    tmp : List;

  begin
    if index /= adl'first
     then tmp := adl(adl'first);
          adl(adl'first) := adl(index);
          adl(index) := tmp;
    end if;
  end Interchange2;

  function Permute ( perm : Integer_Vectors.Vector; p : Laur_Sys )
                   return laur_Sys is

  -- DESCRIPTION :
  --   Returns a permuted polynomial system.

    res : Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := p(perm(i));
    end loop;
    return res;
  end Permute;

  function Initial_Degrees ( p : Poly ) return Degrees is

   -- DESCRIPTION :
   --   Returns the initial degrees of the polynomial p.

    init : Degrees;

    procedure Init_Term ( t : in Term; cont : out boolean ) is
    begin
      init := new Integer_Vectors.Vector'(t.dg.all);
      cont := false;
    end Init_Term;
    procedure Initial_Term is new Visiting_Iterator (Init_Term);

  begin
    Initial_Term(p);
    return init;
  end Initial_Degrees;

  procedure Binomial ( p : in Poly; d : out Integer_Vectors.Link_to_Vector;
                       k : in out integer; c : in out double_complex ) is

   -- DESCRIPTION :
   --   p consists of two terms, this procedure gets the degrees d and
   --   the constant c of the binomial equation.

    first : boolean := true;
    dd : Degrees;

    procedure Scan_Term ( t : in Term; cont : out boolean ) is
    begin
      if first
       then dd := new Integer_Vectors.Vector'(t.dg.all);
	    c := t.cf;
	    first := false;
       else k := dd'first - 1;
	    for i in dd'range loop
	      dd(i) := dd(i) - t.dg(i);
	      if k < dd'first and then dd(i) /= 0
  	       then k := i;
              end if;
            end loop;
	    c := -t.cf/c;
      end if;
      cont := true;
    end Scan_Term;
    procedure Scan_Terms is new Visiting_Iterator (Scan_Term);

  begin
    Scan_Terms(p);
    d := new Integer_Vectors.Vector'(dd.all);
    Integer_Vectors.Clear(Integer_Vectors.Link_to_Vector(dd));
  end Binomial;

  procedure Normalize ( p : in Laur_Sys; dl : in out List; 
			wp : in out Laur_Sys; shifted : out boolean ) is

  -- DESCRIPTION :
  --   Makes sure that the first element of dl contains all zeroes.

  -- REQUIRED :
  --   The list dl is not empty.

  -- ON ENTRY :
  --   p           a Laurent polynomial system;
  --   dl          the power list of p(p'first).

  -- ON RETURN :
  --   dl          the first element of contains all zeroes,
  --                therefore dl has been shifted;
  --   wp          a Laurent polynomial system where
  --                dl is the power list of wp(wp'first);
  --   shifted     is true if dl has been shifted.

    use Integer_Vectors;

    first : Link_to_Vector := Head_Of(dl);
    nullvec : Vector(first'range) := (first'range => 0); 
    shiftvec : Link_to_Vector;

  begin
    if not Is_In(dl,nullvec)
     then shiftvec := Graded_Max(dl);
          Shift(dl,shiftvec);
          Clear(shiftvec);
	  Copy(p,wp);
	  for i in p'range loop
	    Shift(wp(i));
          end loop;
     else wp := p;
          shifted := false;
    end if;
    Move_to_Front(dl,nullvec);
  end Normalize;

  function Evaluate ( p : Poly; x : double_complex; k : integer )
                    return Poly is

  -- DESCRIPTION :
  --   This function returns a polynomial where the kth unknown has
  --   been replaced by x.
  --   It is important to use this function as the term order in p
  --   must remain the same!

    res : Poly;

    procedure Evaluate_Term ( t : in out Term; cont : out boolean ) is
      fac : double_complex;
      pow : natural;
    begin
      if t.dg(k) < 0
       then fac := CMPLX(1.0)/x;
	    pow := -t.dg(k);
       else fac := x;
	    pow := t.dg(k);
      end if;
      for i in 1..pow loop
	t.cf := t.cf * fac;
      end loop;
      declare
        tmp : constant Integer_Vectors.Vector := t.dg.all;
      begin
        Clear(t);
        t.dg := new Integer_Vectors.Vector'(Reduce(tmp,k));
      end;
      cont := true;
    end Evaluate_Term;
    procedure Evaluate_Terms is new Changing_Iterator (Evaluate_Term);

  begin
    Copy(p,res);
    Evaluate_Terms(res);
    return res;
  end Evaluate;

  function Re_Arrange ( p : poly ) return Poly is

  -- DESCRIPTION :
  --   Returns a polynomial whose terms are sorted
  --   in graded lexicographical ordening.
   
    res : Poly := Null_Poly;

    procedure Re_Arrange_Term ( t : in Term; cont : out boolean ) is
    begin
      Plus_Term(res,t);
      cont := true;
    end Re_Arrange_Term;
    procedure Re_Arrange_Terms is new Visiting_Iterator (Re_Arrange_Term);

  begin
    Re_Arrange_Terms(p);
    return res;
  end Re_Arrange;

  function Substitute ( p : Poly; v : Complex_Vectors.Vector;
                        k : integer ) return Poly is

   -- DESCRIPTION :
   --   Substitutes the values in v into the polynomial p,
   --   starting at the last unknowns of p.

    init : Degrees := Initial_Degrees(p);
    index : integer := init'last;
    res,tmp : Poly;
  begin
    if index = k
     then index := index - 1;
	  res := Evaluate(p,v(v'last),index);
     else res := Evaluate(p,v(v'last),index);
    end if;
    for i in reverse v'first..(v'last-1) loop
      index := index - 1;
      if index = k
       then index := index - 1;
      end if;
      tmp := Evaluate(res,v(i),index);
      Clear(res); Copy(tmp,res); Clear(tmp);
    end loop;
    Integer_Vectors.Clear(Integer_Vectors.Link_to_Vector(init));
    return res;
  end Substitute;

  procedure Refine_and_Concat 
	       ( p : in Laur_Sys;
                 newsols,sols,sols_last : in out Solution_List ) is

  -- DESCRIPTION :
  --   This procedure refines the solutions of a given
  --   polynomial system and adds them to the solution list.
  --   This can be very useful for eliminating rounding errors
  --   after transformating the solutions.

    use Complex_Polynomial_Systems;

    pp : Poly_Sys(p'range) := Laurent_to_Polynomial_System(p);
    numb : natural := 0;

  begin
    Silent_Root_Refiner(pp,newsols,10.0**(-12),10.0**(-12),10.0**(-8),numb,5);
    Concat(sols,sols_last,newsols);
    Clear(pp); Shallow_Clear(newsols);
  end Refine_and_Concat;

  procedure Refine_and_Concat 
	       ( file : in file_type; p : in Laur_Sys;
	         newsols,sols,sols_last  : in out Solution_List ) is

  -- DESCRIPTION :
  --   This procedure refines the solutions of a given
  --   polynomial system and adds them to the solution list.
  --   This can be very useful for eliminating rounding errors
  --   after transformating the solutions.

    use Complex_Polynomial_Systems;

    pp : Poly_Sys(p'range) := Laurent_to_Polynomial_System(p);
    numb : natural := 0;

  begin
    Reporting_Root_Refiner
      (file,pp,newsols,10.0**(-12),10.0**(-12),10.0**(-8),numb,5,false);
    Concat(sols,sols_last,newsols);
    Clear(pp); Shallow_Clear(newsols);
  end Refine_and_Concat;

  procedure Write_Direction ( file : in file_type;
                              v : in Integer_Vectors.Link_to_Vector ) is
  begin
    put(file,"++++  considering direction "); put(file,v);
    put_line(file,"  ++++");
  end Write_Direction;

-- INTERMEDIATE LAYER :

  procedure Mixed_Continue
             ( file : in file_type; p : in Laur_Sys;
               k : in integer; m : in Integer_Vectors.Vector;
               sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   This continuation routine computes a part of the solution list 
  --   of a Laurent polynomial system.

  -- ON ENTRY :
  --   file      a file where the intermediate results are written;
  --   p         the transformed Laurent polynomial system to be solved;
  --   k         the index;
  --   m         m(k) = p1(v), m(k/=l) = Maximal_Degree(p(l),k);
  --   sols      the start solutions.

  -- ON RETURN :
  --   sols      the computed solutions.

    h : Laur_Sys(p'range);
    hp : Poly_Sys(h'range);
    hpe : Eval_Poly_Sys(hp'range);
   -- he : Eval_Laur_Sys(h'range);
    j : Jacobi(p'range,p'first..p'last+1);
    je : Eval_Jacobi(j'range(1),j'range(2));
   -- j : Jacobi(h'range,h'first..h'last+1);
   -- je : Eval_Jacobi(j'range(1),j'range(2));

    function Construct_Homotopy
                 ( p : Laur_Sys; k : integer; m : Integer_Vectors.Vector )
                 return Laur_Sys is

      res : Laur_Sys(p'range);
      ran : double_complex;
      re : boolean;
      zeroes : Degrees := new Integer_Vectors.Vector'(p'range => 0);

      function Construct_First_Polynomial
                 ( pp : Poly; max : integer ) return Poly is

        r : Poly := Null_Poly;

        procedure Construct_Term ( t : in Term; cont : out boolean ) is

          rt : Term;

        begin
          rt.cf := t.cf;
          rt.dg := new Integer_Vectors.Vector(t.dg'first..t.dg'last+1);
          rt.dg(t.dg'range) := t.dg.all;
          rt.dg(k) := -t.dg(k) + max;
          if Equal(t.dg,zeroes)
           then rt.dg(rt.dg'last) := 0;
                re := ( IMAG_PART(rt.cf) + 1.0 = 1.0 );
           else rt.dg(rt.dg'last) := rt.dg(k);
          end if;
          Plus_Term(r,rt);
          Clear(rt);
          cont := true;
        end Construct_Term;
        procedure Construct_Terms is new Visiting_Iterator(Construct_Term);

      begin
        Construct_Terms(pp);
        Integer_Vectors.Clear(Integer_Vectors.Link_to_Vector(zeroes));
        return r;
      end Construct_First_Polynomial;

      function Construct_Polynomial ( pp : Poly; max : integer ) return Poly is

        r : Poly := Null_Poly;

        procedure Construct_Term ( t : in Term; cont : out boolean ) is

          rt : Term;

        begin
          rt.cf := t.cf;
          rt.dg := new Integer_Vectors.Vector(t.dg'first..t.dg'last+1);
          rt.dg(t.dg'range) := t.dg.all;
          rt.dg(k) := -t.dg(k) + max;
          rt.dg(rt.dg'last) := rt.dg(k);
          Plus_Term(r,rt);
          Clear(rt);
          cont := true;
        end Construct_Term;
        procedure Construct_Terms is new Visiting_Iterator(Construct_Term);

      begin
        Construct_Terms(pp);
        return r;
      end Construct_Polynomial;

    begin
      res(res'first) := Construct_First_Polynomial(p(p'first),m(m'first));
      for i in p'first+1..p'last loop
        res(i) := Construct_Polynomial(p(i),m(i));
      end loop;
      if re
       then for i in res'range loop
              ran := Random_Number_Generators.Random;
              Mult_Coeff(res(i),ran);
            end loop;
      end if;
      return res;
    end Construct_Homotopy;

    function To_Be_Removed ( flag : in natural ) return boolean is
    begin
      return ( flag /= 1 );
    end To_Be_Removed;
    procedure Extract_Regular_Solutions is new Solutions.Delete(To_Be_Removed);

  begin
    h := Construct_Homotopy(p,k,m);               -- CONSTRUCTION OF HOMOTOPY
    hp := Laurent_to_Polynomial_System(h);
    hpe := Create(hp);
   -- he := Create(h);
    j := Create(hp);
    je := Create(j);
   -- j := Create(h);
   -- je := Create(j);
    declare                                                   -- CONTINUATION 

      use Complex_Vectors,Complex_Matrices;

      function Eval ( x : Vector; t : double_complex ) return Vector is

        xt : Vector(x'first..x'last+1);
        use Complex_Multivariate_Polynomials;
        --use Complex_Multivariate_Laurent_Polynomials;

      begin
        xt(x'range) := x;
        xt(xt'last) := t;
        return Eval(hpe,xt);
        --return Eval(he,xt);
      end Eval;

      function dHt ( x : Vector; t : double_complex ) return Vector is

        xt : Vector(x'first..x'last+1);
        res : Complex_Vectors.Vector(p'range);
        use Complex_Multivariate_Polynomials;
        --use Complex_Multivariate_Laurent_Polynomials;

      begin
        xt(x'range) := x;
        xt(xt'last) := t;
        for i in res'range loop
          res(i) := Eval(je(i,xt'last),xt);
        end loop;
        return res;
      end dHt;

      function dHx ( x : Vector; t : double_complex ) return Matrix is

        xt : Vector(x'first..x'last+1);
        m : Matrix(x'range,x'range);
        use Complex_Multivariate_Polynomials;
        --use Complex_Multivariate_Laurent_Polynomials;

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

      procedure Invert ( k : in integer; sols : in out Solution_List ) is

      -- DESCRIPTION :
      --   For all solutions s in sols : s.v(k) := 1/s.v(k).
      
      tmp : Solution_List := sols;

      begin
        while not Is_Null(tmp) loop
          declare
	    l : Link_to_Solution := Head_Of(tmp);
          begin
            l.v(k) := CMPLX(1.0)/l.v(k);
            l.t := CMPLX(0.0);
          end;
          tmp := Tail_Of(tmp);
        end loop;
      end Invert;

      procedure Cont is new Reporting_Continue(Norm2,Eval,dHt,dHx);

    begin
      Invert(k,sols);
      Cont(file,sols,false);
      Invert(k,sols);
      Extract_Regular_Solutions(sols);
     -- declare
     --   pp : Poly_Sys(p'range);
     --   numb : natural := 0;
     -- begin
     --   pp := Laurent_to_Polynomial_System(p);
     --   Reporting_Root_Refiner
     --     (file,pp,sols,10.0**(-12),10.0**(-12),10.0**(-8),numb,5,false);
     --   Clear(pp); Extract_Regular_Solutions(sols);
     -- end;
    end;
   -- CLEANING :
    Clear(h); Clear(hp); Clear(hpe);
   -- Clear(he);
    Clear(j); Clear(je);
  end Mixed_Continue;

-- THE SOLVERS : 

  function One_Unknown_Solve ( p : Poly ) return Complex_Vectors.Vector is

  -- DESCRIPTION :
  --   Returns the solution vector of p, a polynomial in one unknown.
  --   p will be solved by using the method of Durand-Kerner.

    p1 : Poly := Re_Arrange(p);
    init : Degrees := Initial_Degrees(p1);
    index : integer := init'first;
    min : integer := -Minimal_Degree(p1,index);
    pv : Complex_Vectors.Vector(0..Degree(p1)+min);
    z,res : Complex_Vectors.Vector(1..pv'last);
    maxsteps : constant natural := 10;
    eps : constant double_float := 10.0**(-10);
    nb : integer := pv'last + 1;

    procedure Store_Coeff ( t : in Term; cont : out boolean ) is
    begin
      nb := nb - 1;
      if t.dg(index) = (nb - min)
       then pv(nb) := t.cf;
       else for i in reverse nb..(t.dg(index)+min+1) loop
	      pv(i) := CMPLX(0.0);
            end loop;
	    nb := t.dg(index) + min;
	    pv(nb) := t.cf;
      end if;
      cont := true;
    end Store_Coeff;
    procedure Polynomial_To_Vector is new Visiting_Iterator (Store_Coeff);

    procedure Write ( step : in natural; z,res : in Complex_Vectors.Vector ) is
    begin
      null;  -- no output desired during the iterations
    end Write;
    procedure DuKe is new Durand_Kerner (Write);

  begin
    Integer_Vectors.Clear(Integer_Vectors.Link_to_Vector(init));
    Polynomial_To_Vector(p1);
    Clear(p1);
    for i in z'range loop
      z(i) := Random_Number_Generators.Random;
      res(i) := z(i);
    end loop;
    DuKe(pv,z,res,maxsteps,eps,nb);
    return z;
  end One_Unknown_Solve;

  procedure One_Unknown_Solve ( p : in Poly; sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   If p is a polynomial in one unknown,
  --   p can be solved efficiently by the application of Durand-Kerner.

    init : Degrees := Initial_Degrees(p);

    function Make_Solutions ( x : in Complex_Vectors.Vector )
                            return Solution_List is
      res,res_last : Solution_List;
      s : Solution(1);
    begin
      s.m := 1;
      s.t := CMPLX(0.0);
      for i in x'range loop
	s.v(init'first) := x(i);
	Append(res,res_last,s);
      end loop;
      return res;
    end Make_Solutions;

  begin
    sols := Make_Solutions(One_Unknown_Solve(p));
    Integer_Vectors.Clear(Integer_Vectors.Link_to_Vector(init));
  end One_Unknown_Solve;

  procedure Two_Terms_Solve
               ( file : in file_type; p : in Laur_Sys;
                 tv : in Tree_of_Vectors; bkk : out natural;
                 sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   The first polynomial of p consists of two terms.
  --   A binomial system can be solved efficiently by 
  --   transforming and using de Moivre's rule.

    d : Integer_Vectors.Link_to_Vector;
    c : double_complex := CMPLX(0.0);
    k : natural := 0;
    sols_last : Solution_List;

  begin
   -- put_line(file,"Applying Two_Terms_Solve on "); Write(file,p);
    Binomial(p(p'first),d,k,c);
    if k < d'first
     then bkk := 0;
	  Integer_Vectors.Clear(d); return;
     elsif ( c + CMPLX(1.0) = CMPLX(1.0) )
         then bkk := 0;
	      Integer_Vectors.Clear(d); return;
         elsif d(k) < 0
             then Integer_Vectors.Min_Vector(d);
                  c := CMPLX(1.0)/c;
    end if;
  --  if p'first = p'last
  --   then
  --     declare
  --       sol : Solution(1);
  --     begin
  --       sol.m := 1;
  --       sol.t := CMPLX(0.0);
  --       for i in 1..d(k) loop
  --         sol.v(sol.v'first) := de_Moivre(d(k),i,c);
  --         Append(sols,sols_last,sol);
  --       end loop;
  --       bkk := d(k);
  --     end;
  --   else
    declare
      t : Transfo := Rotate(d,k);
      tmp_bkk : natural := 0;
    begin
      Apply(t,d);
      declare
        v : Complex_Vectors.Vector(1..d(k));
        tmp : Poly;
        rtp : Laur_Sys(p'first..(p'last-1));
        rtp_bkk : natural;
        rtp_sols : Solution_List;
      begin
        for i in v'range loop
          v(i) := de_Moivre(d(k),i,c);
          for j in rtp'range loop
            tmp := Transform(t,p(j+1));
            rtp(j) := Evaluate(tmp,v(i),k);
            Clear(tmp);
          end loop;
          Solve(file,rtp,tv,rtp_bkk,rtp_sols);
          Clear(rtp);
          tmp_bkk := tmp_bkk + rtp_bkk;
          Insert(v(i),k,rtp_sols);
          Transform(t,rtp_sols);
         --Concat(sols,sols_last,rtp_sols);
          Refine_and_Concat(file,p,rtp_sols,sols,sols_last);
        end loop;
      end;
      Clear(t);
      bkk := tmp_bkk;
    end;
  --  end if;
    Integer_Vectors.Clear(d);
   -- put_line(file,"The solutions found : ");  put(file,sols);
  end Two_Terms_Solve;

  function Project_on_First_and_Solve
                  ( p : Poly; k : integer; sols : Solution_List )
                  return Solution_List is

  -- ON ENTRY :
  --   p          a Laurent polynomial in n unknowns x1,..,xk,..,xn;
  --   sols       contains values for x1,..,xn, except xk.

   -- ON RETURN :
   --   a solution list for p, obtained after substition of the values
   --   for x1,..,xn into the polynomial p.

    tmp,res,res_last : Solution_List;
    init : Degrees := Initial_Degrees(p);

  begin
   -- put_line(file,"Calling Project_on_First_and_Solve");
   -- put_line(file," with polynomial : ");
   -- put(file,Laurent_Polynomial_to_Polynomial(p)); new_line(file);
    tmp := sols;
    while not Is_Null(tmp) loop
      declare
	p1 : Poly := Substitute(p,Head_Of(tmp).v,k);
	sols1 : Solution_List;
      begin
       -- put(file,"k : "); put(file,k,1); new_line(file);
       -- put(file,"v : "); put(file,Head_Of(tmp).v,3,3,3); new_line(file);
       -- put_line(file,"After substitution : "); Write(file,p1);
	if Number_of_Terms(p1) < 2
	 then null;
	 elsif Number_of_Terms(p1) = 2
	     then
               declare
	         d : Integer_Vectors.Link_to_Vector;
	         l : natural := 0;
	         c : double_complex := CMPLX(0.0);
               begin
	         Binomial(p1,d,l,c);
	         if l < init'first 
	          then null;
	          elsif ( c + CMPLX(1.0) = CMPLX(1.0) )
	              then null;
	              else
                        if d(l) < 0
		         then d(l) := -d(l); c := CMPLX(1.0)/c;
                        end if;
                        declare
	                  v : Complex_Vectors.Vector(1..d(l));
                        begin
	                  for i in v'range loop
	                    v(i) := de_Moivre(d(l),i,c);
                          end loop;
	                  sols1 := Insert(v,k,Head_Of(tmp).all);
                         -- put_line(file,"Sols1 after Binomial :");
                         -- put(file,sols1);
	     	          Concat(res,res_last,sols1);
			  Shallow_Clear(sols1);
                        end;
                 end if;
		 Integer_Vectors.Clear(d);
               end;
	     else
               sols1 := Insert(One_Unknown_Solve(p1),k,Head_Of(tmp).all);
	       Concat(res,res_last,sols1);
              -- put_line(file,"Sols1 :"); put(file,sols1);
               Shallow_Clear(sols1);
        end if;
        Clear(p1);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Integer_Vectors.Clear(Integer_Vectors.Link_to_Vector(init));
    return res;
  end Project_on_First_and_Solve;

  procedure Project_and_Solve
                ( file : in file_type; p : in Laur_Sys; k : in integer;
                  m : in out Integer_Vectors.Vector;
                  nd : in node; bkk : out natural;
                  sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Solves the projected start system of p along a direction v.

  -- ON ENTRY :
  --   file       a file to write intermediate results on;
  --   p          a Laurent polynomial system;
  --   k          entry in the degrees of p;
  --   m          m(m'first) equals Maximal_Support(p(p'first),v) > 0;
  --   nd         a node in the tree of useful directions.

  -- ON RETURN :
  --   m          m(m'first+1..m'last) contains the maximal degrees
  --               of the last n-1 equations of p in xk;
  --   bkk        the BKK bound of the projected system;
  --   sols       the solutions of the projected system.

    g_v : Laur_Sys(p'first..(p'last-1));
    bkk_g_v : natural;
    sols_g_v : Solution_List;

  begin
   -- put_line(file,"Applying Project_and_Solve on"); Write(file,p);
    for i in g_v'range loop
      m(i+1) := Maximal_Degree(p(i+1),k);
      g_v(i) := Face(k,m(i+1),p(i+1));
      Reduce(k,g_v(i));
    end loop;
    if (nd.ltv = null) or else Is_Null(nd.ltv.all)
     then Solve(file,g_v,bkk_g_v,sols_g_v);
     else Solve(file,g_v,nd.ltv.all,bkk_g_v,sols_g_v);
    end if;
   -- put(file,"After Solve (without tv) bkk_g_v = "); put(file,bkk_g_v,1);
   -- new_line(file);
    declare
      p0 : Poly := Re_Arrange(p(p'first));
      p1 : Poly := Face(k,m(m'first),p0);
      cnst : Term;
    begin
      cnst.dg := new Integer_Vectors.Vector'(p'range => 0);
      if Coeff(p1,cnst.dg) = CMPLX(0.0)
       then cnst.cf := Coeff(p0,cnst.dg);
            Plus_Term(p1,cnst);
      end if;
      Integer_Vectors.Clear(Integer_Vectors.Link_to_Vector(cnst.dg));
      Clear(p0);
      sols := Project_on_First_and_Solve(p1,k,sols_g_v);
     -- sols := Project_on_First_and_Solve(file,p1,k,sols_g_v);
      Set_Continuation_Parameter(sols,CMPLX(0.0));
      Clear(p1);
    end;
    bkk := m(m'first)*bkk_g_v;
    Clear(sols_g_v);
    Clear(g_v);
   -- put_line(file,"The solutions found : "); put(file,sols);
  end Project_and_Solve;

  procedure Unmixed_Solve
                ( file : in file_type; p : in Laur_Sys; dl : in List;
                  tv : in Tree_of_Vectors;
                  bkk : out natural; sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Solves a Laurent polynomial system where all polytopes are the same.

  -- ON ENTRY :
  --   file       a file to write intermediate results on;
  --   p          a Laurent polynomial system;
  --   dl         the list of powers for p;
  --   tv         the tree of degrees containing the useful directions.

  -- ON RETURN :
  --   bkk        the bkk bound;
  --   sols       the list of solutions.

    sols_last : Solution_List;
    tmp_bkk : natural := 0;
    tmp : Tree_of_Vectors := tv;

  begin
    tmp_bkk := 0;
    tmp := tv;
    while not Is_Null(tmp) loop
      declare
        nd : node := Head_Of(tmp);
        v : Integer_Vectors.Link_to_Vector := nd.d;
        i : integer := Pivot(v);
        pv : integer := Maximal_Support(dl,v.all);
        t : Transfo := Build_Transfo(v,i);
        tp : Laur_Sys(p'range) := Transform(t,p);
        bkk_tp : natural;
        sols_tp : Solution_List;
        max : Integer_Vectors.Vector(p'range);
      begin
        Write_Direction(file,v);
        max(max'first) := pv;
       -- if (nd.ltv = null) or else Is_Null(nd.ltv.all)
       --  then Projected_Solve(file,tp,i,max,bkk_tp,sols_tp);
       --  else Projected_Solve
       --           (file,tp,i,max,nd.ltv.all,bkk_tp,sols_tp);
       -- end if;
        Project_and_Solve(file,tp,i,max,nd,bkk_tp,sols_tp);
        Mixed_Continue(file,tp,i,max,sols_tp);
        tmp_bkk := tmp_bkk + bkk_tp;
        Transform(t,sols_tp);
        --Concat(sols,sols_last,sols_tp);
        Refine_and_Concat(file,p,sols_tp,sols,sols_last);
        Clear(t); Clear(tp);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    bkk := tmp_bkk;
  end Unmixed_Solve;

  procedure Unmixed_Solve
                ( file : in file_type; p : in Laur_Sys; dl : in List;
                  bkk : out natural; sols : in out Solution_List ) is

    tv : Tree_of_Vectors;

  begin
    Volumes.Volume(p'last,dl,tv,bkk);
    Unmixed_Solve(file,p,dl,tv,bkk,sols);
    Clear(tv);
  end Unmixed_Solve;

  procedure Mixed_Solve 
               ( file : in file_type; p : in Laur_Sys;
                 adl : in out Array_of_Lists; tv : in Tree_of_Vectors;
                 bkk : out natural; sols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Computes the solutions of the Laurent polynomial system p,
  --   where p has more than one equation.

  -- NOTE :
  --   This procedure mirrors the procedure Minkowski_Sum in the body
  --   of the package Volumes.

    tmp_bkk,len : natural;
    tmp : Tree_of_Vectors;
    index : integer := Index2(adl);
    wp : Laur_Sys(p'range);
    sols_last : Solution_List;
    shifted : boolean;
    perm,mix : Integer_Vectors.Link_to_Vector;

  begin
    Interchange2(adl,index);
    len := Length_Of(adl(adl'first));
   -- put_line(file,"Applying Mixed_Solve on"); Write(file,p);
    if len = 2
     then wp := Interchange2(p,index);
          Two_Terms_Solve(file,wp,tv,bkk,sols);
     else
       if len > 2
        then -- INITIALIZATION :
          Mixture(adl,perm,mix);
          wp := Permute(perm.all,p);
          declare
            zeroes : Degrees := new Integer_Vectors.Vector'(p'range => 0);
            tmpwpi : Poly;
          begin
            if Coeff(wp(wp'first),zeroes) = CMPLX(0.0)
             then shifted := true;
                 -- wp(wp'first) := Shift(p(p'first));
                  Copy(p(index),tmpwpi); wp(wp'first) := tmpwpi;
                  Shift(wp(wp'first));
             else shifted := false;
            end if;
            Integer_Vectors.Clear(Integer_Vectors.Link_to_Vector(zeroes));
          end;
         -- MIXED HOMOTOPY CONTINUATION :
          tmp_bkk := 0;
          tmp := tv;
          while not Is_Null(tmp) loop
            declare
              nd : node := Head_Of(tmp);
              v : Integer_Vectors.Link_to_Vector := nd.d;
              k : integer := Pivot(v);
              pv : integer := Maximal_Support(wp(wp'first),v);
              t : Transfo := Build_Transfo(v,k);
              twp : Laur_Sys(wp'range) := Transform(t,wp);
              bkk_twp : natural;
              sols_twp : Solution_List;
              m : Integer_Vectors.Vector(wp'range);
            begin
              Write_Direction(file,v);
              m(m'first) := pv;
             -- if (nd.ltv = null) or else Is_Null(nd.ltv.all)
             --  then Projected_Solve(file,twp,k,m,bkk_twp,sols_twp);
             --  else Projected_Solve
             --           (file,twp,k,m,nd.ltv.all,bkk_twp,sols_twp);
             -- end if;
              Project_and_Solve(file,twp,k,m,nd,bkk_twp,sols_twp);
              Mixed_Continue(file,twp,k,m,sols_twp);
              tmp_bkk := tmp_bkk + bkk_twp;
              Transform(t,sols_twp);
             --Concat(sols,sols_last,sols_twp);
              Refine_and_Concat(file,wp,sols_twp,sols,sols_last);
              Clear(t); Clear(twp);
            end;
            tmp := Tail_Of(tmp);
          end loop;
          bkk := tmp_bkk;
         -- CLEANING UP :
          Integer_Vectors.Clear(perm); Integer_Vectors.Clear(mix);
          if shifted
           then Clear(wp(wp'first));
          end if;
        else -- len < 2
          bkk := 0;
       end if;
    end if;
  end Mixed_Solve;

  procedure Mixed_Solve
               ( file : in file_type; p : in Laur_Sys;
                 adl : in out Array_of_Lists; bkk : out natural;
                 sols : in out Solution_List ) is

    tv : Tree_of_Vectors;

  begin
    Volumes.Mixed_Volume(adl'last,adl,tv,bkk);
    Mixed_Solve(file,p,adl,tv,bkk,sols);
    Clear(tv);
  end Mixed_Solve;

  procedure Solve ( file : in file_type; p : in Laur_Sys;
                    bkk : out natural; sols : in out Solution_List ) is

    al : Array_of_Lists(p'range) := Construct_Power_Lists(p);
    tv : Tree_of_Vectors;

  begin
    Volumes.Mixed_Volume(p'last,al,tv,bkk);
    Solve(file,p,tv,bkk,sols);
    Deep_Clear(al); Clear(tv);
  end Solve;

  procedure Solve ( file : in file_type; p : in Laur_Sys;
                    tv : in Tree_of_Vectors; bkk : out natural;
                    sols : in out Solution_List ) is

  -- NOTE :
  --   This procedure mirrors the procedure Volumes.Mixed_Volume,
  --   with a tree of useful directions on entry.

  begin
    if p'first = p'last
     then One_Unknown_Solve(p(p'first),sols);
          bkk := Length_Of(sols);
     else --if Is_Fewnomial_System(p)
          -- then
          --   declare
          --     fail : boolean;
          --   begin
          --     Fewnomials.Solve(p,sols,fail);
          --     if fail
          --      then bkk := 0;  Clear(sols);
          --      else bkk := Length_Of(sols);
          --     end if;
          --   end;
          -- else
      declare
        adl : Array_of_Lists(p'range) := Construct_Power_Lists(p);
      begin
        if All_Equal(adl)
         then 
           for i in (adl'first+1)..adl'last loop
             Deep_Clear(adl(i));
           end loop;
           declare
             wp : Laur_Sys(p'range);
             shifted : boolean;
           begin
             Normalize(p,adl(adl'first),wp,shifted);
             if Is_Null(tv)
              then Unmixed_Solve(file,wp,adl(adl'first),bkk,sols);
              else Unmixed_Solve(file,wp,adl(adl'first),tv,bkk,sols);
             end if;
             if shifted
              then Clear(wp);
             end if;
           end;
         elsif Is_Null(tv)
             then Mixed_Solve(file,p,adl,bkk,sols);
             else Mixed_Solve(file,p,adl,tv,bkk,sols);
        end if;
      end;
    end if;
  end Solve;

end Mixed_Homotopy_Continuation;
