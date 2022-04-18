with Lists;
with unchecked_deallocation;
with Floating_Point_Numbers;            use Floating_Point_Numbers;
with Natural_Vectors;
with Complex_Vectors,Complex_Numbers;   use Complex_Vectors,Complex_Numbers;
with Complex_Matrices,Integer_Vectors;  use Complex_Matrices;
with Complex_Linear_System_Solvers;     use Complex_Linear_System_Solvers;

package body Random_Product_System is

-- DATA STRUCTURES :

  package List_of_Vectors is new Lists(Link_to_Vector);
  type Equation_List is new List_of_Vectors.List;
  
  type Equation is record
    first,last : Equation_List;
  end record;

  type Equations is array(positive range <>) of Equation;
  type Link_To_Equations is access Equations;
  procedure free is new unchecked_deallocation (Equations,Link_To_Equations);

-- INTERNAL DATA :
 
  rps : Link_To_Equations;

  Bad_Condition : constant double_float := 10.0**(-12);

-- OPERATIONS :

  procedure Init ( n : in natural ) is
  begin
    rps := new Equations(1..n);
  end Init;

  procedure Clear ( eql : in out Equation_List ) is

    temp : Equation_List := eql;
    lv : Link_To_Vector;

  begin
    while not Is_Null(temp) loop
      lv := Head_Of(temp);
      Clear(lv);
      temp := Tail_of(temp);
    end loop;
    List_Of_Vectors.Clear(List_Of_Vectors.List(eql));
  end Clear;

  procedure Clear ( eq : in out Equation ) is
  begin
    Clear(eq.first);
    -- eq.last is just a pointer to the last element of eq.first;
    -- if eq.first disappears, then also eq.last does
  end Clear;

  procedure Clear ( eqs : in out Equations ) is
  begin
    for i in eqs'range loop
      Clear(eqs(i));
    end loop;
  end Clear;

  procedure Clear is
  begin
    if rps /= null
     then for i in rps'range loop
            Clear(rps(i));
          end loop;
          free(rps);
    end if;
  end Clear;

  procedure Add_Hyperplane ( i : in natural; h : in Vector ) is

    eqi : Equation renames rps(i);
    lh : Link_To_Vector := new Vector'(h);

  begin
    if Is_Null(eqi.first)
     then Construct(lh,eqi.first);
          eqi.last := eqi.first;
     else declare 
            temp : Equation_List;
          begin
            Construct(lh,temp);
            Swap_Tail(eqi.last,temp);
            eqi.last := Tail_Of(eqi.last);
          end;
    end if;
  end Add_Hyperplane;

  function Dimension return natural is
  begin
    if rps = null
     then return 0;
     else return rps'last;
    end if;
  end Dimension;

  function Number_of_Hyperplanes ( i : natural ) return natural is
  begin
    if rps = null
     then return 0;
     else return Length_Of(rps(i).first);
    end if;
  end Number_Of_Hyperplanes;
 
  function Get_Hyperplane ( i,j : in natural ) return Link_to_Vector is

    nulvec : Link_to_Vector := null;

  begin
    if rps = null
     then return nulvec;
     else declare
            eqi : Equation_List := rps(i).first;
            count : natural := 1;
          begin
            while not Is_Null(eqi) loop
              if count = j
               then return Head_Of(eqi);
               else count := count + 1;
                    eqi := Tail_Of(eqi);
              end if;
            end loop;
          end;
          return nulvec;
    end if;
  end Get_Hyperplane;

  function Get_Hyperplane ( i,j : in natural ) return Vector is

    lres : Link_to_Vector := Get_Hyperplane(i,j);
    nulvec : Vector(0..0);
 
  begin
    if lres = null
     then return nulvec;
     else return lres.all;
    end if;
  end Get_Hyperplane;

  procedure Change_Hyperplane ( i,j : in natural; h : in Vector ) is
  begin
    if rps = null
     then return;
     else declare
            eqi : Equation_List := rps(i).first;
	    lv : Link_To_Vector;
            count : natural := 1;
          begin
            while not Is_Null(eqi) loop
              if count = j
               then lv := Head_Of(eqi);
		    for k in h'range loop
		      lv(k) := h(k);
                    end loop;
                    return;
               else count := count + 1;
                    eqi := Tail_Of(eqi);
              end if;
            end loop;
          end;
    end if;
  end Change_Hyperplane;

  procedure Solve ( i,n : in natural; sols,sols_last : in out Solution_List;
                    a : in out Matrix; b : in out Vector; 
                    nb : in out natural ) is
  begin
    if i > n
     then declare
            aa : Matrix(a'range(1),a'range(2));
            bb : Vector(b'range);
            rcond : double_float;
            ipvt : Integer_Vectors.Vector(b'range);
          begin
            for k in a'range(1) loop
              for l in a'range(2) loop
                aa(k,l) := a(k,l);
              end loop;
              bb(k) := b(k);
            end loop;
            lufco(aa,n,ipvt,rcond);
            nb := nb + 1;
            if abs(rcond) > Bad_Condition
             then lusolve(aa,n,ipvt,bb);
                  declare
                    s : Solution(n);
                  begin
                    s.t := CMPLX(0.0);
                    s.m := 1;
                    s.v := bb;
                    s.err := 0.0; s.rco := rcond; s.res := 0.0;
                    Append(sols,sols_last,s);
                  end;
            end if;
          exception
            when numeric_error => return;
          end;
     else declare
            eqi : Equation_List := rps(i).first;
            h : Vector(0..n);
            count : natural := 0;
          begin
            while not Is_Null(eqi) loop
              count := count + 1;
              h := Head_Of(eqi).all;
              b(i) := -h(0);
              for j in 1..n loop
                a(i,j) := h(j);
              end loop;
              Solve(i+1,n,sols,sols_last,a,b,nb);
              eqi := Tail_Of(eqi);
            end loop;
          end;
    end if;
  end Solve;

  procedure Solve ( sols : in out Solution_List; nl : out natural ) is

    n : natural := rps'last;
    m : Matrix(1..n,1..n);
    v : Vector(1..n);
    num : natural := 0;
    last : Solution_List;

  begin
    for i in 1..n loop
      for j in 1..n loop
        m(i,j) := CMPLX(0.0);
      end loop;
      v(i) := CMPLX(0.0);
    end loop;
    Solve(1,n,sols,last,m,v,num);
    nl := num;
  end Solve;

  procedure Solve ( sols : in out Solution_List; nl : out natural;
                    l : in List ) is

    n : natural := rps'last;
    m : Matrix(1..n,1..n);
    v : Vector(1..n);
    num : natural := 0;
    temp : List := l;
    pos : Integer_Vectors.Vector(1..n);
    stop : boolean := false;
    last : Solution_List;

    procedure PSolve ( i,n : in natural; sols,sols_last : in out Solution_List;
                       a : in out Matrix; b : in out Vector;
		       nb : in out natural ) is
    begin
      if i > n
       then declare
              aa : Matrix(a'range(1),a'range(2));
              bb : Vector(b'range);
              rcond : double_float;
              ipvt : Integer_Vectors.Vector(b'range);
            begin
              for k in a'range(1) loop
                for l in a'range(2) loop
                  aa(k,l) := a(k,l);
                end loop;
                bb(k) := b(k);
              end loop;
              lufco(aa,n,ipvt,rcond);
              nb := nb + 1;
              if abs(rcond) > Bad_Condition
               then lusolve(aa,n,ipvt,bb);
                    declare
                      s : Solution(n);
                    begin
                      s.t := CMPLX(0.0);
                      s.m := 1;
                      s.v := bb;
                      s.err := 0.0; s.rco := rcond; s.res := 0.0;
                      Append(sols,sols_last,s);
                    end;
              end if;
              if Is_Null(temp)
               then stop := true;
               else pos := Head_Of(temp).all;
                    temp := Tail_Of(temp);
              end if;
            exception
              when numeric_error => return;
            end;
       else declare
              eqi : Equation_List := rps(i).first;
              h : Vector(0..n);
              count : natural := 0;
            begin
              while not Is_Null(eqi) loop
                count := count + 1;
    	        if count = pos(i)
                 then h := Head_Of(eqi).all;
                      b(i) := -h(0);
                      for j in 1..n loop
                        a(i,j) := h(j);
                      end loop;
                      --put("eq"); put(i,1); put(count,1); put(" ");
                      PSolve(i+1,n,sols,sols_last,a,b,nb);
                end if;
	      	exit when stop;
                eqi := Tail_Of(eqi);
              end loop;
            end;
      end if;
    end PSolve;

  begin
    if not Is_Null(temp)
     then pos := Head_Of(temp).all;
          temp := Tail_Of(temp);
          for i in 1..n loop
            for j in 1..n loop
              m(i,j) := CMPLX(0.0);
            end loop;
            v(i) := CMPLX(0.0);
          end loop;
          PSolve(1,n,sols,last,m,v,num);
          nl := num;
    end if;
  end Solve;

  function Polynomial ( h : Vector ) return Poly is
 
    res : Poly;
    t : Term;
    n : natural := h'last;

  begin
    for j in 0..n loop
      if h(j) /= CMPLX(0.0)
       then t.dg := new Natural_Vectors.Vector'(1..n => 0);
            t.cf := h(j);
            if j /= 0
             then t.dg(j) := 1;
            end if;
            Plus_Term(res,t);
            Natural_Vectors.Clear(Natural_Vectors.Link_to_Vector(t.dg));
      end if;
    end loop;
    return res;
  end Polynomial;

  function Create ( i : in natural ) return Poly is

    eql : Equation_List := rps(i).first;
    hyp,res : Poly := Null_Poly;

  begin
    while not Is_Null(eql) loop
      hyp := Polynomial(Head_Of(eql).all);
      if res = Null_Poly
       then Copy(hyp,res);
       else Mult_Poly(res,hyp);
      end if;
      Clear(hyp);
      eql := Tail_Of(eql);
    end loop;
    return res;
  end Create;

  function Polynomial_System return Poly_Sys is

    res : Poly_Sys(rps'range);

  begin
    for i in rps'range loop
      res(i) := Create(i);
    end loop;
    return res;
  end Polynomial_System;

end Random_Product_System;
