with Floating_Point_Numbers;           use Floating_Point_Numbers;
with Complex_Vectors,Complex_Norms;    use Complex_Vectors,Complex_Norms;
with Complex_Matrices;                 use Complex_Matrices;
with Complex_Polynomial_Systems_io;    use Complex_Polynomial_Systems_io;
with Polynomial_Randomizers;           use Polynomial_Randomizers;

with Scaling;                          use Scaling;
with Reduction_of_Polynomial_Systems;  use Reduction_of_Polynomial_Systems;  

with Homotopy;
with Start_Systems;                    use Start_Systems;
with BKK_Bound_Computations;           use BKK_Bound_Computations;

with Continuation_Parameters;
with Increment_and_Fix_Continuation;   use Increment_and_Fix_Continuation; 
with Root_Refiners;                    use Root_Refiners;

package body PHCPACK is

-- 1. PRE-PROCESSING : SCALING AND REDUCTION

  procedure Equation_Scaling
                ( file : in file_type; p : in Poly_Sys; s : out Poly_Sys ) is

    res : Poly_Sys(p'range);

  begin
    Copy(p,res);
    Scale(res);
    put(file,res);
    s := res;
  end Equation_Scaling;

  procedure Linear_Reduction
                ( file : in file_type; p : in Poly_Sys; r : out Poly_Sys ) is

    res : Poly_Sys(p'range);
    success,inconsistent,infinite : boolean := false;

  begin
    Copy(p,res);
    reduce(res,success,inconsistent,infinite);
    if success 
     then if inconsistent
           then put_line(file,"system is inconsistent");
          end if;
          if infinite
           then put_line(file,"system has infinite number of solutions");
          end if;
    end if;
    put(file,res);
    r := res;
  end Linear_Reduction;

-- 2. ROOT COUNTING AND CONSTRUCTION OF START SYSTEM

  procedure Total_Degree
                ( file : in file_type; p : in Poly_Sys; d : out natural ) is
  begin
    d := Total_Degree(p);
  end Total_Degree;

  procedure Total_Degree
                ( file : in file_type; p : in Poly_Sys; d : out natural;
                  q : out Poly_Sys; qsols : out Solution_List ) is

    qq : Poly_Sys(p'range);
    qqsols : Solution_List;

  begin
    d := Total_Degree(p);
    Start_System(p,qq,qqsols);
    q := qq; qsols := qqsols;
  end Total_Degree;

  procedure Implicit_Lifting
                 ( file : in file_type; p : in Poly_Sys; mv : out natural ) is
  begin
    mv := BKK_by_Implicit_Lifting(p);
  end Implicit_Lifting;

  procedure Implicit_Lifting
                 ( file : in file_type; p : in Poly_Sys; mv : out natural;
                   q : out Poly_Sys; qsols : out Solution_List ) is

    qq : Poly_Sys(p'range) := Complex_Randomize1(p);
    qqsols : Solution_List := Solve_by_Implicit_Lifting(file,qq);

  begin
    mv := Length_Of(qqsols);
    Set_Continuation_Parameter(qqsols,CMPLX(0.0));
    q := qq; qsols := qqsols;
  end Implicit_Lifting;

  procedure Static_Lifting
                 ( file : in file_type; p : in Poly_Sys; mv : out natural ) is
  begin
    mv := BKK_by_Static_Lifting(file,p);
  end Static_Lifting;

  procedure Static_Lifting
                 ( file : in file_type; p : in Poly_Sys; mv : out natural;
                   q : out Poly_Sys; qsols : out Solution_List ) is

    qq : Poly_Sys(p'range) := Complex_Randomize1(p);
    qqsols : Solution_List := Solve_by_Static_Lifting(file,qq);

  begin
    mv := Length_Of(qqsols);
    Set_Continuation_Parameter(qqsols,CMPLX(0.0));
    q := qq; qsols := qqsols;
  end Static_Lifting;

-- 3. POLYNOMIAL CONTINUATION

  procedure Artificial_Parameter_Continuation
                 ( file : in file_type; p,q : in Poly_Sys;
                   sols : in out Solution_List;
                   k : in natural := 2;
                   a : in double_complex := CMPLX(1.0); 
                   target : in double_complex := CMPLX(1.0) ) is

    procedure Cont is 
      new Reporting_Continue(Norm1,Homotopy.Eval,Homotopy.Eval,Homotopy.Diff);

  begin
    Homotopy.Create(p,q,k,a);
    Continuation_Parameters.Tune(0);
    Cont(file,sols,false,target);
    Homotopy.Clear;
  end Artificial_Parameter_Continuation;

  procedure Natural_Parameter_Continuation
                 ( file : in file_type; h : in Poly_Sys; k : in natural;
                   t0,t1 : in double_complex; sols : in out Solution_List ) is
  begin
    null;
  end Natural_Parameter_Continuation;

-- 4. POST-PROCESSING : VALIDATION

  procedure Refine_Roots
                 ( file : in file_type; p : in Poly_Sys;
                   sols : in out Solution_List ) is

    epsxa,epsfa : constant double_float := 10.0**(-8);   -- defaults
    tolsing : constant double_float := 10.0**(-8);
    maxit : constant natural := 3;
    numit : natural := 0;

  begin
    Reporting_Root_Refiner(file,p,sols,epsxa,epsfa,tolsing,numit,maxit,false);
  end Refine_Roots;

end PHCPACK;
