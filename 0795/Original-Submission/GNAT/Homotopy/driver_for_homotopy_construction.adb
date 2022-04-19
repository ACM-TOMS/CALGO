with integer_io;                        use integer_io;
with Communications_with_User;          use Communications_with_User;
with Floating_Point_Numbers;            use Floating_Point_Numbers;
with Complex_Numbers_io,Numbers_io;     use Complex_Numbers_io,Numbers_io;
with Complex_Multivariate_Polynomials;  use Complex_Multivariate_Polynomials;

with Homotopy;
with Random_Number_Generators;

with Projective_Transformations;        use Projective_Transformations;
with Homogenization;                    use Homogenization;

procedure Driver_for_Homotopy_Construction
             ( file : in file_type; p,q : in out Poly_Sys; 
               qsols : in out Solution_List; target : out double_complex ) is

  procedure Info_on_Power is

    i : array(1..6) of string(1..65);

  begin
    i(1):="  The  homotopy  parameter  k  determines  the   power   of   the";
    i(2):="continuation  parameter  t.  Taking k>1 has as effect that at the";
    i(3):="beginning and at the end of the continuation, changes in t do not";
    i(4):="change  the homotopy as much as with a homotopy that is linear in";
    i(5):="t so that paths are better to follow.  The default value  k=2  is";
    i(6):="recommended.                                                     ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Info_on_Power;

  procedure Info_on_Constant is

    i : array(1..3) of string(1..65);

  begin
    i(1):="  The homotopy parameter a ensures the regularity of the solution";
    i(2):="paths, i.e.: by choosing a random complex number for a, all paths";
    i(3):="are regular and do not diverge for t<1.                          ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Info_on_Constant;

  procedure Info_on_Target is

    i : array(1..7) of string(1..65);

  begin
    i(1):="  The target value for the continuation parameter t is by default";
    i(2):="1.   To  create  stepping stones in the continuation stage, it is";
    i(3):="possible to let the continuation stop at t<1, for instance at t =";
    i(4):="0.9  or even at a complex value for t.  The solutions at t<1 will";
    i(5):="serve at start solutions to take up the continuation in  a  later";
    i(6):="stage.   In this stage, the same homotopy parameters k and a must";
    i(7):="be used.                                                         ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Info_on_Target;

  procedure Info_on_Exit is

    i : constant string 
      :="By typing 0, the current settings are used in the homotopy.";

  begin
    put_line(i);
  end Info_on_Exit;

  procedure Info_on_Projective_Transformation is

    i : array(1..5) of string(1..65);

  begin
    i(1):="  A projective transformation of the homotopy and start solutions";
    i(2):="makes  the equations in the polynomials homogeneous and adds an a";
    i(3):="random hyperplane.   The  vectors  of  the  start  solutions  are";
    i(4):="extended  with an additional unknown.  For solutions at infinity,";
    i(5):="this additional unknown equals zero.                             ";
    for j in i'range loop
      put_line(i(j));
    end loop;
  end Info_on_Projective_Transformation;

  procedure Read_Complex ( x : in out double_complex ) is

  -- DESCRIPTION :
  --   User friendly reading of a complex number.

    re,im : double_float;

  begin
    put("  Give a real number for real part : ");      Read_Double_Float(re);
    put("  Give a real number for imaginary part : "); Read_Double_Float(im);
    x := CMPLX(re,im);
  end Read_Complex;

  procedure Driver is

    ans : character;
    choice : string(1..2) := "  ";
    k : positive := 2;
    a : double_complex := Random_Number_Generators.Random1;
    t : double_complex := CMPLX(1.0);
    prt : boolean := false;

  begin
    new_line;
    put_line
       ("Homotopy is H(x,t) = a*(1-t)^k * Q(x) + t^k * P(x) = 0, t in [0,1],");
    put_line
       ("      with Q(x) = 0 a start system, and P(x) = 0 the target system.");
    loop
      new_line;
      put_line("MENU with settings for homotopy :");
      put("  relaxation constant k : "); put(k,2); new_line;
      put("  regularity constant a : "); put(a); new_line;
      put("  target value for t    : "); put(t); new_line;
      put("  projective transformation : ");
           if prt then put("yes"); else put("no"); end if; new_line;
      put("Type k, a, t, or p to change, preceded by i for info." 
          & "  Type 0 to exit : ");
      Ask_Alternative(choice,"katp0",'i');
      exit when choice(1) = '0';
      case choice(1) is
        when 'k' => 
          put("Give a value (1,2, or 3) for the relaxation constant k : ");
          Read_Positive(k);
        when 'a' =>
          put_line("Reading a complex value for the regularity constant a.");
          loop
            Read_Complex(a);
            exit when modulus(a) /= 0.0;
            put_line("The value 0 for a is not allowed! Try again.");
          end loop;
        when 't' =>
          put_line("Reading the (complex) target value for t.");
          Read_Complex(t);
        when 'p' =>
          put("Do you want a projective transformation? ");
          Ask_Yes_or_No(ans); prt := (ans ='y');
        when 'i' =>
          new_line;
          case choice(2) is
            when 'k' => Info_on_Power;
            when 'a' => Info_on_Constant;
            when 't' => Info_on_Target;
            when 'p' => Info_on_Projective_Transformation;
            when '0' => Info_on_Exit;
            when others => null;
          end case;
        when others => null;
      end case;
    end loop;

    new_line(file);
    put_line(file,"HOMOTOPY PARAMETERS :");
    put(file,"  k : "); put(file,k,2); new_line(file);
    put(file,"  a : "); put(file,a);   new_line(file);
    put(file,"  t : "); put(file,t);   new_line(file);
    if prt
     then put_line(file,"  projective transformation");
     else put_line(file,"  no projective transformation");
    end if;

    target := t;

    if prt
     then Projective_Transformation(p);
          Projective_Transformation(q);
          Projective_Transformation(qsols);
          declare
            pp,sysp : Poly_Sys(p'first..p'last+1);
          begin
            pp := Add_Random_Hyperplanes(p,1,true);
            sysp := Add_Standard_Hyperplanes(q,1);
            Homotopy.Create(pp,sysp,k,a);
            Clear(pp); Clear(sysp);
          end;
     else Homotopy.Create(p,q,k,a);
    end if;
  end Driver;

begin
  Driver; 
end Driver_for_Homotopy_Construction;
