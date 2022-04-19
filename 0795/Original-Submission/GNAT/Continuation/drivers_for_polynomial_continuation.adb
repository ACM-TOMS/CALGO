with integer_io,Numbers_io;               use integer_io,Numbers_io;
with Communications_with_User;            use Communications_with_User;
with Timing_Package;                      use Timing_Package;
with File_Operations;                     use File_Operations;

with Floating_Point_Numbers;              use Floating_Point_Numbers;
with Complex_Norms,Complex_Numbers_io;    use Complex_Norms,Complex_Numbers_io;
with Complex_Multivariate_Polynomials;
with Symbol_Table,Symbol_Table_io;        use Symbol_Table;
with Complex_Polynomial_Systems_io;       use Complex_Polynomial_Systems_io;

with Solutions_io;                        use Solutions_io;
with Continuation_Parameters;
with Continuation_Parameters_io;
with Homotopy;
with Projective_Transformations;          use Projective_Transformations;
with Increment_and_Fix_Continuation;      use Increment_and_Fix_Continuation;

with Process_io;                          use Process_io;
with Driver_for_Homotopy_Construction;
with Drivers_for_Path_Directions;         use Drivers_for_Path_Directions;

package body Drivers_for_Polynomial_Continuation is

-- AUXILIARIES :

  procedure Continue ( file : in file_type; sols : in out Solution_List;
                       proj,report : in boolean;
                       target : in double_complex ) is

  -- DECRIPTION :
  --   Performs a normal continuation, without toric compactifications.
   
    timer : timing_widget;

    procedure Sil_Cont is
      new Silent_Continue(Norm1,Homotopy.Eval,Homotopy.Diff,Homotopy.Diff);
    procedure Rep_Cont is
      new Reporting_Continue(Norm1,Homotopy.Eval,Homotopy.Diff,Homotopy.Diff);

  begin
    tstart(timer);
    if report
     then Rep_Cont(file,sols,proj,target);
     else Sil_Cont(sols,proj,target);
    end if;
    tstop(timer);
    new_line(file); print_times(file,timer,"continuation");
  end Continue;

  procedure Ask_Symbol is

  -- DESCRIPTION :
  --   This procedure asks for the symbol to display the additional unknown.

    sb : Symbol;

  begin
    put("Give symbol to display additional unknown : ");
    sb := (sb'range => ' ');
    Symbol_Table.Enlarge(1);
    Symbol_Table_io.Get(sb);
    Symbol_Table.Add(sb);
  end Ask_Symbol;

-- TARGET ROUTINES :

  procedure Driver_for_Process_io ( file : in file_type; oc : out natural ) is

    ans : character;
    m : array(0..8) of string(1..65);

  begin
    put_line("MENU for Output Information during Continuation : ");
    m(0):="  0 : no intermediate output information during continuation     ";
    m(1):="  1 : only the final solutions at the end of the paths           ";
    m(2):="  2 : intermediate solutions at each step along the paths        ";
    m(3):="  3 : information of the predictor: t and step length            ";
    m(4):="  4 : information of the corrector: corrections and residuals    ";
    m(5):="  5 : intermediate solutions and information of the predictor    ";
    m(6):="  6 : intermediate solutions and information of the corrector    ";
    m(7):="  7 : information of predictor and corrector                     ";
    m(8):="  8 : intermediate solutions, info of predictor and corrector    ";
    for i in m'range loop
      put_line(m(i));
    end loop;
    put("Type a number beween 0 and 8 to select output information : ");
    Ask_Alternative(ans,"012345678");
    new_line(file);
    put_line(file,"OUTPUT INFORMATION DURING CONTINUATION :");
    case ans is
      when '0' => Set_output_code(nil); oc := 0; put_line(file,m(0));
      when '1' => Set_output_code(nil); oc := 1; put_line(file,m(1));
      when '2' => Set_output_code(s);   oc := 2; put_line(file,m(2));
      when '3' => Set_output_code(p);   oc := 3; put_line(file,m(3));
      when '4' => Set_output_code(c);   oc := 4; put_line(file,m(4));
      when '5' => Set_output_code(sp);  oc := 5; put_line(file,m(5));
      when '6' => Set_output_code(sc);  oc := 6; put_line(file,m(6));
      when '7' => Set_output_code(pc);  oc := 7; put_line(file,m(7));
      when '8' => Set_output_code(spc); oc := 8; put_line(file,m(8));
      when others => null;
    end case;
  end Driver_for_Process_io;

  procedure Driver_for_Continuation_Parameters ( file : in file_type ) is

    nb : natural := 0;

    procedure Begin_Banner ( file : in file_type ) is
    begin
      put_line(file,"****************** CURRENT CONTINUATION PARAMETERS "
                     & "*****************");
    end Begin_Banner;

    procedure End_Banner ( file : in file_type ) is
    begin
      put_line(file,"***************************************************"
                     & "*****************");
    end End_Banner;

  begin
    loop
      Begin_Banner(Standard_Output);
      Continuation_Parameters_io.put;
      End_Banner(Standard_Output);
      Continuation_Parameters_io.get(nb);
      exit when (nb = 0);
    end loop;
    new_line(file);
    Begin_Banner(file);
    Continuation_Parameters_io.put(file);
    End_Banner(file);
  end Driver_for_Continuation_Parameters;

  procedure Check_Continuation_Parameter
                ( sols : in out Solution_List ) is

    ans : character;
    tre,tim : double_float;

  begin
    if not Is_Null(sols)
     then
       if Head_Of(sols).t = CMPLX(1.0)
        then put_line("The first solution has continuation parameter t = 1.0.");
             put("Do you want to change t ? (y/n) "); Ask_Yes_or_No(ans);
             if ans = 'y'
              then put("Give real part of t : "); Read_Double_Float(tre);
                   put("Give imaginary part of t : "); Read_Double_Float(tim);
                   Set_Continuation_Parameter(sols,CMPLX(tre,tim));
             end if;
       end if;
    end if;
  end Check_Continuation_Parameter;

  procedure Driver_for_Polynomial_Continuation 
                ( file : in file_type; p : in Poly_Sys; 
                  sols : out Solution_List; target : out double_complex ) is

    infile : file_type;
    pp,q : Poly_Sys(p'range);
    t : double_complex;
    qsols : Solution_List;
    found,proj : boolean;

    use Complex_Multivariate_Polynomials;

    procedure Read_Start_System is
    begin
      put_line("Reading the name of the file for start system.");
      Read_Name_and_Open_File(infile);
      get(infile,q);
    exception
      when others => put("The system on the file is not correct.");
                     put_line("  Try again..."); Close(infile);
                     Read_Start_System;
    end Read_Start_System;

  begin
    new_line;
    Read_Start_System;
    put_line(file,"THE START SYSTEM : "); put(file,q); new_line(file);
    Scan_and_Skip(infile,"SOLUTIONS",found);
    if found 
     then get(infile,qsols);
     else new_line; Read(qsols);
    end if;
    Close(infile);
    Check_Continuation_Parameter(qsols);
    put_line(file,"THE START SOLUTIONS : ");
    put(file,Length_Of(qsols),Head_Of(qsols).n,qsols); new_line(file);
    Copy(p,pp);
    Driver_for_Homotopy_Construction(file,pp,q,qsols,t);
    proj := (Number_of_Unknowns(q(q'first)) > q'last);
    if proj
     then Ask_Symbol;
    end if;
    new_line;
    Driver_for_Polynomial_Continuation(file,qsols,proj,t);
   -- Homotopy.Clear;  --> clearing here creates difficulties for root refiner
    sols := qsols;
    target := t;
  end Driver_for_Polynomial_Continuation;

  procedure Driver_for_Polynomial_Continuation
                ( file : in file_type; p : in Poly_Sys; k : in natural;
                  target : in double_complex; sols : out Solution_list ) is

    qsols : Solution_List;

  begin
    new_line; Read(qsols);
    put_line(file,"THE START SOLUTIONS :");
    put(file,Length_Of(qsols),Head_Of(qsols).n,qsols); new_line(file);
    Homotopy.Create(p,k);
    put_line(file,"HOMOTOPY PARAMETERS :");
    put(file,"  k : "); put(file,k,2); new_line(file);
    put(file,"  a : "); put(file,target); new_line(file);
    Driver_for_Polynomial_Continuation(file,qsols,false,target);
   -- Homotopy.Clear; --> clearing here creates difficulties for root refiner
    sols := qsols;
  end Driver_for_Polynomial_Continuation;

  procedure Driver_for_Polynomial_Continuation
                ( file : in file_type; sols : in out Solution_List;
                  proj : in boolean;
                  target : double_complex := CMPLX(1.0) ) is

    oc : natural;
    timer : timing_widget;
    report : boolean;
    n : constant natural := Head_Of(sols).n;
    nv : constant natural := Length_Of(sols);
    v : Float_Vectors_of_Vectors.Link_to_Vector;
    errv : Float_Vectors.Link_to_Vector;

  begin
    new_line;
    Driver_for_Continuation_Parameters(file);
    if Continuation_Parameters.endext_order > 0
     then Init_Path_Directions(n,nv,v,errv);
    end if;
    new_line;
    Driver_for_Process_io(file,oc);
    report := (oc /= 0);
    new_line;
    put_line("No more input expected.  See output file for results.");
    new_line;
    if Continuation_Parameters.endext_order > 0
     then Toric_Continue(file,sols,proj,report,v.all,errv.all,target);
          Write_Directions(file,v.all,errv.all);
     else Continue(file,sols,proj,report,target);
    end if;
  end Driver_for_Polynomial_Continuation;

end Drivers_for_Polynomial_Continuation;
