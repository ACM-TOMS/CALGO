with text_io,integer_io;                use text_io,integer_io;
with Communications_with_User;          use Communications_with_User;
with Timing_Package;                    use Timing_Package;
with Float_Vectors_of_Vectors;
with Complex_Numbers,Complex_Vectors;   use Complex_Numbers,Complex_Vectors;
with Complex_Polynomial_Systems;        use Complex_Polynomial_Systems;
with Complex_Polynomial_Systems_io;     use Complex_Polynomial_Systems_io;
with Complex_Multivariate_Polynomials;
with Homotopy,Solutions;                use Solutions;

with Drivers_for_Scaling;               use Drivers_for_Scaling;
with Drivers_for_Reduction;             use Drivers_for_Reduction;
with Driver_for_Root_Counts;
with Driver_for_Homotopy_Construction;
with Driver_for_Root_Refining;
with Drivers_for_Polynomial_Continuation;
 use Drivers_for_Polynomial_Continuation;

with Bye_Bye_Message;

procedure mainphc ( infilename,outfilename : in string ) is

  outpt : file_type;
  lp : Link_to_Poly_Sys := null;
  timer : timing_widget;

  procedure Display_Options is

  -- DESCRIPTION :
  --   Displays an overview of all options on screen.

    o : array(1..7) of string(1..65);

  begin
    put_line("Running in full mode.  Note also the following options:");
    o(1) := "  phc -s : Equation and variable Scaling on system and solutions ";
    o(2) := "  phc -d : Linear and nonlinear Reduction w.r.t. the total degree";
    o(3) := "  phc -r : Root counting and Construction of start systems       ";
    o(4) := "  phc -m : Mixed-Volume Computation by four lifting strategies   ";
    o(5) := "  phc -p : Polynomial Continuation by a homotopy in one parameter";
    o(6) := "  phc -v : Validation, refinement and purification of solutions  ";
    o(7) := "  phc -b : Batch or black-box processing                         ";
    for i in o'range loop
      put_line(o(i));
    end loop;
  end Display_Options;

  procedure Read_System ( filename : in string ) is

    file : file_type;
    n : natural;

  begin
    if filename /= ""
     then Open_Input_File(file,filename);
          get(file,n);
          lp := new Poly_Sys(1..n);
          get(file,n,lp.all);
          Close(file);
    end if;
  exception
    when others =>
      new_line;
      put("Could not open file with name "); put_line(filename);
      lp := null; return;
  end Read_System;

begin
  new_line;
  Display_Options;

-- READ THE INPUT SYSTEM AND OUTPUT FILE

  Read_System(infilename);
  if lp = null
   then new_line; get(lp);
  end if;

  Create_Output_File(outpt,outfilename);
  put(outpt,lp.all);

  tstart(timer);

  declare

    p,q,scalp,projp : Poly_Sys(lp'range);
    target : double_complex;
    basis,roco : natural;
    scalvec : Link_to_Vector;
    sols : Solution_List;
    proj : boolean;
    use Complex_Multivariate_Polynomials;

  begin

    p := lp.all;  Copy(p,scalp);

   -- PREPROCESSING : SCALING AND REDUCTION

    Driver_for_Scaling(outpt,scalp,basis,scalvec);
    Driver_for_Reduction(outpt,scalp,roco,true);

   -- APPLY ROOT COUNTING METHODS TO CONSTRUCT A START SYSTEM

    Copy(scalp,projp);
    Driver_for_Root_Counts(outpt,projp,q,true,sols,roco);

    if Length_Of(sols) > 0
     then

      -- CONSTRUCTION OF THE HOMOTOPY

       Driver_for_Homotopy_Construction(outpt,projp,q,sols,target);

      -- CONTINUATION

       proj := (Number_of_Unknowns(p(p'first)) > p'last);
       if Head_Of(sols).t /= CMPLX(0.0)
        then Set_Continuation_Parameter(sols,CMPLX(0.0));
       end if;
       Driver_for_Polynomial_Continuation(outpt,sols,proj,target);

      -- ROOT REFINING

       Driver_for_Root_Refining(outpt,scalp,p,basis,scalvec,sols);
     end if;

  end;

  tstop(timer);
  new_line(outpt);
  print_times(outpt,timer,"solving the polynomial system");
  new_line(outpt);
  put(outpt,Bye_Bye_Message);
  Close(outpt);

end mainphc;
