with text_io,integer_io;              use text_io,integer_io;
with Communications_with_User;        use Communications_with_User;
with File_Operations;                 use File_Operations;
with Complex_Vectors;                 use Complex_Vectors;
with Complex_Vectors_io;              use Complex_Vectors_io;
with Complex_Numbers_io;              use Complex_Numbers_io;
with Complex_Polynomial_Systems;      use Complex_Polynomial_Systems;
with Complex_Polynomial_Systems_io;   use Complex_Polynomial_Systems_io;
with Solutions,Solutions_io;          use Solutions,Solutions_io;
with Scaling;
with Drivers_for_Scaling;             use Drivers_for_Scaling;

procedure mainscal ( infilename,outfilename : in string ) is

  lp : Link_to_Poly_Sys;
  infile : file_type;
  outfile : file_type;
  n,basis : natural;
  scalvec : Link_to_Vector;
  ans : character;
  sysonfile : boolean;

  procedure Read_System ( file : in out file_type; filename : in string ) is
  begin
    if filename /= ""
     then Open(file,in_file,filename);
          new_line;
          get(file,lp);
          n := lp'length;
          sysonfile := true;
     else sysonfile := false;
    end if;
  exception
    when others =>
      new_line;
      put("Could not open file with name "); put_line(filename);
      sysonfile := false; lp := null; return;
  end Read_System;

  procedure Separate_File ( p : in Poly_Sys ) is

    scafile : file_type;

  begin
    new_line;
    put("Do you want the scaled system on separate file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then put_line("Reading the name of the output file.");
          Read_Name_and_Create_File(scafile);
          put(scafile,lp.all);
          if basis /= 0
           then new_line(scafile);
                put_line(scafile,"SCALING COEFFICIENTS :");
                new_line(scafile);
                put(scafile,basis,1); new_line(scafile);
                for i in scalvec'range loop
                  put(scafile,scalvec(i)); new_line(scafile);
                end loop;
          end if;
          Close(scafile);
    end if;
  end Separate_File;

  procedure Rescale is
  
    sols : Solution_List;
    found : boolean;
    m : natural;

  begin
    if sysonfile                             -- scan for scaling coefficients
     then Scan_and_Skip(infile,"SCALING COEFFICIENTS",found);
          if found
           then get(infile,basis);
                scalvec := new vector(1..2*n);
                get(infile,scalvec.all);
          end if;
     else found := false;
    end if;
    if not found
     then put("Give the basis : "); get(basis);
          put("Give "); put(2*n,1); put_line(" complex scaling numbers : ");
          scalvec := new vector(1..2*n);
          for i in scalvec'range loop
            get(scalvec(i));
          end loop;
    end if;
    if sysonfile                                    -- scan for the solutions
     then Reset(infile);
          Scan_and_Skip(infile,"SOLUTIONS",found);
          if found
           then get(infile,sols);
          end if;
          Close(infile);
     else found := false;
    end if;
    if not found
     then put_line("Reading the name of the file for the solutions.");
          Read_Name_and_Open_File(infile);
          get(infile,sols);
          Close(infile);
    end if;
    put_line(outfile,"THE SCALING COEFFICIENTS : ");
    new_line(outfile);
    put(outfile,basis,1); new_line(outfile);
    for i in scalvec'range loop
      put(outfile,scalvec(i)); new_line(outfile);
    end loop;
    new_line(outfile);
    Scaling.Scale(basis,scalvec.all,sols);
    m := Length_Of(sols);
    if m > 0
     then put_line(outfile,"THE DE-SCALED SOLUTIONS : ");
          new_line(outfile);
          put(outfile,m,Head_Of(sols).n,sols);
    end if;
    Close(outfile);
  end Rescale;

  procedure Display_and_Dispatch_Menu
               ( file : in file_type; p : in out Poly_Sys ) is

  -- DESCRIPTION :
  --   Displays the menu and returns a choice, corresponding to one of the
  --   three available scaling procedures.

  begin
    loop
      new_line;
      put_line("MENU for Scaling Polynomial Systems :");
      put_line("  1 : Equation Scaling : divide by average coefficient      ");
      put_line("  2 : Variable Scaling : change of variables, as z = (2^c)*x");
      put_line("  3 : Solution Scaling : back to original coordinates       ");
      put("Type 1, 2, or 3 to select scaling, or i for info : ");
      Ask_Alternative(ans,"123i");
      if ans = 'i'
       then new_line; Drivers_for_Scaling.Display_Info; new_line;
      end if;
      exit when ans /= 'i';
    end loop;
    case ans is
      when '1' => Equation_Scaling(file,p); basis := 0;
      when '2' => Variable_Scaling(file,p,basis,scalvec);
      when '3' => Rescale;
      when others => null;
    end case;
    case ans is
      when '1' | '2' => Write_Results(file,p,basis,scalvec);
      when others    => null;
    end case;
    if ans /= '3'
     then Separate_File(p);
    end if;
  end Display_and_Dispatch_Menu;

begin
  Read_System(infile,infilename);
  if lp = null
   then loop
          new_line;
          put("Is the system on a file ? (y/n/i=info) ");
          Ask_Alternative(ans,"yni");
          if ans = 'i'
           then new_line;
                Complex_Polynomial_Systems_io.Display_Format;
                new_line;
          end if;
          exit when ans /= 'i';
        end loop;
        new_line;
        if ans = 'y'
         then put_line("Reading the name of the input file.");
              Read_Name_and_Open_File(infile);
              get(infile,lp);
              sysonfile := true;
              n := lp'length;
         else put("Give the dimension : "); get(n);
              lp := new Poly_Sys(1..n);
              put("Give "); put(n,1); put(" "); put(n,1);
              put_line("-variate polynomials :");
              get(n,lp.all);
              skip_line;  -- skip end_of_line symbol
              sysonfile := false;
        end if;
  end if;
  Create_Output_File(outfile,outfilename);
  put(outfile,lp.all); new_line(outfile);
  Display_and_Dispatch_Menu(outfile,lp.all);
end mainscal;
