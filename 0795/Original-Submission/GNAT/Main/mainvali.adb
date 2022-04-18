with text_io,integer_io;                use text_io,integer_io;
with Floating_Point_Numbers;            use Floating_Point_Numbers;
with Timing_Package;                    use Timing_Package;
with Communications_with_User;          use Communications_with_User;
with File_Operations,Numbers_io;        use File_Operations,Numbers_io;
with Complex_Numbers,Complex_Vectors;   use Complex_Numbers,Complex_Vectors;
with Complex_Vectors_io;                use Complex_Vectors_io;
with Complex_Numbers_io;                use Complex_Numbers_io;
with Complex_Polynomial_Systems;        use Complex_Polynomial_Systems;
with Complex_Polynomial_Systems_io;     use Complex_Polynomial_Systems_io;
with Solutions,Solutions_io;            use Solutions,Solutions_io;

with Symmetry_Group;                    use Symmetry_Group;
with Symbolic_Symmetry_Group_io;        use Symbolic_Symmetry_Group_io;
with Drivers_for_Symmetry_Group_io;     use Drivers_for_Symmetry_Group_io;
with Drivers_for_Orbits_of_Solutions;   use Drivers_for_Orbits_of_Solutions;

with Root_Refiners;                     use Root_Refiners;
with Driver_for_Winding_Numbers;
with valipoco;

with Bye_Bye_Message;

procedure mainvali ( infilename,outfilename : in string ) is

  use Floating_Point_Numbers.double_float_io;

  procedure Display_Validation_Info is

  -- DESCRIPTION :
  --   Displays information about available validation methods on screen.

    i : array(1..9) of string(1..65);

  begin
    i(1):="Basic validation consists in the application of  Newton's  method";
    i(2):="on  the  list  of solutions.  There are facilities to extract the";
    i(3):="generating solutions when the symmetry group is submitted.       ";
    i(4):="  Winding  numbers  can  be  computed  by  homotopy  continuation";
    i(5):="methods.   The user must provide a start system with solutions at";
    i(6):="t < 1.                                                           ";
    i(7):="  Polyhedral validation is based on  the  output  file  of  poco,";
    i(8):="where  the  polyhedral  end  game was turned on.  This validation";
    i(9):="puts up a frequency table of computed path directions.           ";
    for k in i'range loop
      put_line(i(k));
    end loop;
  end Display_Validation_Info;

-- READING THE INPUT :

  procedure Scan_System ( file : in out file_type; filename : in string;
                          lp : in out Link_to_Poly_Sys;
                          sysonfile : out boolean ) is

  -- DESCRIPTION :
  --   Checks whether the given file name corresponds to a file with
  --   a polynomial system in a correct format.
  --   If this is the case, then sysonfile is true on return and lp
  --   contains the system.

  begin
    if filename /= ""
     then Open(file,in_file,filename);
          get(file,lp);
          sysonfile := true;
     else sysonfile := false;
    end if;
  exception
    when others =>
      new_line;
      put("Could not open file with name "); put_line(filename);
      lp := null;
      sysonfile := false;
      return;
  end Scan_System;

  procedure Read_System ( file : in out file_type; filename : in string;
                          lp : in out Link_to_Poly_Sys;
                          sysonfile : out boolean ) is

  -- DESCRIPTION :
  --   Searches first the system on file, using the given filename.
  --   If necessary other files will be openend.

    ans : character;
    n : natural;

  begin
    Scan_System(file,filename,lp,sysonfile);
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
                Read_Name_and_Open_File(file);
                get(file,lp);
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
  end Read_System;

  procedure Scan_Solutions
               ( file : in out file_type; sysonfile : in boolean;
                 sols : in out Solution_List; found : out boolean ) is

    fnd : boolean := false;

  begin
    if sysonfile
     then Scan_and_Skip(file,"SOLUTIONS",fnd);
          if fnd
           then get(file,sols);
          end if;
          Close(file);
     else fnd := false;
    end if;
    found := fnd;
  exception
    when others
      => put_line("Something is wrong with the solutions, will ignore...");
         Close(file);
         found := false;
  end Scan_Solutions;

  procedure Read_Solutions
              ( file : in out file_type; sysonfile : in boolean;
                sols : in out Solution_List ) is

    found : boolean;

  begin
    Scan_Solutions(file,sysonfile,sols,found);
    if not found
     then new_line;
          put_line("Reading the name of the file for the solutions.");
          Read_Name_and_Open_File(file);
          get(file,sols);
          Close(file);
    end if;
  end Read_Solutions;

-- ROOT REFINING AUXILIARIES :

  procedure Default_Root_Refining_Parameters
              ( epsxa,epsfa,tolsing : out double_float;
                maxit : out natural; wout : out boolean ) is

  -- DESCRIPTION :
   --  Defines the default values for the root refining parameters.

  begin
    epsxa := 10.0**(-8);    -- precision for correction on x
    epsfa := 10.0**(-8);    -- precision for residual 
    tolsing := 10.0**(-8);  -- tolerance on inverse condition numbers
    maxit := 3;             -- maximal number of Newton iterations
    wout := false;          -- if intermediate output is wanted
  end Default_Root_Refining_Parameters;

  procedure Put_Root_Refining_Parameters
              ( file : in file_type; epsxa,epsfa,tolsing : in double_float;
                maxit : in natural; wout : in boolean ) is

  -- DESCRIPTION :
  --   Writes the parameters for the root refiner on file.

  begin
    put(file,"  1. output during the iterations    : ");
    if wout
     then put(file," yes"); new_line(file);
     else put(file," no"); new_line(file);
    end if;
    put(file,"  2. tolerance for error on the root : ");
    put(file,epsxa,2,3,3); new_line(file);
    put(file,"  3. tolerance for the residual      : ");
    put(file,epsfa,2,3,3); new_line(file);
    put(file,"  4. tolerance for singular roots    : ");
    put(file,tolsing,2,3,3); new_line(file);
    put(file,"  5. maximum number of iterations    : ");
    put(file,maxit,2); new_line(file);
  end Put_Root_Refining_Parameters;

  procedure Menu_Root_Refining_Parameters
              ( file : in file_type; epsxa,epsfa,tolsing : in out double_float;
                maxit : in out natural; wout : in out boolean ) is

  -- DESCRIPTION :
  --   The user can set the parameters of the root refiner by following
  --   the menu's.

    ans : character;

  begin
    new_line;
    loop
      put_line("MENU with current Settings for the Root Refiner :");
      Put_Root_Refining_Parameters(standard_output,
                                   epsxa,epsfa,tolsing,maxit,wout);
      put("Type 1,2,3,4, or 5 to change, type 0 to exit : ");
      Ask_Alternative(ans,"012345");
      exit when ans = '0';
      case ans is
        when '1' => put("Do you want output during the iterations ? (y/n) ");
                    Ask_Yes_or_No(ans); wout := (ans = 'y');
        when '2' => put("Give new tolerance for error on the root : ");
                    Read_Double_Float(epsxa);
        when '3' => put("Give new tolerance for residual : ");
                    Read_Double_Float(epsfa);
        when '4' => put("Give new tolerance for singular roots : ");
                    Read_Double_Float(tolsing);
        when '5' => put("Give new maximum number of iterations : ");
                    Read_Natural(maxit);
        when others => null;
      end case;
    end loop;
    new_line(file);
    put_line(file,"ROOT REFINING PARAMETERS : ");
    Put_Root_Refining_Parameters(file,epsxa,epsfa,tolsing,maxit,wout);
  end Menu_Root_Refining_Parameters;

  procedure Refine_Roots
                 ( file : in file_type; p : in Poly_Sys;
                   solsfile,invar,allperms,signsym : in boolean;
                   v : in List_of_Permutations;
                   epsxa,epsfa,tolsing : in double_float;
                   maxit : in natural; wout : in boolean;
                   sols,refsols : in out Solution_List ) is

  -- DESCRIPTION :
  --   Refines the roots and computes generating solutions when required.

  -- ON ENTRY :
  --   file        for writing results on;
  --   p           the polynomial system under consideration;
  --   solsfile    whether refined solution have to go to separate file;
  --   invar       whether generating solutions have to be computed;
  --   allperms    whether invariant under all permutations;
  --   signsym     whether there is sign-symmetry;
  --   v           group representation, only needed when invar;
  --   sols        solutions that need to be refined.

  -- ON RETURN :
  --   sols        solutions after applying some Newton iteration;
  --   refsols     refined solutions, with the exception of failures and
  --               the non-generating solutions.

    numit : natural := 0;

  begin
    if solsfile or invar
     then Reporting_Root_Refiner
            (file,p,sols,refsols,epsxa,epsfa,tolsing,numit,maxit,wout);
          if invar
           then Driver_for_Orbits_of_Solutions
                  (file,refsols,v,allperms,signsym,epsxa);
          end if;
     else Reporting_Root_Refiner
            (file,p,sols,epsxa,epsfa,tolsing,numit,maxit,wout);
    end if;
  end Refine_Roots;

  procedure Refine_Roots
                 ( file : in file_type; p : in Poly_Sys;
                   solsfile : in boolean;
                   epsxa,epsfa,tolsing : in double_float;
                   maxit : in natural; wout : in boolean;
                   sols,refsols : in out Solution_List ) is

  -- DESCRIPTION : 
  --   Root refinement without computing of generating solutions.

    numit : natural := 0;

  begin
    if solsfile
     then Reporting_Root_Refiner
            (file,p,sols,refsols,epsxa,epsfa,tolsing,numit,maxit,wout);
     else Reporting_Root_Refiner
            (file,p,sols,epsxa,epsfa,tolsing,numit,maxit,wout);
    end if;
  end Refine_Roots;

  procedure End_of_Input_Message is
  begin
    new_line;
    put_line("No more input expected.  See output file for results.");
    new_line;
  end End_of_Input_Message;

-- VALIDATION PROCEDURES :

  procedure Winding_Validation is

  -- DESCRIPTION :
  --   Validation by computing winding numbers by homotopy continuation.

    lp : Link_to_Poly_Sys;
    timer : timing_widget;
    infile,solsft,outfile : file_type;
    ans : character;
    sysonfile,solsfile,wout : boolean;
    sols,refsols: Solution_List;
    epsxa,epsfa,tolsing : double_float;
    maxit : natural;

  begin
    Read_System(infile,infilename,lp,sysonfile);
    Create_Output_File(outfile,outfilename);
    put(outfile,lp.all);
    Read_Solutions(infile,sysonfile,sols);
    new_line;
    put("Do you want the refined solutions on separate file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then solsfile := true;
          put_line("Reading the name of the file to write the solutions on.");
          Read_Name_and_Create_File(solsft);
     else solsfile := false;
    end if;
    Default_Root_Refining_Parameters(epsxa,epsfa,tolsing,maxit,wout);
    Menu_Root_Refining_Parameters(outfile,epsxa,epsfa,tolsing,maxit,wout);
    Driver_for_Winding_Numbers(outfile,lp.all,sols);
    tstart(timer);
    Refine_Roots(outfile,lp.all,solsfile,
                 epsxa,epsfa,tolsing,maxit,wout,sols,refsols);
    tstop(timer);
    if solsfile
     then put(solsft,Length_Of(refsols),Head_Of(refsols).n,refsols);
          Close(solsft);
    end if;
    new_line(outfile);
    print_times(outfile,timer,"Root Refinement");
    Close(outfile);
  end Winding_Validation;

  procedure Weeding_Validation is

  -- DESCRIPTION :
  --   Validation by refining the roots and weeding out the solution set.

    lp : Link_to_Poly_Sys;
    timer : timing_widget;
    infile,solsft,outfile : file_type;
    n,maxit : natural;
    ans : character;
    sysonfile,solsfile,wout : boolean;
    invar,allperms,signsym,allsigns : boolean;
    g,v : List_of_Permutations;
    sols,refsols: Solution_List;
    epsxa,epsfa,tolsing : double_float;

  begin
    Read_System(infile,infilename,lp,sysonfile);
    Create_Output_File(outfile,outfilename);
    put(outfile,lp.all);
    Read_Solutions(infile,sysonfile,sols);
    new_line;
    put("Is the system invariant under group actions ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then invar := true; n := lp'length;
          Read_Symmetry_Group(n,g,v,allperms,signsym,allsigns);
          new_line(outfile);
          put_line(outfile,"THE SYMMETRY GROUP : ");
          new_line(outfile);
          Symbolic_Symmetry_Group_io.put(outfile,v);
          new_line(outfile);
     else invar := false;
    end if;
    new_line;
    put("Do you want the refined solutions on separate file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then solsfile := true;
          put_line("Reading the name of the file to write the solutions on.");
          Read_Name_and_Create_File(solsft);
     else solsfile := false;
    end if;
    Default_Root_Refining_Parameters(epsxa,epsfa,tolsing,maxit,wout);
    Menu_Root_Refining_Parameters(outfile,epsxa,epsfa,tolsing,maxit,wout);
    End_of_Input_Message;
    tstart(timer);
    Refine_Roots(outfile,lp.all,solsfile,invar,allperms,signsym,v,
                 epsxa,epsfa,tolsing,maxit,wout,sols,refsols);
    tstop(timer);
    if solsfile
     then put(solsft,Length_Of(refsols),Head_Of(refsols).n,refsols);
          Close(solsft);
    end if;
    new_line(outfile);
    print_times(outfile,timer,"Root Refinement");
    Close(outfile);
  end Weeding_Validation;

  procedure Polyhedral_End_Game_Validation is

  -- DESCRIPTION :
  --   Validation of the polyhedral end game.

    pocofile,resultfile : file_type;

  begin
    new_line;
    put_line("Reading name of the output file of poco.");
    Read_Name_and_Open_File(pocofile);
    new_line;
    put_line("Reading name of output file.");
    Read_Name_and_Create_File(resultfile);
    End_of_Input_Message;
    valipoco(pocofile,resultfile);
    Close(pocofile);
    new_line(resultfile);
    put(resultfile,Bye_Bye_Message);
    Close(resultfile);
  end Polyhedral_End_Game_Validation;

  procedure Display_and_Dispatch_Menu is

    ans : character;
    timer : timing_widget;

  begin
    loop
      new_line;
      put_line("MENU with Validation Methods : ");
      put_line
         ("  1. Basic Validation : refining and weeding out the solution set");
      put_line
         ("  2. Winding-Number Computation by homotopy continuation");
      put_line
         ("  3. Polyhedral Validation : frequency table of path directions");
      put("Type 1, 2, or 3 to select validation method, or i for info : ");
      Ask_Alternative(ans,"123i");
      case ans is
        when 'i' => new_line;
                    Display_Validation_Info;
        when '1' => Weeding_Validation;
        when '2' => Winding_Validation;
        when '3' => Polyhedral_End_Game_Validation;
        when others => null;
      end case;
      exit when ans /= 'i';
    end loop;
  end Display_and_Dispatch_Menu;

begin
  Display_and_Dispatch_Menu;
end mainvali;
