with text_io,integer_io;              use text_io,integer_io;
with Floating_Point_Numbers;          use Floating_Point_Numbers;
with Timing_Package;                  use Timing_Package;
with Communications_with_User;        use Communications_with_User;
with File_Operations,Numbers_io;      use File_Operations,Numbers_io;
with Complex_Vectors;                 use Complex_Vectors;
with Complex_Vectors_io;              use Complex_Vectors_io;
with Complex_Numbers_io;              use Complex_Numbers_io;
with Complex_Polynomial_Systems;      use Complex_Polynomial_Systems;
with Complex_Polynomial_Systems_io;   use Complex_Polynomial_Systems_io;
with Solutions,Solutions_io;          use Solutions,Solutions_io;
with Root_Refiners;                   use Root_Refiners;

procedure bablvali ( infilename,outfilename : in string ) is

  use Floating_Point_Numbers.double_float_io;

  timer : timing_widget;
  lp : Link_to_Poly_Sys;
  infile,outfile : file_type;
  ans : character;
  n : natural;
  sysonfile,found : boolean;
  sols : Solution_List;

  procedure Read_System ( file : in out file_type; filename : in string ) is
  begin
    if filename /= ""
     then Open(file,in_file,filename);
          get(file,lp);
          sysonfile := true;
     else sysonfile := false;
    end if;
  exception
    when others => put_line("Something is wrong with argument file...");
                   sysonfile := false;
                   lp := null; return;
  end Read_System;

  procedure Refine_Roots is

    epsxa,epsfa : constant double_float := 10.0**(-8);
    tolsing : constant double_float := 10.0**(-8);
    maxit : constant natural := 3;
    numb : natural := 0;
    refsols : Solution_List;

  begin
    new_line(outfile);
    put_line(outfile,"ROOT REFINING PARAMETERS");
    put(outfile,"  tolerance for error on the root : ");
    put(outfile,epsxa,2,3,3); new_line(outfile);
    put(outfile,"  tolerance for residual          : ");
    put(outfile,epsfa,2,3,3); new_line(outfile);
    put(outfile,"  tolerance for singular roots    : ");
    put(outfile,tolsing,2,3,3); new_line(outfile);
    put(outfile,"  maximum number of iterations    : ");
    put(outfile,maxit,2); new_line(outfile);
    tstart(timer);
    Reporting_Root_Refiner
      (outfile,lp.all,sols,refsols,epsxa,epsfa,tolsing,numb,maxit,false);
    tstop(timer);
    new_line(outfile);
    put_line(outfile,"THE SOLUTIONS :");
    put(outfile,Length_Of(refsols),Head_Of(refsols).n,refsols);
    new_line(outfile);
    print_times(outfile,timer,"Root refining");
  end Refine_Roots;

begin
  Read_System(infile,infilename);
  if lp = null
   then 
     new_line;
     put("Is the system on file ? (y/n) ");
     Ask_Yes_or_No(ans);
     if ans = 'y'
      then put_line("Reading the name of the input file.");
           Read_Name_and_Open_File(infile);
           get(infile,lp);
           sysonfile := true;
      else put("Give the dimension : "); get(n);
           lp := new Poly_Sys(1..n);
           put("Give "); put(n,1); put(" "); put(n,1); 
           put_line("-variate polynomials :");
           get(n,lp.all);
           skip_line;
           sysonfile := false;
     end if;
  end if;

  Create_Output_File(outfile,outfilename);
  put(outfile,lp.all);

  if sysonfile
   then Scan_and_Skip(infile,"SOLUTIONS",found);
        if found
         then get(infile,sols);
        end if;
        Close(infile);
   else found := false;
  end if;
  if not found
   then new_line; Read(sols);
  end if;

  Refine_Roots;
end bablvali;
