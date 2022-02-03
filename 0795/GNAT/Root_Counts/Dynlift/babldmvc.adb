with text_io,integer_io;                use text_io,integer_io;
with Communications_with_User;          use Communications_with_User;
with Timing_Package;                    use Timing_Package;
with Complex_Polynomial_Systems;        use Complex_Polynomial_Systems;
with Complex_Polynomial_Systems_io;     use Complex_Polynomial_Systems_io;
with Solutions,Solutions_io;            use Solutions,Solutions_io;
with Integer_Vectors;                   use Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;    use Arrays_of_Integer_Vector_Lists;
with Integer_Mixed_Subdivisions;        use Integer_Mixed_Subdivisions;
with Black_Box_Mixed_Volume_Computations;
 use Black_Box_Mixed_Volume_Computations;

procedure babldmvc ( infilename,outfilename : in string ) is

  infile,outfile : file_type;
  lp : Link_to_Poly_Sys;

  procedure Read_System ( file : in out file_type; filename : in string ) is
  begin
    if filename /= ""
     then Open_Input_File(file,filename);
          get(file,lp);
    end if;
  exception
    when others => put_line("Something is wrong with argument file...");
                   lp := null; return;
  end Read_System;

  procedure Main ( file : in file_type; p : in Poly_Sys ) is

    timer : Timing_Widget;
    q : Poly_Sys(p'range);
    qsols : Solution_List;
    mix : Link_to_Vector;
    lifsup : Link_to_Array_of_Lists;
    mixsub : Mixed_Subdivision;
    mv : natural;

  begin
    tstart(timer);
    Black_Box_Mixed_Volume_Computation(p,mix,lifsup,mixsub,mv);
    tstop(timer);
    new_line(outfile);
    put(outfile,"mixed volume : "); put(outfile,mv,1); new_line(outfile);
    new_line(outfile);
    print_times(outfile,timer,"Mixed-Volume Computation");
    if mv > 0
     then tstart(timer);
          Black_Box_Polyhedral_Continuation
            (p,mix.all,lifsup.all,mixsub,q,qsols);
          tstop(timer);
          new_line(outfile);
          put_line(outfile,"RANDOM COEFFICIENT START SYSTEM :");
          new_line(outfile);
          put_line(outfile,q);
          new_line(outfile);
          put_line(outfile,"START SOLUTIONS :");
          new_line(outfile);
          put(outfile,Length_Of(qsols),Head_Of(qsols).n,qsols);
          new_line(outfile);
          print_times(outfile,timer,"Polyhedral Continuation");
    end if;
  end Main;

begin
  Read_System(infile,infilename);
  if lp = null
   then new_line;
        get(lp);
  end if;
  Close(infile);
  Create_Output_File(outfile,outfilename);
  put(outfile,lp.all);
  Main(outfile,lp.all);
  Close(outfile);
end babldmvc;
