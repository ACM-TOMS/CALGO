with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Timing_Package;                     use Timing_Package;
with Complex_Polynomial_Systems;         use Complex_Polynomial_Systems;
with Complex_Polynomial_Systems_io;      use Complex_Polynomial_Systems_io;
with Solutions;                          use Solutions;
with Scaling;                            use Scaling;
with Black_Box_Root_Counting;
with Black_Box_Polynomial_Continuations; use Black_Box_Polynomial_Continuations;

procedure bablphc ( infilename,outfilename : in string ) is

  procedure Read_System ( file : in out file_type; filename : in string;
                          p : in out Link_to_Poly_Sys ) is
  begin
    if filename /= ""
     then Open_Input_File(file,filename);
          get(file,p);
    end if;
  exception
    when others => put_line("Something is wrong with argument file...");
                   p := null; return;
  end Read_System;

  procedure Timing_Summary ( file : in file_type;
                             roco,hoco,poco,total : in duration ) is

    b0 : constant string :=
     "  ---------------------------------------------------------------------";
    b1 : constant string :=
     "  |                    TIMING INFORMATION SUMMARY                     |";
    b2 : constant string :=
     "  |   root counts  |  start system  |  continuation  |   total time   |";

  begin
    put_line(file,b0);
    put_line(file,b1);
    put_line(file,b0);
    put_line(file,b2);
    put_line(file,b0);
    put(file,"  | ");
    print_hms(file,roco); put(file," | ");
    print_hms(file,hoco); put(file," | ");
    print_hms(file,poco); put(file," | ");
    print_hms(file,total); put_line(file," |");
    put_line(file,b0);
  end Timing_Summary;

  procedure Main is

    timer : timing_widget;
    infile,outfile : file_type;
    p,q : Link_to_Poly_Sys;
    rc : natural;
    sols : Solution_List;
    roco,hoco,poco,total : duration;

  begin
    Read_System(infile,infilename,p);
    if p = null
     then new_line; get(p);
    end if;
    Create_Output_File(outfile,outfilename);
    put(outfile,p.all);
    q := new Poly_Sys(p'range);
    tstart(timer);
    Black_Box_Root_Counting(outfile,p.all,rc,q.all,sols,roco,hoco);
    if rc /= 0
     then
       Scale(p.all);
       Black_Box_Polynomial_Continuation(infile,outfile,p.all,q.all,sols,poco);
    end if;
    tstop(timer);
    total := Elapsed_User_Time(timer);
    Close(infile);
    new_line(outfile);
    print_times(outfile,timer,"Solving the polynomial system");
    new_line(outfile);
    Timing_Summary(outfile,roco,hoco,poco,total);
    Close(outfile);
  end Main;

begin
  Main;
end bablphc;
