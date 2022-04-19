with text_io;                          use text_io;
with Communications_with_User;         use Communications_with_User;
with Complex_Polynomial_Systems;       use Complex_Polynomial_Systems;
with Complex_Polynomial_Systems_io;    use Complex_Polynomial_Systems_io;
with Solutions;                        use Solutions;
with Black_Box_Root_Counting;

procedure bablroco ( infilename,outfilename : in string ) is

  lp,lq : Link_to_Poly_Sys;
  infile,outfile : file_type;
  rc : natural;
  roco,poco : duration;
  qsols : Solution_List;

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

begin
  Read_System(infile,infilename);
  if lp = null
   then new_line;
        get(lp);
  end if;
  Create_Output_File(outfile,outfilename);
  put(outfile,lp.all);
  lq := new Poly_Sys(lp'range);
  Black_Box_Root_Counting(outfile,lp.all,rc,lq.all,qsols,roco,poco);
end bablroco;
