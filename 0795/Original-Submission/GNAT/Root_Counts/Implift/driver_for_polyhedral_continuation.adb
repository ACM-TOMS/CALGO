with Communications_with_User;        use Communications_with_User;
with Complex_Polynomial_Systems_io;   use Complex_Polynomial_Systems_io;
with Polynomial_Randomizers;          use Polynomial_Randomizers;
with Drivers_for_Polynomial_Continuation;
 use Drivers_for_Polynomial_Continuation;

procedure Driver_for_Polyhedral_Continuation
                ( file : in file_type; p : in Poly_Sys; k : in natural;
                  byebye : in boolean;
                  q : out Poly_Sys; qfile,solsfile : in out file_type;
                  tosolve,ranstart,contrep : out boolean ) is

  ans : character;
  oc : natural;
  qq : Poly_Sys(p'range);

begin
  new_line;
  put_line("MENU for Polyhedral Continuation : ");
  put_line("  0. No polyhedral continuation, leave the menu.");
  put_line("  1. Solve given system by polyhedral continuation.");
  put_line("  2. Create and solve random coefficient system.");
  put("Type 0,1, or 2 to choose : "); Ask_Alternative(ans,"012");
  tosolve := (ans /= '0');
  ranstart := (ans = '2');
  if ans /= '0'
   then 
     if ans = '2'
      then
        put_line("Reading a file name to write random coefficient system.");
        Read_Name_and_Create_File(qfile);
        qq := Complex_Randomize1(p); q := qq;
        if k = 0
         then new_line(file);
              put_line(file,"RANDOM COEFFICIENT START SYSTEM :");
              new_line(file);
              put_line(file,qq);
              put_line(qfile,qq);
        end if;
      else
        q := p;
        put_line("Reading a name of a file to write start solutions on.");
        Read_Name_and_Create_File(solsfile);
     end if;
     new_line;
     Driver_for_Continuation_Parameters(file);
     new_line;
     Driver_for_Process_io(file,oc);
     contrep := (oc /= 0);
  end if;
  if byebye
   then new_line;
        put_line("No more input expected.  See output file for results.");
        new_line;
  end if;
end Driver_for_Polyhedral_Continuation;
