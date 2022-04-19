with Communications_with_User;         use Communications_with_User;
with text_io,integer_io,Solutions;     use text_io,integer_io,Solutions;
with Floating_Point_Numbers;           use Floating_Point_Numbers;
with Complex_Polynomial_Systems;       use Complex_Polynomial_Systems;
with Complex_Polynomial_Systems_io;    use Complex_Polynomial_Systems_io;
with Root_Refiners;                    use Root_Refiners;

with Drivers_for_Implicit_Lifting;     use Drivers_for_Implicit_Lifting;
with Drivers_for_Static_Lifting;       use Drivers_for_Static_Lifting;
with Drivers_for_Dynamic_Lifting;      use Drivers_for_Dynamic_Lifting;
with Drivers_for_Symmetric_Lifting;    use Drivers_for_Symmetric_Lifting;

procedure mainsmvc ( infilename,outfilename : in string ) is

  outft : file_type;
  lp : Link_to_Poly_Sys;
  ans : character;

  procedure Read_System ( filename : in string ) is
  
    file : file_type;
    n : natural;

  begin
    if filename /= ""
     then Open(file,in_file,filename);
          get(file,n);
          lp := new Poly_Sys(1..n);
          get(file,n,lp.all);
          Close(file);
    end if;
  exception
    when others => put_line("Something is wrong with argument file...");
                   lp := null; return;
  end Read_System;

  function Lifting_Strategy return natural is

    choice : string(1..2) := "  ";

  begin
    loop
      new_line;
      put_line("MENU with available Lifting Strategies :");
      put_line("  1. Implicit lifting  : based on recursive formula.");
      put_line("  2. Static lifting    : lift points and prune lower hull.");
      put_line("  3. Dynamic lifting   : incrementally add the points.");
      put_line
          ("  4. Symmetric lifting : points in same orbit get same lifting.");
      put("Type 1, 2, 3, or 4 to select lifting,"
                & " eventually preceded by i for info : ");
      Ask_Alternative(choice,"1234",'i');
      exit when choice(1) /= 'i';
      new_line;
      case choice(2) is
        when '1' => Implicit_Lifting_Info; new_line;
                    put("Do you want to apply implicit lifting ? (y/n) ");
                    Ask_Yes_or_No(ans);
                    if ans = 'y'
                     then choice(1) := '1';
                    end if;
        when '2' => Static_Lifting_Info; new_line;
                    put("Do you want to apply static lifting ? (y/n) ");
                    Ask_Yes_or_No(ans);
                    if ans = 'y'
                     then choice(1) := '2';
                    end if;
        when '3' => Dynamic_Lifting_Info; new_line;
                    put("Do you want to apply dynamic lifting ? (y/n) ");
                    Ask_Yes_or_No(ans);
                    if ans = 'y'
                     then choice(1) := '3';
                    end if;
        when '4' => Symmetric_Lifting_Info; new_line;
                    put("Do you want to apply implicit lifting ? (y/n) ");
                    Ask_Yes_or_No(ans);
                    if ans = 'y'
                     then choice(1) := '4';
                    end if;
        when others => put_line("No information available.");
      end case;
      exit when choice(1) /= 'i';
    end loop;
    case choice(1) is
      when '1'    => return 1;
      when '2'    => return 2;
      when '3'    => return 3;
      when others => return 4;
    end case;
  end Lifting_Strategy;

begin
  Read_System(infilename);
  if lp = null
   then new_line; get(lp);
  end if;
  declare
    q : Poly_Sys(lp'range);
    qsols : Solution_List;
    mv : natural;
    strategy : natural;
  begin
    Create_Output_File(outft,outfilename);
    put(outft,lp.all);
    strategy := Lifting_Strategy;
    new_line(outft);
    case strategy is
      when 1 => put_line(outft,"IMPLICIT LIFTING");
                Driver_for_Mixture_Bezout_BKK(outft,lp.all,true,q,qsols,mv);
      when 2 => put_line(outft,"STATIC LIFTING");
                Driver_for_Mixed_Volume_Computation
                                             (outft,lp.all,true,q,qsols,mv);
      when 3 => put_line(outft,"DYNAMIC LIFTING");
                Driver_for_Dynamic_Mixed_Volume_Computation
                                             (outft,lp.all,true,q,qsols,mv);
      when others => put_line(outft,"SYMMETRIC LIFTING"); 
                Driver_for_Symmetric_Mixed_Volume_Computation
                                             (outft,lp.all,true,q,qsols,mv);
    end case;
    if Length_Of(qsols) > 0
     then declare
            epsxa,epsfa : constant double_float := 10.0**(-8);
            tolsing : constant double_float := 10.0**(-8);
            nb : natural := 0;
          begin
            new_line(outft);
            Reporting_Root_Refiner
              (outft,q,qsols,epsxa,epsfa,tolsing,nb,5,false);
          end;
    end if;
    Close(outft);
  end;
end mainsmvc;
