with integer_io;   use integer_io;

package body Orbits_of_Solutions_io is

  procedure put ( lorb : in List_of_Orbits ) is

    tmp : List_of_Orbits;
    totgen,totsols : natural;

  begin
    tmp := lorb;
    if not Is_Null(lorb)
     then put_line("------------------------------------------------");
          put_line("|    ORBIT INFORMATION OF LIST OF SOLUTIONS    |");
	  put_line("------------------------------------------------");
	  put_line("|       TYPE        | NB <> | NB GEN | NB SOLS |");
	  put_line("------------------------------------------------");
	  totgen := 0;
	  totsols := 0;
          while not Is_Null(tmp) loop
            declare
	      n : natural := Head_Of(tmp).n;
	      orbi : Orbit(n) := Head_Of(tmp).all;
            begin
	      put("|  ");
	      for i in orbi.orb'range loop
	        put(orbi.orb(i),1);
	        put(' ');
	      end loop;
	      for i in 1..(17-2*n) loop
		put(' ');
              end loop;
	      put("|  ");
	      put(orbi.nbdiff,2);
	      put("   |   ");
	      put(orbi.nbgen,2);  totgen := totgen + orbi.nbgen;
	      put("   |   ");
	      put(orbi.nbsols,3); totsols := totsols + orbi.nbsols;
	      put_line("   |");
            end;
            tmp := Tail_Of(tmp);
          end loop;
	  put_line("------------------------------------------------");
	  put("| Total number of generating solutions : ");
	  put(totgen,4); put_line("  |");
	  put("| Total number of solutions generated  : ");
	  put(totsols,4); put_line("  |");
	  put_line("------------------------------------------------");
    end if;
  end put;

  procedure put ( file : in file_type; lorb : in List_of_Orbits ) is

    tmp : List_of_Orbits;
    totgen,totsols : natural;

  begin
    tmp := lorb;
    if not Is_Null(lorb)
     then put_line(file,"------------------------------------------------");
          put_line(file,"|    ORBIT INFORMATION OF LIST OF SOLUTIONS    |");
	  put_line(file,"------------------------------------------------");
	  put_line(file,"|       TYPE        | NB <> | NB GEN | NB SOLS |");
	  put_line(file,"------------------------------------------------");
	  totgen := 0;
	  totsols := 0;
          while not Is_Null(tmp) loop
            declare
	      n : natural := Head_Of(tmp).n;
	      orbi : Orbit(n) := Head_Of(tmp).all;
            begin
	      put(file,"|  ");
	      for i in orbi.orb'range loop
	        put(file,orbi.orb(i),1);
	        put(file,' ');
	      end loop;
	      for i in 1..(17-2*n) loop
		put(file,' ');
              end loop;
	      put(file,"|  ");
	      put(file,orbi.nbdiff,2);
	      put(file,"   |   ");
	      put(file,orbi.nbgen,2);  totgen := totgen + orbi.nbgen;
	      put(file,"   |   ");
	      put(file,orbi.nbsols,3); totsols := totsols + orbi.nbsols;
	      put_line(file,"   |");
            end;
            tmp := Tail_Of(tmp);
          end loop;
	  put_line(file,"------------------------------------------------");
	  put(file,"| Total number of generating solutions : ");
	  put(file,totgen,4); put_line(file,"  |");
	  put(file,"| Total number of solutions generated  : ");
	  put(file,totsols,4); put_line(file,"  |");
	  put_line(file,"------------------------------------------------");
    end if;
  end put;

end Orbits_of_Solutions_io;
