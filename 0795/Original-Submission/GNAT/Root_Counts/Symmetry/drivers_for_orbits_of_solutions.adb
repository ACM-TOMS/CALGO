with integer_io;                          use integer_io;
with Symmetry_Group_io,Solutions_io;      use Solutions_io;
with Orbits_of_Solutions;                 use Orbits_of_Solutions;

package body Drivers_for_Orbits_of_Solutions is

  procedure Driver_for_Orbits_of_Solutions
                  ( file : in file_type; sols : in out Solution_List;
                    v : in List_of_Permutations; allperms,signsym : in boolean;
                    tol : in double_float ) is

  begin
    if not Is_Null(sols)
     then
       declare
         n : constant natural := Head_Of(sols).n;
         orb : Permutation(1..n);
       begin
         Driver_for_Orbits_of_Solutions(file,sols,v,allperms,signsym,tol,orb);
       end;
    end if;
  end Driver_for_Orbits_of_Solutions;

  procedure Driver_for_Orbits_of_Solutions
                  ( file : in file_type; sols : in out Solution_List;
                    v : in List_of_Permutations; allperms,signsym : in boolean;
                    tol : in double_float; orbi : out Permutation ) is

  begin
    if not Is_Null(sols)
     then declare
            n : constant natural := Head_Of(sols).n;
            orb : Permutation(1..n);
            len : natural;
          begin
            if allperms
             then sols := Generating(sols,signsym,tol);
             else Analyze(v,signsym,tol,sols);
            end if;
            len := Length_Of(sols);
            new_line(file);
            put(file,"The number of generating solutions : ");
            put(file,len,1); new_line(file);
            orb := Orbits(sols,tol);
            put(file,"The orbits : "); Symmetry_Group_io.put(file,orb);
            new_line(file);
            new_line(file);
            put_line(file,"THE GENERATING SOLUTIONS : ");
            put(file,len,Head_Of(sols).n,sols); new_line(file);
            orbi := orb;
          end;
    end if;
  end Driver_for_Orbits_of_Solutions;

end Drivers_for_Orbits_of_Solutions;
