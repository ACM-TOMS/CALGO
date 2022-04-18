with Communications_with_User;         use Communications_with_User;
with integer_io,Numbers_io;            use integer_io,Numbers_io;
with Timing_Package;                   use Timing_Package;
with Complex_Polynomial_Systems_io;    use Complex_Polynomial_Systems_io;
with Natural_Vectors;
with Solutions_io;                     use Solutions_io;
with Set_Structure,Set_Structure_io;
with Degree_Sets_Tables;               use Degree_Sets_Tables;
with Random_Product_System;
with Construct_Random_Product_Start_System;
 use Construct_Random_Product_Start_System;

package body Drivers_for_Set_Structures is

  procedure Set_Structure_Info is

    i : array(1..18) of string(1..65);

  begin
    i( 1):="  A generalized Bezout  number  is  based  on  a  supporting  set";
    i( 2):="structure.   A  set  structure is a tuple of arrays of subsets of";
    i( 3):="unknowns.                                                        ";
    i( 4):="  The corresponding start system is a linear-product system:  the";
    i( 5):="i-th  equation  is  the  product  of linear equations with random";
    i( 6):="coefficient in the unknowns of the set of the  i-th  array.   The";
    i( 7):="number  of  factors  in  the product for the i-th equation of the";
    i( 8):="start system equals the number of subsets in the  i-th  array  of";
    i( 9):="the set structure.                                               ";
    i(10):="  A set structure is supporting for a polynomial system if  every";
    i(11):="monomial  in  the system also occurs in the corresponding linear-";
    i(12):="product start system.                                            ";
    i(13):="  Given a supporting set structure, the generalized Bezout number";
    i(14):="equals  the  number  of  solutions  of  the corresponding linear-";
    i(15):="product start system.   Before  the  construction  of  the  start";
    i(16):="system, a generalized Bezout number is first computed in a formal";
    i(17):="way as a generalized permanent of a degree matrix.   A  heuristic";
    i(18):="procedure is available for generating a supporting set structure.";
    for k in i'range loop
      put_line(i(k));
    end loop;
  end Set_Structure_Info;

  procedure Driver_for_Set_Structure
               ( file : in file_type; p : in Poly_Sys; 
                 b : in out natural; lpos : in out List;
                 q : out Poly_Sys; qsols : out Solution_List ) is

    procedure Write_Results ( file : in file_type; bb : in natural ) is
    begin
      new_line(file);
      put(file,"  generalized Bezout number is "); put(file,bb,1);
      new_line(file);
      put_line(file,"  based on the set structure :");
      Set_Structure_io.put(file);
    end Write_Results;

    procedure Save_Results ( qq : in Poly_Sys; qqsols : in Solution_List ) is

      qqfile : file_type;

    begin
      if not Is_Null(qqsols)
       then new_line;
            put_line("Reading file name to write start system.");
            Read_Name_and_Create_File(qqfile);
            put_line(qqfile,qq);
            new_line(qqfile);
            put_line(qqfile,"THE SOLUTIONS : ");
            new_line(qqfile);
            put(qqfile,Length_Of(qqsols),Head_Of(qqsols).n,qqsols);
            Close(qqfile);
      end if;
    end Save_Results;

    procedure Display_Menu ( choice : out character; bb : in natural ) is

      ans : character;

    begin
      new_line;
      put_line("MENU for generalized Bezout Numbers based on Set Structures :");
      put     ("  0. exit - current Bezout number is "); put(bb,1); new_line;
      put_line("  1. Apply heuristic constructor for set structure");
      put_line("  2. Evaluate your own set structure");
      put("Type 0, 1, or 2 to make your choice : ");
      Ask_Alternative(ans,"012"); choice := ans;
    end Display_Menu;

    procedure Dispatch_Menu ( file : in file_type;
                              choice : in character; bb : in out natural ) is
    begin
      case choice is
        when '1' =>  
          Construct_Random_Product_Start_System.Build_Set_Structure(p);
          bb := Permanent(Degree_Sets_Tables.Create);
        when '2' => 
          declare
            ns : Natural_Vectors.Vector(p'range);
          begin
            for i in ns'range loop
              put("  Give the number of sets for polynomial ");
              put(i,1); put(" : "); Read_Natural(ns(i));
            end loop;
            Set_Structure.Init(ns);
            put_line("Give the set structure :");
            Set_Structure_io.get;
          end;
          bb := Permanent(Degree_Sets_Tables.Create);
        when others => null;
      end case;
      Write_Results(Standard_Output,bb); Write_Results(file,bb);
    end Dispatch_Menu;

    procedure Driver_for_Bezout_Number
                  ( file : in file_type; bb : in out natural ) is

      method : character;
      timer : timing_widget;

    begin
      new_line(file);
      put_line(file,"SET STRUCTURE ANALYSIS :");
      tstart(timer);
      loop
        Display_Menu(method,bb);
        exit when method = '0';
        Dispatch_Menu(file,method,bb);
      end loop;
      tstop(timer);
      new_line(file);
      print_times(file,timer,"set structure analysis");
    end Driver_for_Bezout_Number;

    procedure Driver_for_Start_System
                  ( file : in file_type; bb : in natural ) is

      ans : character;
      timer : timing_widget;
      qq : Poly_Sys(p'range);
      qqsols : Solution_List;

    begin
      new_line;
      put("Do you want a start system based on the set structure ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then
         declare
           nl : natural;
           n : natural := p'length;
         begin
          -- new_line;
          -- put("Solving "); put(bb,1); put(" linear systems...");
           tstart(timer);
           Random_Product_System.Init(n);
           Build_Random_Product_System(n);
           Set_Structure.Clear;
           qq := Random_Product_System.Polynomial_System;
          -- Random_Product_System.Solve(qqsols,nl,lpos);
           Random_Product_System.Solve(qqsols,nl);
           Random_Product_System.Clear;
           tstop(timer);
           Save_Results(qq,qqsols);
           q := qq; qsols := qqsols;
           new_line(file);
           put_line(file,"RANDOM LINEAR-PRODUCT START SYSTEM : ");
           put_line(file,qq);
           new_line(file);
           put_line(file,"THE SOLUTIONS :");
           new_line(file);
           put(file,Length_Of(qqsols),Head_Of(qqsols).n,qqsols);
           new_line(file);
           print_times(file,timer,"constructing and solving the start system");
         end;
       else
         Set_Structure.Clear;
      -- Clear(lpos);
      end if;
    end Driver_for_Start_System;

    procedure Main_Driver is

      bb : natural := b;

    begin
      Driver_for_Bezout_Number(file,bb);
      if not Set_Structure.Empty
       then b := bb;
            Driver_for_Start_System(file,bb);
      end if;
    end Main_Driver;

  begin
    Main_Driver;
  end Driver_for_Set_Structure;

end Drivers_for_Set_Structures;
