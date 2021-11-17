with Communications_with_User;          use Communications_with_User;
with integer_io,Numbers_io;             use integer_io,Numbers_io;
with Timing_Package;                    use Timing_Package;
with Complex_Polynomial_Systems_io;     use Complex_Polynomial_Systems_io;
with Solutions_io;                      use Solutions_io;
with Natural_Vectors,Degree_Structure;  use Natural_Vectors,Degree_Structure;
with Sets_of_Unknowns;                  use Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns;    use Partitions_of_Sets_of_Unknowns;
with Partitions_of_Sets_of_Unknowns_io; use Partitions_of_Sets_of_Unknowns_io;
with Random_Product_Start_Systems;      use Random_Product_Start_Systems;

package body Drivers_for_Multi_Homogenization is

  procedure Multi_Homogenization_Info is

    i : array(1..17) of string(1..65);

  begin
    i( 1):="  A multi-homogeneous Bezout  number  is  based  on  a  tuple  of";
    i( 2):="partitions  of  the set of unknowns.  For every polynomial in the";
    i( 3):="system, a different partition can model its structure.           ";
    i( 4):="  The corresponding start system is a linear-product system:  the";
    i( 5):="i-th  equation  is  the  product  of linear equations with random";
    i( 6):="coefficients in the unknowns of the set of  the  partition.   The";
    i( 7):="number  of  factors  in  the product for the i-th equation of the";
    i( 8):="start system equals the  product  of  the  degrees  of  the  i-th";
    i( 9):="polynomial  in  the  original  system  w.r.t.  every  set  in the";
    i(10):="partition.                                                       ";
    i(11):="  Given a  tuple  of  partitions,  the  multi-homogeneous  Bezout";
    i(12):="number  equals  the  number  of  solutions  of  the corresponding";
    i(13):="linear-product start system.   Before  the  construction  of  the";
    i(14):="start system, a multi-homogeneous Bezout number is first computed";
    i(15):="in a formal way as a generalized permanent of a degree matrix.  A";
    i(16):="heuristic  procedure  is  available  for  generating  a  tuple of";
    i(17):="partitions.                                                      ";
    for k in i'range loop
      put_line(i(k));
    end loop;
  end Multi_Homogenization_Info;

  procedure Driver_for_Multi_Homogenization
                ( file : in file_type; p : in Poly_Sys; b : in out natural;
                  q : out Poly_Sys; qsols : out Solution_List ) is

    procedure Write_Results ( file : in file_type; gb : in natural ) is

      m : natural;

    begin
      new_line(file);
      put(file,"  multi-homogeneous Bezout number is ");
      put(file,gb,1); new_line(file);
      put_line(file,"  with partitions :");
      for i in p'range loop
        m := Degree_Structure.Get(i);
        declare
          z : partition(1..m);
          dg : Natural_Vectors.Vector(1..m);
        begin
          Degree_Structure.Get(i,z,dg);
          put(file,"     partition for equation "); put(file,i,2);
          put(file," : "); put(file,z); new_line(file);
          Clear(z);
        end;
      end loop;
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

    procedure Display_Menu ( choice : out character; gb : in natural ) is

      ans : character;

    begin
      new_line;
      put_line("MENU for Multi-Homogeneous Bezout Numbers :");
      put     ("  0. exit - current Bezout number is "); put(gb,1); new_line;
      put_line("  1. Apply heuristic partitioner");
      put_line("  2. Evaluate your own tuple of partitions.");
      put("Type 0, 1, or 2 to make your choice : "); 
      Ask_Alternative(ans,"012"); choice := ans;
    end Display_Menu;

    procedure Dispatch_Menu ( choice : in character; gb : in out natural ) is

      m : natural;

    begin
      case choice is
        when '1' => gb := Degree_Structure.Generalized_Bezout_Number(p);
        when '2' => for i in p'range loop
                      put("Give the number of sets in partition ");
                      put(i,1); put(" : "); Read_Natural(m);
                      put("Give "); put(m,1); put(" sets : ");
                      declare
                        zz : Partition(1..m);
                      begin
                        Create(zz,p'length); get(zz);
                        Degree_Structure.Put(p,i,m,zz);
                        Clear(zz);
                      end;
                    end loop;
                    gb := Degree_Structure.Generalized_Bezout_Number;
        when others => null;
      end case;
      Write_Results(Standard_Output,gb); Write_Results(file,gb);
    end Dispatch_Menu;

    procedure Driver_for_Partitions
                  ( file : in file_type; gb : in out natural ) is

      timer : timing_widget;
      method : character;

    begin
      new_line(file);
      put_line(file,"MULTI-HOMOGENIZATION :");
      tstart(timer);
      loop
        Display_Menu(method,gb);
        exit when method = '0';
        Dispatch_Menu(method,gb);
      end loop;
      tstop(timer);
      b := gb;
      new_line(file);
      print_times(file,timer,"Computation of multi-homogeneous Bezout number");
    end Driver_for_Partitions;

    procedure Driver_for_Start_System ( file : in file_type ) is

      ans : character;
      timer : timing_widget;
      qq : Poly_Sys(p'range);
      qqsols : Solution_List;

    begin
      new_line;
      put("Do you want a multi-homogeneous start system ? (y/n) ");
      Ask_Yes_or_No(ans);
      if ans = 'y'
       then 
         tstart(timer);
         GBQ(p,qq,qqsols);
         tstop(timer);
         Save_Results(qq,qqsols);
         new_line(file);
         put_line(file,"MULTI-HOMOGENEOUS START SYSTEM :");
         put_line(file,qq);
         q := qq; qsols := qqsols;
         new_line(file);
         put_line(file,"THE SOLUTIONS :");
         new_line(file);
         put(file,Length_Of(qqsols),Head_Of(qqsols).n,qqsols);
         new_line(file);
         print_times(file,timer,
                     "Construction of multi-homogeneous start system");
      end if;
    end Driver_for_Start_System;

    procedure Main_Driver is

      gb : natural := b;

    begin
      Driver_for_Partitions(file,gb);
      if not Degree_Structure.Empty
       then Driver_for_Start_System(file);
      end if;
    end Main_Driver;

  begin
    Main_Driver;
  end Driver_for_Multi_Homogenization;

end Drivers_for_Multi_Homogenization;
