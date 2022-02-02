with integer_io,Solutions_io;           use integer_io,Solutions_io;
with Communications_with_User;          use Communications_with_User;
with File_Operations;                   use File_Operations;
with Complex_Polynomial_Systems_io;     use Complex_Polynomial_Systems_io;
with Drivers_for_Polynomial_Continuation;
 use Drivers_for_Polynomial_Continuation;

procedure Driver_for_Own_Start_System
             ( file : in file_type; p : in Poly_Sys;
               q : out Poly_Sys; qsols : in out Solution_List ) is

  qfile : file_type;
  qq : Poly_Sys(p'range);
  found : boolean;

begin
  new_line;
  put_line("Reading the name of the file that contains the start system.");
  Read_Name_and_Open_File(qfile);
  get(qfile,qq);
  Scan_and_Skip(qfile,"SOLUTIONS",found);
  if found
   then get(qfile,qsols);
   else declare
          sfile : file_type;
        begin
          put_line("Reading the name of the file for the solutions.");
          Read_Name_and_Open_File(sfile);
          get(sfile,qsols);
          Close(sfile);
        end;
  end if;
  Close(qfile);
  Check_Continuation_Parameter(qsols);
  q := qq;
  new_line(file);
  put_line(file,"Start system delivered by user : ");
  put(file,qq);
  new_line(file);
  put_line(file,"with start solutions : "); new_line(file);
  put(file,qsols);
  new_line(file);
end Driver_for_Own_Start_System;
