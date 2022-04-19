with text_io,integer_io,Numbers_io;   use text_io,integer_io,Numbers_io;
with Communications_with_User;        use Communications_with_User;
with Complex_Polynomial_Systems;      use Complex_Polynomial_Systems;
with Complex_Polynomial_Systems_io;   use Complex_Polynomial_Systems_io;
with Solutions,Solutions_io;          use Solutions,Solutions_io;
with Driver_for_Root_Counts;
with Bye_Bye_Message;

procedure mainroco ( infilename,outfilename : in string ) is
 
  n : natural;
  inft,outft : file_type;
  lp : Link_to_Poly_Sys;

  procedure Read_System ( filename : in string ) is

    file : file_type;

  begin
    if filename /= ""
     then Open(file,in_file,filename);
          get(file,n);
          lp := new Poly_Sys(1..n);
          get(file,n,lp.all);
          Close(file);
    end if;
  exception
    when others =>
      new_line;
      put("Could not open file with name "); put_line(filename);
      lp := null; return;
  end Read_System;

begin
  Read_System(infilename);
  if lp = null
   then new_line; get(lp); new_line;
  end if;
  Create_Output_File(outft,outfilename);
  put(outft,lp.all);
  declare
    q : Poly_Sys(lp'range);
    qsols : Solution_List;
    rc : natural;
  begin
    Driver_for_Root_Counts(outft,lp.all,q,false,qsols,rc);
  end; 
  new_line(outft);
  put(outft,Bye_Bye_Message);
  Close(outft);
end mainroco;
