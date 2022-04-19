with text_io;                            use text_io;
with Communications_with_User;           use Communications_with_User;
with Black_Box_Polynomial_Continuations; use Black_Box_Polynomial_Continuations;

procedure bablpoco ( infilename,outfilename : in string ) is

  infile,outfile : file_type;
  poco : duration;

begin
  if infilename /= ""
   then Open_Input_File(infile,infilename);
   else new_line;
        put_line("Reading the name of the input file.");
        Read_Name_and_Open_File(infile);
  end if;
  Create_Output_File(outfile,outfilename);
  Black_Box_Polynomial_Continuation(infile,outfile,poco);
  Close(infile); Close(outfile);
end bablpoco;
