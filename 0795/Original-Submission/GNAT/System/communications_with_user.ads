with text_io;  use text_io;

package Communications_with_User is

-- DESCRIPTION :
--   This package provides some routines which make the
--   communication with the user more easy.
--   They help in catching typing errors.

  function Read_String return string;

  -- DESCRIPTION :
  --   Reads a string from standard input and returns it to the caller.
  --   The length of the input string must be smaller than 80 characters.

  generic
    with function Valid_Alternative ( alt : character ) return boolean;
  procedure Ask ( ans : out character );

  -- DESCRIPTION :
  --   This procedure keeps reading a character from standard input,
  --   until a valid one has been given.
  --   The function Valid_Alternative decides whether a certain
  --   character corresponds to a valid answer to the question.

  procedure Ask_Yes_or_No ( ans : out character );

  -- DESCRIPTION :
  --   Keeps reading a character from standard output,
  --   until the user gives a 'y' or a 'n'.

  procedure Ask_Alternative ( ans : out character; alternatives : in string );

  -- DESCRIPTION :
  --   This procedure keeps reading a character from standard input,
  --   until a character that belongs to the string s has been given.

  procedure Ask_Alternative
                ( ans : in out string; alternatives : string;
                  prefix : in character );

  -- DESCRIPTION :
  --   Ask the user to give a character that occurs in the string of
  --   alternatives, eventually preceded by the given prefix character.
  --   This procedure keeps reading till a valid choice has been made.

  -- REQUIRED : ans'range = 1..2.

  procedure Read_Name_and_Open_File ( file : in out file_type );

  -- DESCRIPTION :
  --   This procedure reads a name from standard input and
  --   tries to open this file for input.
  --   If this is unsuccesful, then another name will be asked.

  procedure Read_Name_and_Create_File ( file : in out file_type );

  -- DESCRIPTION :
  --   This procedure reads a name from standard input and
  --   a file with this name will be created for output.
  --   If a file with the given name already exists,
  --   the user will be asked if the existing file may be destroyed.

  procedure Open_Input_File
               ( file : in out file_type; filename : in string );

  -- DESCRIPTION :
  --   Tries to open a file for input, starting with the given file name.
  --   If the opening of the file with the given name is not succesful,
  --   then the procedure `Read_Name_and_Open_File' will be invoked.

  procedure Create_Output_File
               ( file : in out file_type; filename : in string );

  -- DESCRIPTION :
  --   This procedure creates an output file, starting with the given
  --   file name.  If the creation of a file with this name is unsuccesful,
  --   then the procedure `Read_Name_and_Create_File' will be invoked.

end Communications_with_User;
