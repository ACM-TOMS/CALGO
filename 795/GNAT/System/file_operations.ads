with text_io;    use text_io;

package File_Operations is

-- DESCRIPTION :
--   This package collects some useful operations on text files.

  procedure Scan ( file : in file_type; ch : in character;
                   found : out boolean );

  -- DESCRIPTION :
  --   Scans the file on the search of a character.
  --   When the character has been found, then the procedure stops
  --   and sets found to true, otherwise, found will be false and
  --   End_of_File(file) will be true.

  procedure Scan ( file : in file_type; banner : in string;
                   found : out boolean );

  -- DESCRIPTION :
  --   Scans the file on the search of a text banner.

  -- ON INPUT :
  --   file        a file opened for input;
  --   banner      a string.

  -- ON RETURN :
  --   found       true if the banner has been found, false otherwise.

  procedure Scan_and_Skip ( file : in file_type; banner : in string;
                            found : out boolean );

  -- DESCRIPTION :
  --   The line where the banner has been found will be skipped.

end File_Operations;
