package body Communications_with_User is

-- AUXILIARY :

  function Is_In ( s : string; ch : character ) return boolean is

  -- DESCRIPTION :
  --   Returns true if the character occurs in the string, false otherwise.

  begin
    for i in s'range loop
      if s(i) = ch
       then return true;
      end if;
    end loop;
    return false;
  end Is_In;

-- TARGET ROUTINES :

  function Read_String return string is

    temp : string(1..80);
    cnt : natural;

  begin
    put("Give a string of characters : ");
    get_line(temp,cnt);
    return temp(1..cnt);
  end Read_String;

  procedure Ask ( ans : out character ) is

    ch : character;

  begin
    loop
      get(ch); skip_line;
      exit when Valid_Alternative(ch);
      put("Invalid alternative.  Please try again : ");
    end loop;
    ans := ch;
  end Ask;

  procedure Ask_Yes_or_No ( ans : out character ) is

    function Yes_or_No ( alt : character ) return boolean is
    begin
      if alt = 'y' or else alt = 'n'
       then return true;
       else return false;
      end if;
    end Yes_or_No;
    procedure Yes_or_No_Ask is new Ask (Yes_or_No);

  begin
    Yes_or_No_Ask(ans);
  end Ask_Yes_or_No;

  procedure Ask_Alternative ( ans : out character; alternatives : in string ) is

    function Is_Valid ( alt : character ) return boolean is
    begin
      return Is_In(alternatives,alt);
    end Is_Valid;
    procedure Ask_Alt is new Ask ( Is_Valid );

  begin
    Ask_Alt(ans);
  end Ask_Alternative;

  procedure Ask_Alternative
                ( ans : in out string; alternatives : string;
                  prefix : in character ) is

    ok : boolean := false;
    tmp : string(1..10);
    ind,cnt : natural;

  begin
    loop
      get_line(tmp,cnt);
      ans := "  ";
      ind := 1;
      while (ans(1) = ' ') and (ind <= cnt) loop
        ans(1) := tmp(ind);
        ind := ind+1;
      end loop;
      if ans(1) = prefix
       then while (ans(2) = ' ') and (ind <= cnt) loop
              ans(2) := tmp(ind);
              ind := ind+1;
            end loop;
            if Is_In(alternatives,ans(2))
             then ok := true;
             else put("Invalid alternative.  Please try again : ");
            end if;
       else if Is_In(alternatives,ans(1))
             then ok := true;
             else put("Invalid alternative.  Please try again : ");
            end if;
      end if;
      exit when ok;
    end loop;
  end Ask_Alternative;

  procedure Read_Name_and_Open_File ( file : in out file_type ) is

    name : constant string := Read_String;

  begin
    Open(file,in_file,name);
  exception
    when NAME_ERROR => 
       put_line("The file could not be located, please try again...");
       Read_Name_and_Open_File(file);
    when USE_ERROR =>
       put_line("File is not readable, please try again...");
       Read_Name_and_Open_File(file);
  end Read_Name_and_Open_File;

  procedure Read_Name_and_Create_File ( file : in out file_type ) is

    filename : constant string := Read_String;
    ans : character;
    temp : file_type;

    procedure Retry is
    begin
      Create(file,out_file,filename);
    exception
      when USE_ERROR =>
        put_line("Could not create file, file already in use.");
        put_line("Please, try again...");
        Read_Name_and_Create_File(file);
      when NAME_ERROR =>
        put_line("Could not create file, perhaps wrong directory ?");
        put_line("Please, try again...");
        Read_Name_and_Create_File(file);
    end Retry;

  begin
    Open(temp,in_file,filename);
    Close(temp);
    put("There exists already a file named "); put_line(filename);
    put("Do you want to destroy this file ? (y/n) ");
    Ask_Yes_or_No(ans);
    if ans = 'y'
     then create(file,out_file,filename);
     else Read_Name_and_Create_File(file);
    end if;
  exception
    when others => Retry;
  end Read_Name_and_Create_File;

  procedure Open_Input_File
               ( file : in out file_type; filename : in string ) is
  begin
    Open(file,in_file,filename);
  exception
    when NAME_ERROR =>
       put("The file "); put(filename);
       put_line(" could not be located, please try again...");
       Read_Name_and_Open_File(file);
    when USE_ERROR =>
       put("The file "); put(filename);
       put_line(" is not readable, please try again...");
       Read_Name_and_Open_File(file);
  end Open_Input_File;

  procedure Create_Output_File
                 ( file : in out file_type; filename : in string ) is

    ans : character;
    temp : file_type;

    procedure Retry is
    begin
      Create(file,out_file,filename);
    exception
      when USE_ERROR =>
        put("Could not create file "); put(filename);
        put_line(", file already in use.");
        put_line("Please, try again...");
        Read_Name_and_Create_File(file);
      when NAME_ERROR =>
        put("Could not create file "); put(filename);
        put_line(", perhaps wrong directory ?");
        put_line("Please, try again...");
        Read_Name_and_Create_File(file);
    end Retry;

  begin
    if filename = ""
     then new_line;
          put_line("Reading the name of the output file.");
          Read_Name_and_Create_File(file);
     else Open(temp,in_file,filename); Close(temp);
          new_line;
          put("There exists already a file named "); put_line(filename);
          put("Do you want to destroy this file ? (y/n) ");
          Ask_Yes_or_No(ans);
          if ans = 'y'
           then create(file,out_file,filename);
           else Read_Name_and_Create_File(file);
          end if;
    end if;
  exception
    when others => Retry;
  end Create_Output_File;

end Communications_with_User;
