package body File_Operations is

  procedure Scan ( file : in file_type; ch : in character;
                   found : out boolean ) is

    c : character;

  begin
    while not End_of_File(file) loop
      get(file,c);
      if c = ch
       then found := true;
            return;
      end if;
    end loop;
    found := false;
  end Scan;

  procedure Scan ( file : in file_type; banner : in string;
                   found : out boolean ) is

    index : natural := banner'first-1;
    ch : character;

  begin
    while not End_of_File(file) loop
      get(file,ch);
      if index < banner'first
       then
         if ch = banner(banner'first)
          then index := banner'first+1;
         end if;
       else
         if ch = banner(index)
          then index := index + 1;
          else index := banner'first-1;
         end if;
      end if;
      exit when index > banner'last;
    end loop;
    if index > banner'last
     then found := true;
     else found := false;
    end if;
  exception
    when others => found := false; return;
  end Scan;

  procedure Scan_and_Skip ( file : in file_type; banner : in string;
                            found : out boolean ) is

    fnd : boolean;

  begin
    Scan(file,banner,fnd);
    if fnd
     then skip_line(file);
    end if;
    found := fnd;
  end Scan_and_Skip;

end File_Operations;
