package body Strings_to_Natural_Numbers is

  function Convert ( c : character ) return natural is
  begin
    case c is
      when '0' => return 0;
      when '1' => return 1;
      when '2' => return 2;
      when '3' => return 3;
      when '4' => return 4;
      when '5' => return 5;
      when '6' => return 6;
      when '7' => return 7;
      when '8' => return 8;
      when '9' => return 9;
      when others => return 10;
    end case;
  end convert;

  function Convert ( s : string ) return integer is

    res : integer := -1;
    num : natural;

  begin
    for i in s'range loop
      num := convert(s(i));
      if num < 10
       then if res < 0
             then res := num;
             else res := 10*res + num;
            end if;
      end if;
      exit when (res > 0) and (num = 10);
    end loop;
    return res;
  end Convert;

end Strings_to_Natural_Numbers;
