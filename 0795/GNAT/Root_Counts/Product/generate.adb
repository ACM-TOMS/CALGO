procedure generate ( k,first,last : in integer ) is

  acc : boolean_array(first..last) := (first..last => false);
  cont : boolean := true;

  procedure generate (k,start,last : in integer;
                      acc : in out boolean_array) is
  begin
    if k = 0
     then process(acc,cont);
     elsif k > last - start + 1
         then return;
         else for i in start..last loop
                acc(i) := true;
                generate(k-1,i+1,last,acc);
                exit when not cont;
                acc(i) := false;
              end loop;
    end if;
  end generate;

begin
  generate(k,first,last,acc);
end generate;
