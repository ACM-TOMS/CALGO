with Mathematical_Functions;       use Mathematical_Functions;
with Complex_Numbers;              use Complex_Numbers;

package body Complex_Norms is

  function Norm1 ( x : Vector ) return double_float is

    max,temp : double_float;

  begin
    max := modulus(x(x'first));
    for i in (x'first+1)..x'last loop
      temp := modulus(x(i));
      if temp > max
       then max := temp;
      end if;
    end loop;
    return max;
  end Norm1;

  function Norm2 ( x : Vector ) return double_float is

    sum : double_float := 0.0;

  begin
    for i in x'range loop
      sum := sum + modulus(x(i))**2;
    end loop;
    return SQRT(sum);
  end Norm2;

end Complex_Norms;
