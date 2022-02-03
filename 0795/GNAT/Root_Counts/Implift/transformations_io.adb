with integer_io;               use integer_io;
with Integer_Vectors;
with Integer_Vectors_of_Vectors;

package body Transformations_io is

  procedure get ( n : in natural; t : out Transfo ) is
  begin
    get(Standard_Input,n,t);
  end get;

  procedure get ( file : in file_type; n : in natural; t : out Transfo ) is

    v : Integer_Vectors_of_Vectors.Vector(1..n);

  begin
    for i in v'range loop
      v(i) := new Integer_Vectors.Vector(1..n);
      for j in v(i)'range loop
	get(file,v(i)(j));
      end loop;
    end loop;
    t := Create(v);
    Integer_Vectors_of_Vectors.Clear(v);
  end get;

  procedure put ( t : in Transfo ) is
  begin
    put(Standard_Output,t);
  end put;

  procedure put ( file : in file_type; t : in Transfo ) is

    n : natural := Dimension(t);
    v : Integer_Vectors_of_Vectors.Vector(1..n);

  begin
    for i in v'range loop
      v(i) := new Integer_Vectors.Vector'(1..n => 0);
      v(i)(i) := 1;
      Apply(t,v(i).all);
    end loop;
    for i in v'range loop
      for j in v(i)'range loop
	put(file,' '); put(file,v(j)(i),1);
      end loop;
      new_line(file);
    end loop;
    Integer_Vectors_of_Vectors.Clear(v);
  end put;

end Transformations_io;
