with integer_io;                       use integer_io;
with Integer_Vectors_io;               use Integer_Vectors_io;
with Integer_Vectors_of_Vectors;       use Integer_Vectors_of_Vectors;
with Integer_Vectors_of_Vectors_io;    use Integer_Vectors_of_Vectors_io;

package body Simplices_io is

  procedure get ( s : in out Simplex ) is

    n : natural;

  begin
    get(n);
    declare
      v : Vector(1..n);
    begin
      get(n,v);
      s := Create(v);
    end;
  end get;

  procedure get ( n : in natural; s : in out Simplex ) is
  
    v : Vector(1..n);

  begin
    get(n,v);
    s := Create(v);
  end get;

  procedure get ( file : in file_type; s : in out Simplex ) is

    n : natural;

  begin
    get(file,n);
    declare
      v : Vector(1..n);
    begin
      get(file,n,v);
      s := Create(v);
    end;
  end get;

  procedure get ( file : in file_type; n : in natural; s : in out Simplex ) is
  
    v : Vector(1..n);

  begin
    get(file,n,v);
    s := Create(v);
  end get;

  procedure put ( s : in Simplex ) is
  begin
    put(Normal(s)); new_line;
    put(Normal(s)'last,1); new_line;
    put(Vertices(s));
  end put;

  procedure put ( file : in file_type; s : in Simplex ) is
  begin
    put(file,Normal(s)); new_line(file);
    put(file,Normal(s)'last,1); new_line(file);
    put(file,Vertices(s));
  end put;

end Simplices_io;
