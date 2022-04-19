with text_io,integer_io;       use text_io,integer_io;

package body Numbers_io is

  use Floating_Point_Numbers.single_float_io;
  use Floating_Point_Numbers.double_float_io;

  procedure Read_Positive ( p : out positive ) is
  begin
    get(p); skip_line;
  exception
    when DATA_ERROR | CONSTRAINT_ERROR =>
      skip_line;  -- skip the rubbish
      put("This is not a positive number, try again : ");
      Read_Positive(p);
  end Read_Positive;

  procedure Read_Natural ( n : out natural ) is
  begin
    get(n); skip_line;
  exception
    when DATA_ERROR | CONSTRAINT_ERROR =>
      skip_line;  -- skip the rubbish
      put("This is not a natural number, try again : ");
      Read_Natural(n);
  end Read_Natural;

  procedure Read_Integer ( i : out integer ) is
  begin
    get(i); skip_line;
  exception
    when DATA_ERROR | CONSTRAINT_ERROR =>
      skip_line;  -- skip the rubbish
      put("This is not an integer number, try again : ");
      Read_Integer(i);
  end Read_Integer;

  procedure Read_Single_Float ( f : out single_float ) is
  begin
    get(f); skip_line;
  exception
    when DATA_ERROR | CONSTRAINT_ERROR =>
      skip_line;  -- skip the rubbish
      put("This is not a floating point number, try again : ");
      Read_Single_Float(f);
  end Read_Single_Float;

  procedure Read_Double_Float ( f : out double_float ) is
  begin
    get(f); skip_line;
  exception
    when DATA_ERROR | CONSTRAINT_ERROR =>
      skip_line;  -- skip the rubbish
      put("This is not a floating point number, try again : ");
      Read_Double_Float(f);
  end Read_Double_Float;

end Numbers_io;
