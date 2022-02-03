with Floating_Point_Numbers;     use Floating_Point_Numbers;

package body Complex_Numbers_io is

  use Floating_Point_Numbers.double_float_io;

  procedure get ( c : out double_complex ) is
    x,y : double_float;
  begin
    get(x); get(y);
    c := CMPLX(x,y);
  end get;

  procedure get ( f : in file_type; c : out double_complex ) is
    x,y : double_float;
  begin
    get(f,x); get(f,y);
    c := CMPLX(x,y);
  end get;

  procedure put ( c : in double_complex ) is
  begin
    put(REAL_PART(c));
    put("  ");
    put(IMAG_PART(c));
  end put;

  procedure put ( f : in file_type; c : in double_complex ) is
  begin
    put(f,REAL_PART(c));
    put(f,"  ");
    put(f,IMAG_PART(c));
  end put;

  procedure put ( c : in double_complex; fore,aft,exp : in natural ) is
  begin
    put(REAL_PART(c),fore,aft,exp);
    put("  ");
    put(IMAG_PART(c),fore,aft,exp);
  end put;

  procedure put ( f : in file_type;
                  c : in double_complex; fore,aft,exp : in natural ) is
  begin
    put(f,REAL_PART(c),fore,aft,exp);
    put(f,"  ");
    put(f,IMAG_PART(c),fore,aft,exp);
  end put;

end Complex_Numbers_io;
