package body Complex_Matrices is

-- MATRIX-MATRIX OPERATIONS :

  function "+" ( a,b : Matrix ) return Matrix is
    r : Matrix(a'range(1),a'range(2));
  begin
    for i in r'range(1) loop
      for j in r'range(2) loop
        r(i,j) := a(i,j) + b(i,j);
      end loop;
    end loop;
    return r;
  end "+";

  function "-" ( a,b : Matrix ) return Matrix is
    r : Matrix(a'range(1),a'range(2));
  begin
    for i in r'range(1) loop
      for j in r'range(2) loop
        r(i,j) := a(i,j) - b(i,j);
      end loop;
    end loop;
    return r;
  end "-";

  function "-" ( a   : Matrix ) return Matrix is
    r : Matrix(a'range(1),a'range(2));
  begin
    for i in r'range(1) loop
      for j in r'range(2) loop
        r(i,j) := -a(i,j);
      end loop;
    end loop;
    return r;
  end "-";

  function "*" ( a,b : Matrix ) return Matrix is
    r : Matrix(a'range(1),b'range(2));
  begin
    for i in r'range(1) loop
      for j in r'range(2) loop
        r(i,j) := CMPLX(0.0);
        for k in a'range(2) loop
          r(i,j) := r(i,j) + a(i,k)*b(k,j);
        end loop;
      end loop;
    end loop;
    return r;
  end "*";

  procedure Mult1 ( a : in out Matrix; b : in Matrix ) is
    temp : Vector(a'range(2));
  begin
    for i in a'range(1) loop
      for j in b'range(2) loop
        temp(j) := CMPLX(0.0);
        for k in a'range(2) loop
          temp(j) := temp(j) + a(i,k)*b(k,j);
        end loop;
      end loop;
      for j in a'range(2) loop
        a(i,j) := temp(j);
      end loop;
    end loop;
  end Mult1;

  procedure Mult2 ( a : in Matrix; b : in out Matrix ) is
    temp : Vector(a'range(1));
  begin
    for i in b'range(2) loop
      for j in a'range(1) loop
        temp(j) := CMPLX(0.0);
        for k in a'range(1) loop
          temp(j) := temp(j) + a(j,k)*b(k,i);
        end loop;
      end loop;
      for j in b'range(1) loop
        b(j,i) := temp(j);
      end loop;
    end loop;
  end Mult2;

-- MATRIX-VECTOR OPERATIONS :

  function "*" ( a : Matrix; v : Vector ) return Vector is
    r : Vector(a'range(1));
  begin
    for i in r'range loop
      r(i) := CMPLX(0.0);
      for j in a'range(2) loop
        r(i) := r(i) + a(i,j)*v(j);
      end loop;
    end loop;
    return r;
  end "*";

  procedure Mult ( a : in Matrix; v : in out Vector ) is
    iv : Vector(v'range);
  begin
    for i in iv'range loop
      iv(i) := CMPLX(0.0);
      for j in a'range(2) loop
        iv(i) := iv(i) + a(i,j)*v(j);
      end loop;
    end loop;
    for i in v'range loop
      v(i) := iv(i);
    end loop;
  end Mult;

end Complex_Matrices;
