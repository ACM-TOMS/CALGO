package body Exponent_Vectors is

-- CREATORS :

  function Create ( p : Complex_Multivariate_Laurent_Polynomials.Poly )
                  return Integer_Vectors_of_Vectors.Vector is

    use Complex_Multivariate_Laurent_Polynomials;
    res : Integer_Vectors_of_Vectors.Vector(1..Number_of_Terms(p));
    ind : natural := 0;
    
    procedure Add_Exponent ( t : in Term; continue : out boolean ) is
    begin
      ind := ind + 1;
      res(ind) := new Integer_Vectors.Vector(t.dg'range);
      for i in t.dg'range loop
        res(ind)(i) := t.dg(i);
      end loop;
      continue := true;
    end Add_Exponent;
    procedure Add_Exponents is new Visiting_Iterator(Add_Exponent);

  begin
    Add_Exponents(p);
    return res;
  end Create;

  function Create ( p : Complex_Multivariate_Polynomials.Poly )
                  return Integer_Vectors_of_Vectors.Vector is

    use Complex_Multivariate_Polynomials;
    res : Integer_Vectors_of_Vectors.Vector(1..Number_of_Terms(p));
    ind : natural := 0;
    
    procedure Add_Exponent ( t : in Term; continue : out boolean ) is
    begin
      ind := ind + 1;
      res(ind) := new Integer_Vectors.Vector(t.dg'range);
      for i in t.dg'range loop
        res(ind)(i) := t.dg(i);
      end loop;
      continue := true;
    end Add_Exponent;
    procedure Add_Exponents is new Visiting_Iterator(Add_Exponent);

  begin
    Add_Exponents(p);
    return res;
  end Create;

  function Create ( p : Poly_Sys ) return Exponent_Vectors_Array is

    res : Exponent_Vectors_Array(p'range);

  begin
    for i in p'range loop
      declare
        cpi : constant Integer_Vectors_of_Vectors.Vector := Create(p(i));
      begin
        res(i) := new Integer_Vectors_of_Vectors.Vector(cpi'range);
        for j in cpi'range loop
          res(i)(j) := cpi(j);
        end loop;
      end;  -- a detour for GNAT 3.07
     -- res(i) := new Integer_Vectors_of_Vectors.Vector'(Create(p(i)));
    end loop;
    return res;
  end Create;

  function Create ( p : Laur_Sys ) return Exponent_Vectors_Array is

    res : Exponent_Vectors_Array(p'range);

  begin
    for i in p'range loop
      declare
        cpi : constant Integer_Vectors_of_Vectors.Vector := Create(p(i));
      begin
        res(i) := new Integer_Vectors_of_Vectors.Vector(cpi'range);
        for j in cpi'range loop
          res(i)(j) := cpi(j);
        end loop;
      end;  -- a detour for GNAT 3.07
     -- res(i) := new Integer_Vectors_of_Vectors.Vector'(Create(p(i)));
    end loop;
    return res;
  end Create;

-- SELECTOR :

  function Position ( ev : Integer_Vectors_of_Vectors.Vector;
                      v : Integer_Vectors.Vector ) return integer is
  begin
    for i in ev'range loop
      if Integer_Vectors.Equal(ev(i).all,v)
       then return i;
      end if;
    end loop;
    return ev'last+1;
  end Position;

-- DESTRUCTORS :

  procedure Clear ( v : in out Exponent_Vectors_Array ) is
  begin
    for i in v'range loop
      Integer_Vectors_of_Vectors.Clear(v(i));
    end loop;
  end Clear;

end Exponent_Vectors;
