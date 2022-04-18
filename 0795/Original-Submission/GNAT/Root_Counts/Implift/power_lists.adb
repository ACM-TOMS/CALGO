with Integer_Vectors;           use Integer_Vectors;

package body Power_Lists is

  function Construct_Power_List ( p : Complex_Multivariate_Polynomials.Poly )
                                return List is

    use Complex_Multivariate_Polynomials;
    res,res_last : List;

    procedure Fill_In_Term (t : in Term; cont : out boolean) is

      h : Link_to_Vector;

    begin
      h := new Integer_Vectors.Vector(t.dg'range);
      for j in h'range loop
        h(j) := t.dg(j);
      end loop;
      Append(res,res_last,h);
      cont := true;
    end Fill_In_Term;
    procedure Fill_In is new Visiting_Iterator(Fill_In_Term);

  begin
    Fill_In(p);
    return res;
  end Construct_Power_List;

  function Construct_Power_List
               ( p : Complex_Multivariate_Laurent_Polynomials.Poly )
	       return List is

    use Complex_Multivariate_Laurent_Polynomials;
    res,res_last : List;

    procedure Fill_In_Term (t : in Term; cont : out boolean) is

      h : Link_to_Vector;

    begin
      h := new Integer_Vectors.Vector(t.dg'range);
      for j in h'range loop
        h(j) := t.dg(j);
      end loop;
      Append(res,res_last,h);
      cont := true;
    end Fill_In_Term;
    procedure Fill_In is new Visiting_Iterator(Fill_In_Term);

  begin
    Fill_In(p);
    return res;
  end Construct_Power_List;

  function Select_Terms 
              ( p : Complex_Multivariate_Polynomials.Poly; l : List )
              return Complex_Multivariate_Polynomials.Poly is

    use Complex_Multivariate_Polynomials;
    res : Poly := Null_Poly;

    procedure Select_Term ( t : in Term; cont : out boolean ) is

      v : Integer_Vectors.Vector(t.dg'range);

    begin
      for i in v'range loop
        v(i) := t.dg(i);
      end loop;
      if Is_In(l,v)
       then Plus_Term(res,t);
      end if;
      cont := true;
    end Select_Term;
    procedure Select_Poly is new Visiting_Iterator ( Select_Term );

  begin
    Select_Poly(p);
    return res;
  end Select_Terms;

  function Select_Terms 
              ( p : Complex_Multivariate_Laurent_Polynomials.Poly; l : List )
              return Complex_Multivariate_Laurent_Polynomials.Poly is

    use Complex_Multivariate_Laurent_Polynomials;
    res : Poly := Null_Poly;

    procedure Select_Term ( t : in Term; cont : out boolean ) is
    begin
      if Is_In(l,t.dg.all)
       then Plus_Term(res,t);
      end if;
      cont := true;
    end Select_Term;
    procedure Select_Poly is new Visiting_Iterator ( Select_Term );

  begin
    Select_Poly(p);
    return res;
  end Select_Terms;

  function Construct_Power_Lists ( p : Poly_Sys ) return Array_of_Lists is

    res : Array_of_Lists(p'range);

  begin
    for i in p'range loop
      res(i) := Construct_Power_List(p(i));
    end loop;
    return res;
  end Construct_Power_Lists;

  function Construct_Power_Lists ( p : Laur_Sys ) 
	        	         return Array_of_Lists is

    res : Array_of_Lists(p'range);

  begin
    for i in p'range loop
      res(i) := Construct_Power_List(p(i));
    end loop;
    return res;
  end Construct_Power_Lists;

  function Select_Terms ( p : Poly_Sys; al : Array_of_Lists ) 
                        return Poly_Sys is

    res : Poly_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Select_Terms(p(i),al(i));
    end loop;
    return res;
  end Select_Terms;

  function Select_Terms ( p : Laur_Sys; al : Array_of_Lists ) 
                        return Laur_Sys is

    res : Laur_Sys(p'range);

  begin
    for i in p'range loop
      res(i) := Select_Terms(p(i),al(i));
    end loop;
    return res;
  end Select_Terms;

end Power_Lists;
