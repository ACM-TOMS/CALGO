with Integer_Vectors,Integer_Matrices;  use Integer_Vectors,Integer_Matrices;
with Integer_Linear_System_Solvers;     use Integer_Linear_System_Solvers;
with Sets_of_Unknowns;                  use Sets_of_Unknowns;
with Degrees_in_Sets_of_Unknowns;       use Degrees_in_Sets_of_Unknowns;
with Complex_Multivariate_Polynomials;  use Complex_Multivariate_Polynomials;

package body m_Homogeneous_Bezout_Numbers is

-- UTILITIES :

  function Create ( p : Poly_Sys ) return Set is

  -- DESCRIPTION :
  --   Returns the set of the unknowns of the polynomial system p.

    s : Set := Create(p'length);

  begin
    for i in p'range loop
      Add(s,i);
    end loop;
    return s;
  end Create;

  function Cardinalities ( z : Partition ) return Vector is

  -- DESCRIPTION :
  --   Returns a vector which contains the cardinality of each set.

    res : Vector(z'range);

  begin
    for i in z'range loop
      res(i) := Extent_Of(z(i));
    end loop;
    return res;
  end Cardinalities;

-- TARGET ROUTINES :

  function Total_Degree ( p : Poly_Sys ) return natural is

    d : natural := 1;

  begin
    for i in p'range loop
      d := d*Degree(p(i));
    end loop;
    return d;
  end Total_Degree;

  function Bezout_Number ( p : Poly_Sys; z : Partition ) return natural is

    k : constant vector := Cardinalities(z);
    d : constant matrix := Degree_Table(p,z);

  begin
    return Per(d,k);  -- returns the permanent of the degree table
  end Bezout_Number;

  function Bezout_Number ( p : Poly_Sys; z : Partition; max : natural )
                         return natural is

  -- DESCRIPTION :
  --   Stops when the Bezout number becomes bigger than max.

    k : constant vector := Cardinalities(z);
    d : constant matrix := Degree_Table(p,z);

  begin
    return Per(d,k,max);  -- returns the permanent of the degree table
  end Bezout_Number;

  function Bezout_Number ( p : Poly_Sys ) return natural is

    s : Set := Create(p);
    res : natural := Total_Degree(p);

    procedure Evaluate ( z : in Partition; cont : out boolean ) is
      b : constant natural := Bezout_Number(p,z,res);
    begin
      if b < res
       then res := b;
      end if;
      cont := true;
    end Evaluate;
    procedure Evaluate_Partitions is new Generate_Partitions(Evaluate);

  begin
    Evaluate_Partitions(s);
    Clear(s);
    return res;
  end Bezout_Number;

  function Bezout_Number ( max : natural; p : Poly_Sys ) return natural is

    s : Set := Create(p);
    res : natural := Total_Degree(p);
    cnt : natural := 0;

    procedure Evaluate ( z : in Partition; cont : out boolean ) is
      b : constant natural := Bezout_Number(p,z,res);
    begin
      if b < res
       then res := b;
      end if;
      cnt := cnt + 1;
      if cnt < max
       then cont := true;
       else cont := false;
      end if;
    end Evaluate;
    procedure Evaluate_Partitions is new Generate_Partitions(Evaluate);

  begin
    Evaluate_Partitions(s);
    Clear(s);
    return res;
  end Bezout_Number;

  function Bezout_Number ( p : Poly_Sys; min : natural ) return natural is

    s : Set := Create(p);
    res : natural := Total_Degree(p);

    procedure Evaluate ( z : in Partition; cont : out boolean ) is
      b : constant natural := Bezout_Number(p,z,res);
    begin
      if b < res
       then res := b;
      end if;
      if res > min
       then cont := true;
       else cont := false;
      end if;
    end Evaluate;
    procedure Evaluate_Partitions is new Generate_Partitions(Evaluate);

  begin
    Evaluate_Partitions(s);
    Clear(s);
    return res;
  end Bezout_Number;

  function Bezout_Number ( max : natural; p : Poly_Sys; min : natural )
                         return natural is

    s : Set := Create(p);
    res : natural := Total_Degree(p);
    cnt : natural := 0;   

    procedure Evaluate ( z : in Partition; cont : out boolean ) is
      b : constant natural := Bezout_Number(p,z,res);
    begin
      if b < res
       then res := b;
      end if;
      cnt := cnt + 1;
      if cnt < max and then res > min
       then cont := true;
       else cont := false;
      end if;
    end Evaluate;
    procedure Evaluate_Partitions is new Generate_Partitions(Evaluate);

  begin
    Evaluate_Partitions(s);
    Clear(s);
    return res;
  end Bezout_Number;

  procedure Bezout_Number 
               ( p : in Poly_Sys; b,m : out natural; z : in out Partition ) is

    s : Set := Create(p);
    tdg : constant natural := Total_Degree(p);
    res : natural := tdg;

    procedure Evaluate ( nz : in Partition; cont : out boolean ) is
      nb : constant natural := Bezout_Number(p,nz,res);
    begin
      if nb < res
       then res := nb;
            m := nz'length; Clear(z);
            z(nz'range) := Create(nz);
      end if;
      cont := true;
    end Evaluate;
    procedure Evaluate_Partitions is new Generate_Partitions(Evaluate);

  begin
    Evaluate_Partitions(s);
    if res = tdg
     then m := 1; z(1) := s;
     else Clear(s);
    end if;
    b := res;
  end Bezout_Number;
 
  procedure Bezout_Number 
               ( max : in natural; p : in Poly_Sys; b,m : out natural;
                 z : in out Partition ) is

    s : Set := Create(p);
    tdg : constant natural := Total_Degree(p);
    res : natural := tdg;
    cnt : natural := 0;

    procedure Evaluate ( nz : in Partition; cont : out boolean ) is
      nb : constant natural := Bezout_Number(p,nz,res);
    begin
      if nb < res
       then res := nb;
            m := nz'length; Clear(z);
            z(nz'range) := Create(nz);
      end if;
      cnt := cnt + 1;
      if cnt < max
       then cont := true;
       else cont := false;
      end if;
    end Evaluate;
    procedure Evaluate_Partitions is new Generate_Partitions(Evaluate);

  begin
    Evaluate_Partitions(s);
    if res = tdg
     then m := 1; z(1) := s;
     else Clear(s);
    end if;
    b := res;
  end Bezout_Number;

  procedure Bezout_Number
               ( p : in Poly_Sys; min : in natural; b,m : out natural;
                 z : in out Partition ) is

    s : Set := Create(p);
    tdg : constant natural := Total_Degree(p);
    res : natural := tdg;

    procedure Evaluate ( nz : in Partition; cont : out boolean ) is
      nb : constant natural := Bezout_Number(p,nz,res);
    begin
      if nb < res
       then res := nb;
            m := nz'length; Clear(z);
            z(nz'range) := Create(nz);
      end if;
      if res > min
       then cont := true;
       else cont := false;
      end if;
    end Evaluate;
    procedure Evaluate_Partitions is new Generate_Partitions(Evaluate);

  begin
    Evaluate_Partitions(s);
    if res = tdg
     then m := 1; z(1) := s;
     else Clear(s);
    end if;
    b := res;
  end Bezout_Number;

  procedure Bezout_Number
               ( max : in natural; p : in Poly_Sys; min : in natural;
                 b,m : out natural; z : in out Partition ) is

    s : Set := Create(p);
    tdg : constant natural := Total_Degree(p);
    res : natural := tdg;
    cnt : natural := 0;

    procedure Evaluate ( nz : in Partition; cont : out boolean ) is
      nb : constant natural := Bezout_Number(p,nz,res);
    begin
      if nb < res
       then res := nb;
            m := nz'length; Clear(z);
            z(nz'range) := Create(nz);
      end if;
      cnt := cnt + 1;
      if cnt < max and then res > min
       then cont := true;
       else cont := false;
      end if;
    end Evaluate;
    procedure Evaluate_Partitions is new Generate_Partitions(Evaluate);

  begin
    Evaluate_Partitions(s);
    if res = tdg
     then m := 1; z(1) := s;
     else Clear(s);
    end if;
    b := res;
  end Bezout_Number;

  function Evaluate ( z : partition; m : natural; p : Poly_Sys )
                    return natural is

    n : natural := p'length;
    d : constant matrix := Degree_Table(p,z);

    function Bezout_number 
                 ( n,m : natural; z: partition; d : matrix ) return natural is

    -- DESCRIPTION : the Bezout number is computed

      type boolean_array is array ( positive range <> ) of boolean;
      Ii,Iacc : boolean_array(1..n) := (1..n => false);
      b,b_mult : natural;

      procedure column ( j,start,number : in natural;
                         Ii,Iacc : in out boolean_array );

      -- DESCRIPTION : the computation of a term coming from a column
      --               of the degree table;

      procedure row ( j : in natural; Ii,Iacc : in out boolean_array );

      -- DESCRIPTION : the computation of a row in the degree table

      procedure column ( j,start,number : in natural;
                         Ii,Iacc : in out boolean_array ) is
      begin
        if number > (n - start + 1)
         then return;
         elsif number = 0
             then row(j,Ii,Iacc);
             else for i in start..n loop
                    if not Ii(i)
                     then Ii(i) := true; Iacc(i) := true;
                          column(j,i+1,number-1,Ii,Iacc);
                          Ii(i) := false; Iacc(i) := false;
                    end if;
                  end loop;
        end if;
      end column;

      procedure row ( j : in natural; Ii,Iacc : in out boolean_array ) is
        temp : natural := 1;
        Iacc1 : boolean_array(1..n) := (1..n => false);
      begin
        for k in 1..n loop
          if Iacc(k)
           then temp := temp * d(k,j);
          end if;
        end loop;
        if (j /= m) and (temp /= 0)
         then b_mult := b_mult * temp;
              column(j+1,1,Extent_Of(z(j+1)),Ii,Iacc1);
              b_mult := b_mult / temp;
         elsif j = m
             then temp := temp * b_mult;
                   b := b + temp;
        end if;
      end row;

    begin
      b := 0; b_mult := 1;
      column(1,1,Extent_Of(z(1)),Ii,Iacc);
      return b;
    end Bezout_number;

  begin
    return Bezout_number(n,m,z,d);
  end Evaluate;

  procedure PB ( p : in Poly_Sys; b,m : out natural; z : in out Partition ) is

    n : natural := p'length;
    wb,b_min : natural;
    wz,z_min : partition(1..n);
    wm,m_min : natural := 0;

    procedure pcopy ( p1,p2 : in out partition; p1n,p2n : in out natural ) is
    -- DESCRIPTION : the partition p1 is copied to p2
    begin
      Clear(p2); p2 := Create(p1);
      p2n := p1n;
    end pcopy;

  begin
    b_min := Total_Degree(p);
    for i in 1..n loop
      for k in 1..wm loop
        Add(wz(k),i); 
        wb := Evaluate(wz,wm,p);
        if (k = 1) or else (wb < b_min)
         then pcopy(wz,z_min,wm,m_min);
              b_min := wb;
        end if;
        Remove(wz(k),i);
      end loop;
      wm := wm + 1;
      wz(wm) := Create(n);
      Add(wz(wm),i); 
      wb := Evaluate(wz,wm,p);
      if wb < b_min
       then pcopy(wz,z_min,wm,m_min);
            b_min := wb;
      end if;
      pcopy(z_min,wz,m_min,wm);
    end loop;
    b := b_min;
    m := m_min;
    pcopy(z_min,z,m_min,wm);
    Clear(wz);
  end PB;

end m_Homogeneous_Bezout_Numbers;
