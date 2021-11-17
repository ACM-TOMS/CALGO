with Floating_Point_Numbers;        use Floating_Point_Numbers;
with Complex_Numbers;               use Complex_Numbers;
with Natural_Vectors;
with Cayley_Embedding;              use Cayley_Embedding;
with Mixed_Volume_Computation;      use Mixed_Volume_Computation;

package body Minkowski_Polynomials is

  function Minkowski_Polynomial ( n,r : natural ) return Poly is

    res : Poly := Null_Poly;
    acc : Degrees := new Natural_Vectors.Vector'(1..r => 0);

    procedure Generate_Monomials
                  ( k,sum : in natural; deg : in out Degrees ) is

    -- DESCRIPTION :
    --   Generates all exponent vectors whose sum equals n.

      t : Term;

    begin
      if k = r
       then t.cf := CMPLX(1.0);
            deg(r) := n-sum;
            t.dg := deg;
            Plus_Term(res,t);
       else for i in 0..(n-sum) loop
              deg(k) := i;
              Generate_Monomials(k+1,sum+i,deg);
            end loop;
      end if;
    end Generate_Monomials;

  begin
    Generate_Monomials(1,0,acc);
    Natural_Vectors.Clear(Natural_Vectors.Link_to_Vector(acc));
    return res;
  end Minkowski_Polynomial;

  function Convert ( dg : Degrees ) return Integer_Vectors.Vector is

  -- DESCRIPTION :
  --   Converts the degrees vector to a vector with integer numbers.

    res : Integer_Vectors.Vector(dg'range);

  begin
    for i in res'range loop
      res(i) := dg(i);
    end loop;
    return res;
  end Convert;

  procedure Minkowski_Polynomial
                  ( p : in out Poly; t : in Triangulation; n : in natural;
                    mix : in Vector; mixsub : out Mixed_Subdivision ) is

    procedure Coefficient_Volume
                  ( submix : in Vector; sub : in Mixed_Subdivision;
                    vol : out natural ) is
    begin
      vol := Mixed_Volume(n,submix,sub);
    end Coefficient_Volume;
    procedure Coefficient_Volumes is
      new Minkowski_Polynomial_Subdivisions(Coefficient_Volume);

  begin
    Coefficient_Volumes(p,t,n,mix,mixsub);
  end Minkowski_Polynomial;

  procedure Minkowski_Polynomial_Subdivisions
                 ( p : in out Poly; t : in Triangulation; n : in natural;
                   mix : in Vector; mixsub : out Mixed_Subdivision ) is

    procedure Coefficient_Volume ( tt : in out Term; cont : out boolean ) is

      wrkmix : Vector(mix'range) := Convert(tt.dg);
      wrksub : Mixed_Subdivision := Extract_Mixed_Cells(n,wrkmix,t);
      vol : natural;

    begin
      Deflate(n,wrksub);
      Process(wrkmix,wrksub,vol);
      tt.cf := CMPLX(double_float(vol));
      if wrkmix = mix
       then mixsub := wrksub;
       else Deep_Clear(wrksub);
      end if;
      cont := true;
    end Coefficient_Volume;
    procedure Coefficient_Volumes is new Changing_Iterator(Coefficient_Volume);

  begin
    Coefficient_Volumes(p);
  end Minkowski_Polynomial_Subdivisions;

end Minkowski_Polynomials;
