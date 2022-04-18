with Complex_Numbers,Integer_Vectors;  use Complex_Numbers,Integer_Vectors;
with Permutations,Permute_Operations;  use Permutations,Permute_Operations;
with Random_Number_Generators;
with Complex_Multivariate_Laurent_Polynomials;

function Symmetric_Randomize ( p : Laur_Sys; v,w : List_of_Permutations )
                             return Laur_Sys is

  res : Laur_Sys(p'range);
  
  use Complex_Multivariate_Laurent_Polynomials;

  procedure Symmetric_Randomize_Terms ( index : in natural; py : in Poly ) is

    tpy : Term;

    procedure Permute_and_Randomize ( t : in Term ) is

      tmpv,tmpw : List_of_Permutations;

    begin
      tmpv := v;  tmpw := w;
      while not Is_Null(tmpv) loop
        declare
          permt : Term := Permutation(Head_Of(tmpv).all)*t;
          indw : natural := Head_Of(tmpw)(index);
        begin
          if Coeff(res(indw),permt.dg) = CMPLX(0.0)
           then Plus_Term(res(indw),permt);
          end if;
          Clear(permt);
        end;
        tmpv := Tail_Of(tmpv); 
        tmpw := Tail_Of(tmpw);
      end loop;
    end Permute_and_Randomize;

    procedure Pick_Term ( t : in Term; cont : out boolean ) is
    begin
      if Coeff(res(index),t.dg) = CMPLX(0.0)
       then Copy(t,tpy);
            tpy.cf := Random_Number_Generators.Random1;
            cont := false;
       else cont := true;
      end if;
    end Pick_Term;
    procedure Pick_A_Term is new Visiting_Iterator(Pick_Term);

  begin
    tpy.cf := CMPLX(0.0);
    Pick_A_Term(py);
    if tpy.cf /= CMPLX(0.0)
     then Permute_and_Randomize(tpy);
          Clear(tpy);
    end if;
  end Symmetric_Randomize_Terms;

begin
  res := (res'range => Null_Poly);
  for k in res'range loop
    while Number_of_Terms(res(k)) < Number_of_Terms(p(k)) loop
      Symmetric_Randomize_Terms(k,p(k));
    end loop;
  end loop;
  return res;
end Symmetric_Randomize;
