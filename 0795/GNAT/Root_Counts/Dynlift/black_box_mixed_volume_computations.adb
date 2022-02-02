with Lists_of_Integer_Vectors;           use Lists_of_Integer_Vectors; 
with Power_Lists,Vertices;               use Power_Lists,Vertices;
with Integer_Lifting_Utilities;          use Integer_Lifting_Utilities;
with Cayley_Trick;                       use Cayley_Trick;
with Triangulations;                     use Triangulations;
with Dynamic_Triangulations;             use Dynamic_Triangulations;
with Triangulations_and_Subdivisions;    use Triangulations_and_Subdivisions;
with Flatten_Mixed_Subdivisions;         use Flatten_Mixed_Subdivisions;
with Mixed_Volume_Computation;           use Mixed_Volume_Computation;
with Complex_Numbers;                    use Complex_Numbers;
with Complex_Vectors;
with Complex_Vectors_of_Vectors;
with Exponent_Vectors;                   use Exponent_Vectors;
with Polynomial_to_Laurent_Converters;   use Polynomial_to_Laurent_Converters;
with Laurent_to_Polynomial_Converters;   use Laurent_to_Polynomial_Converters;
with Complex_Multivariate_Laurent_Polynomials;
 use Complex_Multivariate_Laurent_Polynomials;
with Complex_Laurent_Polynomial_Systems;
 use Complex_Laurent_Polynomial_Systems;
with Laurent_Jacobi_Matrices;            use Laurent_Jacobi_Matrices;
with Polynomial_Randomizers;             use Polynomial_Randomizers;
with Continuation_Parameters;
with Integer_Polyhedral_Continuation;    use Integer_Polyhedral_Continuation;

package body Black_Box_Mixed_Volume_Computations is

  procedure Black_Box_Mixed_Volume_Computation
                 ( p : in Poly_Sys; mix : out Link_to_Vector;
                   lifsup : out Link_to_Array_of_Lists;
                   mixsub : out Mixed_Subdivision; mv : out natural ) is

    n : constant natural := p'length;
    supports : Array_of_Lists(p'range) := Construct_Power_Lists(p);
    verpts : Array_of_Lists(p'range);
    tmix,perms : Integer_Vectors.Link_to_Vector;
    tmixsub : Mixed_Subdivision;

    procedure Collect_Flattening ( t : in Triangulation; l : List ) is

    -- DESCRIPTION :
    --   Updates the subdivision tmixsub with the flattened cells.
    --   The triangulation on entry contains the whole triangulation,
    --   not just the new cells.

      cells : Mixed_Subdivision;

    begin
      if Is_Null(tmixsub)
       then cells := Deep_Create(n,t);
       else cells := Non_Flat_Deep_Create(n,t);
            Construct(Head_Of(tmixsub),cells);
      end if;
      Flatten(cells);
      tmixsub := cells;
    end Collect_Flattening;
    procedure C_Dynamic_Lifting is
      new Dynamic_Triangulations.Dynamic_Lifting_with_Flat(Collect_Flattening);

  begin
    for i in supports'range loop
      verpts(i) := Vertex_Points(supports(i));
    end loop;
    Compute_Mixture(verpts,tmix,perms);
    declare
      pts,lifted : Array_of_Lists(tmix'range);
      last : List;
      t : Triangulation;
      nt : natural;
      lastcells : Mixed_Subdivision;
    begin
      if tmix'length = 1
       then C_Dynamic_Lifting(verpts(1),false,false,0,lifted(1),last,t);
            if Is_Null(tmixsub)
             then tmixsub := Deep_Create(n,t);
             else lastcells := Non_Flat_Deep_Create(n,t);
                  Construct(Head_Of(tmixsub),lastcells);
                  tmixsub := lastcells;
            end if;
            Clear(t);
            Mixed_Volume(n,tmix.all,tmixsub,mv);
       elsif tmix'length <= n/2
           then
             pts := Typed_Lists(tmix.all,verpts);
             Dynamic_Cayley(n,tmix.all,pts,false,false,0,lifted,tmixsub,nt);
             Mixed_Volume(n,tmix.all,tmixsub,mv);
           else
             Mixed_Volume(n,tmix.all,verpts,lifted,tmixsub,mv);
      end if;
      lifsup := new Array_of_Lists'(lifted);
    end;
    mix := tmix; mixsub := tmixsub;
  end Black_Box_Mixed_Volume_Computation;

  procedure Black_Box_Polyhedral_Continuation
                 ( p : in Poly_Sys; mix : in Vector;
                   lifsup : in Array_of_Lists; mixsub : in Mixed_Subdivision;
                   q : in out Poly_Sys; qsols : in out Solution_List ) is

    n : constant natural := p'length;
    lq,llq : Laur_Sys(p'range);
    h : Eval_Coeff_Laur_Sys(q'range);
    c : Complex_Vectors_of_Vectors.Vector(h'range);
    e : Exponent_Vectors_Array(h'range);
    j : Eval_Coeff_Jacobi(h'range,h'first..h'last+1);
    m : Mult_Factors(j'range(1),j'range(2));

  begin
    q := Complex_Randomize1(p);
    lq := Polynomial_to_Laurent_System(q);
    llq := Perform_Lifting(n,mix,lifsup,lq);
    Clear(lq); Clear(q);
    lq := Eval(llq,CMPLX(1.0),n+1);
    q := Laurent_to_Polynomial_System(lq);
    Continuation_Parameters.Tune(0);
   -- Mixed_Solve(llq,mix,mixsub,qsols);   too expensive !!!!
    h := Create(lq);
    for i in c'range loop
      declare
        coeff_lq : constant Complex_Vectors.Vector := Coeff(lq(i));
      begin
        c(i) := new Complex_Vectors.Vector(coeff_lq'range);
        for k in coeff_lq'range loop
          c(i)(k) := coeff_lq(k);
        end loop;
      end;
    end loop;
    e := Create(lq);
    Create(lq,j,m);
    Mixed_Solve(llq,lifsup,h,c,e,j,m,mix,mixsub,qsols);
    Set_Continuation_Parameter(qsols,CMPLX(0.0));
    Clear(lq); Clear(llq);
    Clear(h); Clear(j); Clear(m);
    Complex_Vectors_of_Vectors.Clear(c);
  end Black_Box_Polyhedral_Continuation;

end Black_Box_Mixed_Volume_Computations;
