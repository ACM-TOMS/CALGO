with integer_io;                        use integer_io;
with Integer_Vectors,Float_Vectors;     use Integer_Vectors;
with Lists_of_Integer_Vectors;          use Lists_of_Integer_Vectors;
with Arrays_of_Integer_Vector_Lists;    use Arrays_of_Integer_Vector_Lists;
with Power_Lists;                       use Power_Lists;

with Polynomial_to_Laurent_Converters;  use Polynomial_to_Laurent_Converters;
with Complex_Laurent_Polynomial_Systems;
 use Complex_Laurent_Polynomial_Systems;

with Trees_of_Vectors,Volumes;          use Trees_of_Vectors;
with Trees_of_Vectors_io;               use Trees_of_Vectors_io;
with Mixed_Homotopy_Continuation;

with Integer_Lifting_Utilities;         use Integer_Lifting_Utilities;
with Mixed_Volume_Computation;
with Integer_Mixed_Subdivisions;        use Integer_Mixed_Subdivisions;
with Integer_Mixed_Subdivisions_io;     use Integer_Mixed_Subdivisions_io;
with Integer_Polyhedral_Continuation;   use Integer_Polyhedral_Continuation;

package body BKK_Bound_Computations is

  function BKK_by_Implicit_Lifting ( p : Poly_Sys ) return natural is

    n : constant natural := p'length;
    supports : Array_of_Lists(p'range) := Construct_Power_Lists(p);
    bkk : constant natural := Volumes.Mixed_Volume(n,supports);

  begin
    Deep_Clear(supports);
    return bkk;
  end BKK_by_Implicit_Lifting;

  function BKK_by_Implicit_Lifting ( file : file_type; p : Poly_Sys )
                                   return natural is
    n : constant natural := p'length;
    supports : Array_of_Lists(p'range) := Construct_Power_Lists(p);
    tv : Tree_of_Vectors;
    bkk : natural;

  begin
    Volumes.Mixed_Volume(n,supports,tv,bkk);
    new_line(file);
    put_line(file,"The tree of useful directions : "); put(file,tv);
    new_line(file);
    Deep_Clear(supports); Clear(tv);
    return bkk;
  end BKK_by_Implicit_Lifting;

  function BKK_by_Static_Lifting ( p : Poly_Sys ) return natural is

    n : constant natural := p'length;
    supports : Array_of_Lists(p'range) := Construct_Power_Lists(p);
    bkk : constant natural := Mixed_Volume_Computation.Mixed_Volume(n,supports);

  begin
    Deep_Clear(supports);
    return bkk;
  end BKK_by_Static_Lifting;

  function BKK_by_Static_Lifting ( file : file_type; p : Poly_Sys )
                                 return natural is

    n : constant natural := p'length;
    supports : Array_of_Lists(p'range) := Construct_Power_Lists(p);
    bkk : constant natural
        := Mixed_Volume_Computation.Mixed_Volume(file,n,supports);

  begin
    Deep_Clear(supports);
    return bkk;
  end BKK_by_Static_Lifting;

  function Solve_by_Implicit_Lifting ( p : Poly_Sys ) return Solution_List is

    n : constant natural := p'length;
    supports : Array_of_Lists(p'range) := Construct_Power_Lists(p);
    tv : Tree_of_Vectors;
    bkk : natural;
    lp : Laur_Sys(p'range) := Polynomial_to_Laurent_System(p);
    sols : Solution_List;

  begin
    Volumes.Mixed_Volume(n,supports,tv,bkk);   
    Deep_Clear(supports);
    Mixed_Homotopy_Continuation.Solve(Standard_Output,lp,tv,bkk,sols);
    Clear(tv); Clear(lp);
    return sols;
  end Solve_by_Implicit_Lifting;

  function Solve_by_Implicit_Lifting ( file : file_type; p : Poly_Sys )
                                     return Solution_List is

    n : constant natural := p'length;
    supports : Array_of_Lists(p'range) := Construct_Power_Lists(p);
    tv : Tree_of_Vectors;
    bkk : natural;
    lp : Laur_Sys(p'range) := Polynomial_to_Laurent_System(p);
    sols : Solution_List;

  begin
    Volumes.Mixed_Volume(n,supports,tv,bkk);
    Deep_Clear(supports);
    new_line(file);
    put(file,"The BKK bound equals "); put(file,bkk,1); new_line(file);
    put_line(file," with tree of useful directions : "); put(file,tv);
    new_line(file);
    Mixed_Homotopy_Continuation.Solve(file,lp,tv,bkk,sols);
    Clear(tv); Clear(lp);
    return sols;
  end Solve_by_Implicit_Lifting;
  
  function Solve_by_Static_Lifting ( p : Poly_Sys ) return Solution_List is

    n : constant natural := p'length;
    supports : Array_of_Lists(p'range) := Construct_Power_Lists(p);
    mv : natural;
    sols : Solution_List;
    mix,per : Link_to_Vector;
    mixsub : Mixed_Subdivision;
    permp : Poly_Sys(p'range);
    lp,pp : Laur_Sys(p'range);

    use Mixed_Volume_Computation;

  begin
    Compute_Mixture(supports,mix,per);
    permp := Permute(p,per); Clear(per);
    declare
      lifted : Array_of_Lists(mix'range);
    begin
      Mixed_Volume(n,mix.all,supports,lifted,mixsub,mv);
      pp := Polynomial_to_Laurent_System(p);
      lp := Perform_Lifting(n,mix.all,lifted,pp);
      Mixed_Solve(lp,mix.all,mixsub,sols);
      Deep_Clear(lifted);
    end;
    Clear(pp); Clear(lp);
    Deep_Clear(supports);
    Deep_Clear(mixsub);
    Clear(mix);
    return sols;
  end Solve_by_Static_Lifting;

  function Solve_by_Static_Lifting ( file : file_type; p : Poly_Sys )
                                   return Solution_List is

    n : constant natural := p'length;
    supports : Array_of_Lists(p'range) := Construct_Power_Lists(p);
    mv : natural;
    sols : Solution_List;
    mix,per : Link_to_Vector;
    mixsub : Mixed_Subdivision;
    permp : Poly_Sys(p'range);
    lp,pp : Laur_Sys(p'range);

    use Mixed_Volume_Computation;

  begin
    Compute_Mixture(supports,mix,per);
    permp := Permute(p,per); Clear(per);
    declare
      lifted : Array_of_Lists(mix'range);
    begin
      Mixed_Volume(file,n,mix.all,supports,lifted,mixsub,mv);
      put_line(file,"THE MIXED SUBDIVISION :"); 
      put(file,n,mix.all,mixsub);
      pp := Polynomial_to_Laurent_System(p);
      lp := Perform_Lifting(n,mix.all,lifted,pp);
      new_line(file); put_line(file,"POLYHEDRAL HOMOTOPY CONTINUATION");
      Mixed_Solve(file,lp,mix.all,mixsub,sols);
      Deep_Clear(lifted);
    end;
    Clear(pp); Clear(lp);
    Deep_Clear(supports);
    Deep_Clear(mixsub);
    Clear(mix);
    return sols;
  end Solve_by_Static_Lifting;

end BKK_Bound_Computations;
