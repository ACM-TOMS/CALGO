with Floating_Point_Numbers;            use Floating_Point_Numbers;
with Complex_Numbers,Complex_Vectors;   use Complex_Numbers,Complex_Vectors;
with Complex_Matrices;                  use Complex_Matrices;
with Complex_Linear_System_Solvers;     use Complex_Linear_System_Solvers;

with Arrays_of_Integer_Vector_Lists;    use Arrays_of_Integer_Vector_Lists;
with Natural_Vectors,Volumes;           use Volumes;

with Power_Lists,Substitutors;          use Power_Lists,Substitutors;
with Complex_Multivariate_Polynomials;
with Complex_Laurent_Polynomial_Systems;
 use Complex_Laurent_Polynomial_Systems;
with Polynomial_to_Laurent_Converters;  use Polynomial_to_Laurent_Converters;
with Polynomial_Randomizers;            use Polynomial_Randomizers;

with Mixed_Homotopy_Continuation;       use Mixed_Homotopy_Continuation;
with Set_Structure;
with Random_Product_System;
with Construct_Random_Product_Start_System;

package body Set_Structures_and_Volumes is

-- AUXILIARIES :

  procedure Build_RPS ( k,n : in natural; p : in Poly_Sys ) is

  -- DESCRIPTION :
  --   If the set structure is empty, then a set structure will
  --   be constructed for the first k polynomials of p;
  --   the data in Random_Product_System will be build.

  begin
    if Set_Structure.Empty
     then declare
	    tmp : Poly_Sys(p'range);
	    use Complex_Multivariate_Polynomials;
	    t : Term;
	    cnst : Poly;
          begin
	    tmp(1..k) := p(1..k);
	    t.dg := new Natural_Vectors.Vector'(1..n => 0);
	    t.cf := CMPLX(1.0);
	    cnst := Create(t);
	    tmp(k+1..n) := (k+1..n => cnst);
	    Construct_Random_Product_Start_System.Build_Set_Structure(tmp);
	    Natural_Vectors.Clear(Natural_Vectors.Link_to_Vector(t.dg));
            Clear(cnst);
          end;
    end if;
    Random_Product_System.Init(n);
    Construct_Random_Product_Start_System.Build_Random_Product_System(n);
  end Build_RPS;

  procedure Build_Elimination_Matrix 
	      ( k,n : in natural; ind : in Natural_Vectors.Vector;
	        m : in out matrix; pivots : in out Natural_Vectors.Vector;
		degenerate : out boolean ) is

  -- DESCRIPTION :
  --   builds the elimination matrix m.

  -- ON ENTRY :
  --   k            the number of unknowns to be eliminated;
  --   n            the number of polynomials and unknowns;
  --   ind          entries in Random_Product_System;
  --                ind(l) indicates lth hyperplane;

  -- ON RETURN :
  --   m            contains all k hyperplanes to use for elimination;
  --   pivots       if not degenerate, then m(i,pivots(i)) /= 0;
  --   degenerate   is true if the first k hyperplanes were inconsistent.

    degen : boolean;

  begin
   -- BUILD MATRIX :
    for i in ind'range loop
      declare
	h : vector(0..n) := Random_Product_System.Get_Hyperplane(i,ind(i));
      begin
        for j in 1..n loop
  	  m(i,j) := h(j);
        end loop;
        m(i,n+1) := h(0);
      end;
    end loop;
    diagonalize(m,k,n+1);
   -- CHECK IF DEGENERATE :
    degen := false;
    declare
      eps : constant double_float := 10.0**(-10);
    begin
      for i in pivots'range loop
        for j in 1..n loop
  	  if modulus(m(i,j)) > eps
	   then pivots(i) := j;
          end if;
	  exit when (pivots(i) /= 0);
        end loop;
        degen := (pivots(i) = 0);
        exit when degen;
      end loop;
    end;
    degenerate := degen;
  end Build_Elimination_Matrix;

  procedure Eliminate ( k,n : in natural; p : in Poly_Sys;
  		        m : in matrix; pivots : in Natural_Vectors.Vector;
	                q : in out Poly_Sys ) is
  -- DESCRIPTION :
  --   eliminates k unknowns of the polynomial system p.

  -- ON ENTRY :
  --   k            the number of unknowns to be eliminated;
  --   n            the number of polynomials and unknowns;
  --   p            a polynomial system with random coefficients;
  --   m            the diagonalized matrix, to be used for elimination,
  --                for i in m'range(1) : m(i,pivots(i)) /= 0;
  --   pivots       a vector indicating nonzero entries in m.

  -- ON RETURN :
  --   q            the reduced system.

  begin
    Clear(q); Copy(p,q);
    for i in pivots'range loop
      declare
        h : vector(0..n-i+1);
        piv : natural := pivots(i)-i+1;
      begin
        h(0) := m(i,n+1);
        h(1..n-i+1) := (1..n-i+1 => CMPLX(0.0));
        h(piv) := m(i,pivots(i));
        for j in piv+1..n-i+1 loop
          h(j) := m(i,j+i-1);
        end loop;
        Substitute(piv,h,q);
      end;
    end loop;
  end Eliminate;

  procedure Update ( sols : in out Solution_List; k,n : in natural;
		     m : in matrix; pivots : in Natural_Vectors.Vector ) is

  -- DESCRIPTION :
  --   Based on the elimination matrix m, solution components for x1..xk
  --   can be computed.

  -- ON ENTRY :
  --   sols      a list of solutions, with values for x(n-k+1)..x(n);
  --   k         the number of unknowns that have been eliminated;
  --   n         total number of unknowns;
  --   m         the elimination matrix;
  --   pivots    for i in pivots'range : m(i,pivots(i)) /= 0.

  -- ON RETURN :
  --   sols      the solution list with n-dimensional vectors.

    tmp,res,res_last : Solution_List;

  begin
    tmp := sols;
    while not Is_Null(tmp) loop
      declare
	sol : Solution(n-k) := Head_Of(tmp).all;
	soln : Solution(n);
	ind : natural;
      begin
	soln.t := CMPLX(0.0);
	soln.m := sol.m;
	soln.v := (1..n => CMPLX(0.0));
	for i in 1..k loop
	  soln.v(pivots(i)) := -m(i,n+1);
        end loop;
	ind := 0;
	for i in 1..n loop
	  if soln.v(i) = CMPLX(0.0)
	   then ind := ind + 1;
		soln.v(i) := sol.v(ind);
          end if;
        end loop;
	for i in 1..k loop
	  ind := pivots(i);
	  for j in ind+1..n loop
	    soln.v(ind) := soln.v(ind) - m(i,j)*soln.v(j);
	  end loop;
	  soln.v(ind) := soln.v(ind) / m(i,ind);
        end loop;
	Append(res,res_last,soln);
      end;
      tmp := Tail_Of(tmp);
    end loop;
    Clear(sols); Copy(res,sols);
  end Update;

  procedure Build_Indices ( l,k,n : in natural; p : in Poly_Sys;
			    acc : in out Natural_Vectors.Vector;
			    bb : in out natural ) is

  -- DESCRIPTION :
  --   Builds the indices in the vector acc.

  -- ON ENTRY : 
  --   l           index in the array acc;
  --   k           number of unknowns to eliminate;
  --   n           number of polynomials and unknowns;
  --   p           a polynomial system;
  --   acc(1..l-1) contains entries indicating hyperplanes
  --               in Random_Product_System.

  -- ON RETURN :
  --   acc(1..l)   contains entries indicating hyperplanes
  --               in Random_Product_System;
  --   bb          the mixed Bezout BKK bound.

  begin
    if l > k
     then declare
	    degenerate : boolean;
	    m : matrix(1..k,1..n+1);
	    pivots : Natural_Vectors.Vector(1..k);
          begin
	    pivots := (1..k => 0);
	    Build_Elimination_Matrix(k,n,acc,m,pivots,degenerate);
	    if not degenerate
	     then declare
		    q : Poly_Sys(p'range);
                  begin
		    Eliminate(k,n,p,m,pivots,q);
		    bb := bb + Bezout_BKK(0,n-k,q);
		    Clear(q);
                  end;
            end if;
          end;
     else for j in 1..Random_Product_System.Number_of_Hyperplanes(l) loop
	    acc(l) := j;
	    Build_Indices(l+1,k,n,p,acc,bb);
          end loop;
    end if;
  end Build_Indices;

  procedure Build_Indices ( l,k,n : in natural; p : in Poly_Sys;
			    tv : in Tree_of_Vectors;
			    acc : in out Natural_Vectors.Vector;
			    bb : in out natural ) is

  -- DESCRIPTION :
  --   Builds the indices in the vector acc.

  -- ON ENTRY : 
  --   l           index in the array acc;
  --   k           number of unknowns to eliminate;
  --   n           number of polynomials and unknowns;
  --   p           a polynomial system;
  --   tv          the tree of degrees used to compute the mixed volume;
  --   acc(1..l-1) contains entries indicating hyperplanes
  --               in Random_Product_System.

  -- ON RETURN :
  --   acc(1..l)   contains entries indicating hyperplanes
  --               in Random_Product_System;
  --   bb          the mixed Bezout BKK bound.

  begin
    if l > k
     then declare
	    degenerate : boolean;
	    m : matrix(1..k,1..n+1);
	    pivots : Natural_Vectors.Vector(1..k);
          begin
	    pivots := (1..k => 0);
	    Build_Elimination_Matrix(k,n,acc,m,pivots,degenerate);
	    if not degenerate
	     then declare
		    q : Poly_Sys(p'range);
                  begin
		    Eliminate(k,n,p,m,pivots,q);
		    bb := bb + Bezout_BKK(0,n-k,q,tv);
		    Clear(q);
                  end;
            end if;
          end;
     else for j in 1..Random_Product_System.Number_of_Hyperplanes(l) loop
	    acc(l) := j;
	    Build_Indices(l+1,k,n,p,tv,acc,bb);
          end loop;
    end if;
  end Build_Indices;

  procedure Build_Indices2 ( l,k,n : in natural; p : in Poly_Sys;
			     tv : in out Tree_of_Vectors;
			     acc : in out Natural_Vectors.Vector;
			     bb : in out natural ) is

  -- DESCRIPTION :
  --   Builds the indices in the vector acc.

  -- ON ENTRY : 
  --   l           index in the array acc;
  --   k           number of unknowns to eliminate;
  --   n           number of polynomials and unknowns;
  --   p           a polynomial system;
  --   acc(1..l-1) contains entries indicating hyperplanes
  --               in Random_Product_System.

  -- ON RETURN :
  --   tv          the tree of degrees used to compute the mixed volume;
  --   acc(1..l)   contains entries indicating hyperplanes
  --               in Random_Product_System;
  --   bb          the mixed Bezout BKK bound.

  begin
    if l > k
     then declare
	    degenerate : boolean;
	    m : matrix(1..k,1..n+1);
	    pivots : Natural_Vectors.Vector(1..k);
          begin
	    pivots := (1..k => 0);
	    Build_Elimination_Matrix(k,n,acc,m,pivots,degenerate);
	    if not degenerate
	     then declare
		    q : Poly_Sys(p'range);
		    qtv,tmp : Tree_of_Vectors;
		    mv : natural;
                  begin
		    Eliminate(k,n,p,m,pivots,q);
		    Bezout_BKK(0,n-k,q,qtv,mv);
		    bb := bb + mv;
		    tmp := qtv;
		    while not Is_Null(tmp) loop
		      declare
			nd : node := Head_Of(tmp);
                      begin
			if Is_In(tv,nd.d)
			 then Clear(nd);
			 else Construct(nd,tv);
                        end if;
                      end;
		      tmp := Tail_Of(tmp);
                    end loop;
		    Clear(q);
                  end;
            end if;
          end;
     else for j in 1..Random_Product_System.Number_of_Hyperplanes(l) loop
	    acc(l) := j;
	    Build_Indices2(l+1,k,n,p,tv,acc,bb);
          end loop;
    end if;
  end Build_Indices2;

  procedure Build_Indices
                 ( file : in file_type; p : in Poly_Sys;
                   l,k,n : in natural; acc : in out Natural_Vectors.Vector;
                   bb : in out natural;
                   sols,sols_last : in out Solution_List ) is

  -- DESCRIPTION :
  --   Builds the indices in the vector acc.

  -- ON ENTRY : 
  --   file        file to write intermediate results on;
  --   p           a polynomial system;
  --   l           index in the array acc;
  --   k           number of unknowns to eliminate;
  --   n           number of polynomials and unknowns;
  --   acc(1..l-1) contains entries indicating hyperplanes
  --               in Random_Product_System.

  -- ON RETURN :
  --   acc(1..l)   contains entries indicating hyperplanes
  --               in Random_Product_System;
  --   bb          the mixed Bezout BKK bound;
  --   sols        the solutions of pi;
  --   sols_last   points to the last element of the list sols.

  begin
    if l > k
     then declare
	    degenerate : boolean;
	    m : matrix(1..k,1..n+1);
	    pivots : Natural_Vectors.Vector(1..k);
          begin
	    pivots := (1..k => 0);
	    Build_Elimination_Matrix(k,n,acc,m,pivots,degenerate);
	    if not degenerate
	     then declare
		    q : Poly_Sys(p'range);
                  begin
		    Eliminate(k,n,p,m,pivots,q);
		    declare
		      las : Laur_Sys(q'range);
		      qsols : Solution_List;
		      bkk : natural;
                    begin
		      las := Polynomial_to_Laurent_System(q);
		      Solve(file,las,bkk,qsols);
		      bb := bb + bkk;
		      Update(qsols,k,n,m,pivots);
		      Concat(sols,sols_last,qsols);
                      Clear(las); Shallow_Clear(qsols);
                    end;
		    Clear(q);
                  end;
            end if;
          end;
     else for j in 1..Random_Product_System.Number_of_Hyperplanes(l) loop
	    acc(l) := j;
	    Build_Indices(file,p,l+1,k,n,acc,bb,sols,sols_last);
          end loop;
    end if;
  end Build_Indices;

  procedure Build_Indices
                 ( file : in file_type; p : in Poly_Sys; l,k,n : in natural;
                   acc : in out Natural_Vectors.Vector;
                   tv : in Tree_of_Vectors; bb : in out natural;
                   sols,sols_last : in out Solution_List ) is

  -- DESCRIPTION :
  --   Builds the indices in the vector acc.

  -- ON ENTRY : 
  --   file        file to write intermediate results on;
  --   p           a polynomial system;
  --   l           index in the array acc;
  --   k           number of unknowns to eliminate;
  --   n           number of polynomials and unknowns;
  --   acc(1..l-1) contains entries indicating hyperplanes
  --               in Random_Product_System;
  --   tv          the tree containing useful directions;

  -- ON RETURN :
  --   acc(1..l)   contains entries indicating hyperplanes
  --               in Random_Product_System;
  --   bb          the mixed Bezout BKK bound;
  --   sols        the solutions of pi;
  --   sols_last   points to the last element of the list sols.

  begin
    if l > k
     then declare
	    degenerate : boolean;
	    m : matrix(1..k,1..n+1);
	    pivots : Natural_Vectors.Vector(1..k);
          begin
	    pivots := (1..k => 0);
	    Build_Elimination_Matrix(k,n,acc,m,pivots,degenerate);
	    if not degenerate
	     then declare
		    q : Poly_Sys(p'range);
                  begin
		    Eliminate(k,n,p,m,pivots,q);
		    declare
		      las : Laur_Sys(q'range);
		      qsols : Solution_List;
		      bkk : natural;
                    begin
		      las := Polynomial_to_Laurent_System(q);
		      Solve(file,las,tv,bkk,qsols);
		      bb := bb + bkk;
		      Update(qsols,k,n,m,pivots);
		      Concat(sols,sols_last,qsols);
                      Clear(las); Shallow_Clear(qsols);
                    end;
		    Clear(q);
                  end;
            end if;
          end;
     else for j in 1..Random_Product_System.Number_of_Hyperplanes(l) loop
	    acc(l) := j;
	    Build_Indices(file,p,l+1,k,n,acc,tv,bb,sols,sols_last);
          end loop;
    end if;
  end Build_Indices;

-- TARGET ROUTINES :

  function Bezout_BKK ( k,n : natural; p : Poly_Sys ) return natural is

    res : natural;

  begin
    if k = 0
     then declare
	    adl : Array_of_Lists(p'range);
          begin
	    adl := Construct_Power_Lists(p);
	    res := Mixed_Volume(n,adl);
	    Deep_Clear(adl);
          end;
     elsif k = n
	 then if Set_Structure.Empty
	       then Construct_Random_Product_Start_System.
					      Build_Set_Structure(p);
	            res := Set_Structure.B;
		    Set_Structure.Clear;
               else res := Set_Structure.B;
              end if;
         else Build_RPS(k,n,p);
              declare
                acc : Natural_Vectors.Vector(1..k);
                q : Poly_Sys(1..n-k);
              begin
                acc := (1..k => 0); res := 0;
                for i in q'range loop
                  q(i) := Complex_Randomize1(p(i+k));
                end loop;
                Build_Indices(1,k,n,q,acc,res);
		Clear(q);
              end;
	      Random_Product_System.Clear;
    end if;
    return res;
  end Bezout_BKK;

  function Bezout_BKK ( k,n : natural; p : Poly_Sys; tv : Tree_of_Vectors )
		      return natural is
    res : natural;

  begin
    if k = 0
     then declare
	    adl : Array_of_Lists(p'range);
          begin
	    adl := Construct_Power_Lists(p);
	    res := Mixed_Volume(n,adl,tv);
	    Deep_Clear(adl);
          end;
     elsif k = n
	 then if Set_Structure.Empty
	       then Construct_Random_Product_Start_System.
					      Build_Set_Structure(p);
	            res := Set_Structure.B;
		    Set_Structure.Clear;
               else res := Set_Structure.B;
              end if;
         else Build_RPS(k,n,p);
              declare
                acc : Natural_Vectors.Vector(1..k);
                q : Poly_Sys(1..n-k);
              begin
                acc := (1..k => 0); res := 0;
                for i in q'range loop
                  q(i) := Complex_Randomize1(p(i+k));
                end loop;
                Build_Indices(1,k,n,q,tv,acc,res);
		Clear(q);
              end;
	      Random_Product_System.Clear;
    end if;
    return res;
  end Bezout_BKK;

  procedure Bezout_BKK ( k,n : in natural; p : in Poly_Sys;
			 tv : in out Tree_of_Vectors; bb : out natural ) is
    res : natural;

  begin
    if k = 0
     then declare
	    adl : Array_of_Lists(p'range);
          begin
	    adl := Construct_Power_Lists(p);
	    Mixed_Volume(n,adl,tv,res);
	    Deep_Clear(adl);
          end;
     elsif k = n
	 then if Set_Structure.Empty
	       then Construct_Random_Product_Start_System.
					      Build_Set_Structure(p);
	            res := Set_Structure.B;
		    Set_Structure.Clear;
               else res := Set_Structure.B;
              end if;
         else Build_RPS(k,n,p);
              declare
                acc : Natural_Vectors.Vector(1..k);
                q : Poly_Sys(1..n-k);
              begin
                acc := (1..k => 0); res := 0;
                for i in q'range loop
                  q(i) := Complex_Randomize1(p(i+k));
                end loop;
                Build_Indices2(1,k,n,q,tv,acc,res);
		Clear(q);
              end;
	      Random_Product_System.Clear;
    end if;
    bb := res;
  end Bezout_BKK;

  procedure Mixed_Solve
                 ( file : in file_type; k,n : in natural; p : in Poly_Sys;
                   bb : out natural;
                   g : in out Poly_Sys; sols : in out Solution_List ) is
  begin
    if k = 0
     then declare
            l : Laur_Sys(g'range);
            tmp : Solution_List;
          begin
            l := Polynomial_to_Laurent_System(g);
            Solve(file,l,bb,sols);
            tmp := sols;
            while not Is_Null(tmp) loop
              Head_Of(tmp).t := CMPLX(0.0);
              tmp := Tail_Of(tmp);
            end loop;
            Clear(l);
          end;
     elsif k = n
	 then Construct_Random_Product_Start_System.Construct(p,g,sols);
	      bb := Length_Of(sols);
         else Build_RPS(k,n,p);
              declare
                acc : Natural_Vectors.Vector(1..k);
                q : Poly_Sys(1..n-k);
		rq : Poly_Sys(p'range);
		sols_last : Solution_List;
		res : natural;
              begin
                acc := (1..k => 0); res := 0;
                for i in q'range loop
                  q(i) := g(i+k);
                end loop;
                Build_Indices(file,q,1,k,n,acc,res,sols,sols_last);
		bb := res;
	        rq := Random_Product_System.Polynomial_System;
		g(1..k) := rq(1..k);
              end;
	      Random_Product_System.Clear;
    end if;
  end Mixed_Solve;

  procedure Mixed_Solve
                 ( file : in file_type; k,n : in natural; p : in Poly_Sys;
                   tv : in Tree_of_Vectors; bb : out natural;
                   g : in out Poly_Sys; sols : in out Solution_List ) is
  begin
    if k = 0
     then declare
            l : Laur_Sys(g'range);
            tmp : Solution_List;
          begin
	    l := Polynomial_to_Laurent_System(g);
	    Solve(file,l,tv,bb,sols);
	    tmp := sols;
	    while not Is_Null(tmp) loop
	      Head_Of(tmp).t := CMPLX(0.0);
	      tmp := Tail_Of(tmp);
            end loop;
	    Clear(l);
          end;
     elsif k = n
	 then Construct_Random_Product_Start_System.Construct(p,g,sols);
	      bb := Length_Of(sols);
         else Build_RPS(k,n,p);
              declare
                acc : Natural_Vectors.Vector(1..k);
                q : Poly_Sys(1..n-k);
		rq : Poly_Sys(p'range);
		sols_last : Solution_List;
		res : natural;
              begin
                acc := (1..k => 0); res := 0;
                for i in q'range loop
                  q(i) := g(i+k);
                end loop;
                Build_Indices(file,q,1,k,n,acc,tv,res,sols,sols_last);
		bb := res;
	        rq := Random_Product_System.Polynomial_System;
		g(1..k) := rq(1..k);
              end;
	      Random_Product_System.Clear;
    end if;
  end Mixed_Solve;

end Set_Structures_and_Volumes;
