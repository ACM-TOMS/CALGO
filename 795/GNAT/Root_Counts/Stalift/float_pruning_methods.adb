with Floating_Point_Numbers;              use Floating_Point_Numbers;
with Float_Matrices;                      use Float_Matrices;
with Float_Linear_System_Solvers;         use Float_Linear_System_Solvers;
with Farkas_Lemma;                        use Farkas_Lemma;
with Lists_of_Float_Vectors;              use Lists_of_Float_Vectors;

package body Float_Pruning_Methods is

-- GENERAL AUXILIARIES :

  procedure Normalize ( v : in out Float_Vectors.Vector ) is

  -- DESCRIPTION : Divides every entry by v(v'last).

  begin
    for i in v'range loop
      v(i) := v(i)/v(v'last);
    end loop;
  end Normalize;

  function Convert ( fa : Face_Array ) return Array_of_Lists is

  -- DESCRIPTION :
  --   Converts the array of faces into an array of lists by
  --   converting the first element of each list of faces.

    res : Array_of_Lists(fa'range);

  begin
    for k in fa'range loop
      res(k) := Shallow_Create(fa(k).all);
    end loop;
    return res;
  end Convert;

-- AUXILIARIES FOR THE PRUNING ALGORITHMS :

  procedure Initialize ( n : in natural; mat : out Matrix; 
                         ipvt : out Vector ) is

  -- DESCRIPTION :
  --   Initializes an n*(n+1) matrix with zeroes and the pivoting vector
  --   with the entries 1..n.

  begin
    for i in 1..n loop
      for j in 1..n+1 loop
        mat(i,j) := 0.0;
      end loop;
    end loop;
    for i in 1..n loop
      ipvt(i) := i;
    end loop;
  end Initialize;

  function Number_of_Inequalities
              ( mix : Vector; lifted : Array_of_Lists ) return natural is

  -- DESCRIPTION :
  --   Returns the maximal number of inequalities for pruning.

    res : natural := 0;

  begin
    for k in lifted'range loop
      res := res + Length_Of(lifted(k)) - mix(k) - 1;
    end loop;
    return res;
  end Number_of_Inequalities;

  procedure Ordered_Inequalities ( k : in natural; mat : in out Matrix ) is

  -- DESCRIPTION :
  --   Defines k inequalities mat(k,k) - mat(k+1,k) >= 0.

  begin
    for i in mat'first(1)..k loop
      for j in mat'range(2) loop
        mat(i,j) := 0.0;
      end loop;
      mat(i,i) := 1.0; mat(i,i+1) := -1.0;
    end loop;
  end Ordered_Inequalities;

  procedure Check_and_Update
                ( mic : in Face_Array; lifted : in Array_of_Lists;
                  m : in Matrix; ipvt : in Vector; tol : in double_float;
                  mixsub,mixsub_last : in out Mixed_Subdivision ) is

  -- DESCRIPTION :
  --   Computes the normal to the points in pts, by solving the 
  --   linear system defined by m and ipvt.
  --   If the computed normal is an inner normal w.r.t. the lifted points,
  --   then the mixed subdivision will be updated with a new cell.

    use Float_Vectors;
    v : Float_Vectors.Vector(m'range(2)) := Solve(m,tol,ipvt);
    pts : Array_of_Lists(mic'range);

  begin
    if abs(v(v'last)) > tol
     then Normalize(v);
          if v(v'last) < 0.0
           then Min_Vector(v);
          end if;
          pts := Convert(mic);
          Update(pts,v,mixsub,mixsub_last);
    end if;
  end Check_and_Update;

  procedure Create_Equalities
                ( n : in natural; f : in Face; mat,ineq : in Matrix;
                  tol : in double_float;
                  ipvt : in Vector; row,rowineq : in natural;
                  newmat,newineq : in out Matrix; newipvt : in out Vector;
                  newrow,newrowineq : in out natural; fail : out boolean ) is

  -- DESCRIPTION :
  --   Creates new equalities and uses them to eliminate unknowns in the
  --   matrices of equalities and inequalities.  Failure is reported when
  --   a zero row is encountered.  On entry, all new* variables must be
  --   initialized with the corresponding *-ones.

    shi : Float_Vectors.Vector(1..n+1) := f(f'first).all;
    fl : boolean := false;
    pivot : natural;

  begin
    for i in f'first+1..f'last loop
      newrow := newrow + 1;
      for j in f(i)'range loop
        newmat(newrow,j) := f(i)(j) - shi(j);
      end loop;
      Switch(newipvt,newrow,newmat);
      Upper_Triangulate(newrow,newmat,tol,newipvt,pivot);
      fl := (pivot = 0);
      exit when fl;
      Switch(newrow,pivot,ineq'first,rowineq,newineq);
    end loop;
    fail := fl;
  end Create_Equalities;

  procedure Complementary_Slackness
                  ( tableau : in out matrix;
                    tol : in double_float; feasible : out boolean ) is

    lastcol : constant integer := tableau'last(2)-1;
    rhs,sol : Float_Vectors.Vector(tableau'range(1));
    columns : Vector(sol'range);

  begin
    for i in rhs'range loop
      rhs(i) := double_float(tableau(i,tableau'last(2)));
    end loop;
    Complementary_Slackness(tableau,lastcol,rhs,tol,sol,columns,feasible);
  end Complementary_Slackness;

  function Check_Feasibility ( k,m,n : natural; ineq : Matrix;
                               tol : double_float ) return boolean is
 
  -- DESCRIPTION :
  --   Returns true if -v(n+1) > 0 can be derived, false otherwise.

  -- ON ENTRY :
  --   k      current unknown that has been eliminated,
  --           for all i <= k : ineq(l,i) = 0, for l in ineq'first..m;
  --   m      number of inequalities;
  --   n      dimension;
  --   ineq   matrix of inequalities.

    tableau : Matrix(1..n-k+1,1..m+1);
    feasi : boolean;

  begin
    if m = 0
     then feasi := false;
     else for i in k+1..n+1 loop
            for j in 1..m loop
              tableau(i-k,j) := ineq(j,i);
            end loop;
            tableau(i-k,m+1) := 0.0;
          end loop;
          tableau(n-k+1,m+1) := -1.0;
          Complementary_Slackness(tableau,tol,feasi);
    end if;
    return feasi;
  end Check_Feasibility;

  procedure Update_Inequalities
               ( k,rowmat1,rowmat2,n : in natural;
                 mat : in Matrix; ipvt : in Vector; tol : in double_float;
                 rowineq : in out natural; ineq : in out Matrix;
                 lifted : in Array_of_Lists; mic : in out Face_Array ) is

  -- DESCRIPTION :
  --   The inequalities will be updated w.r.t. the equality
  --   constraints on the inner normal.

    tmp : List;
    pt,shi : Float_Vectors.Link_to_Vector;

  begin
    for i in ineq'first..rowineq loop      -- update the old inequalities
      for j in rowmat1..rowmat2 loop
        Upper_Triangulate(j,mat,tol,i,ineq);
      end loop;
    end loop;
    shi := mic(k)(mic(k)'first);
    tmp := lifted(k);                      -- make new inequalities
    while not Is_Null(tmp) loop
      pt := Head_Of(tmp);
      if not Is_In(mic(k),pt.all)
       then rowineq := rowineq + 1;
            for j in pt'range loop
              ineq(rowineq,j) := pt(j) - shi(j);
            end loop;
            Switch(ipvt,rowineq,ineq);
            for i in 1..rowmat2 loop
              Upper_Triangulate(i,mat,tol,rowineq,ineq);
            end loop;
      end if;
      tmp := Tail_Of(tmp);
    end loop;
  end Update_Inequalities;

-- CONSTRUCTION WITH PRUNING :

  procedure Gen1_Create
               ( n : in natural; mix : in Vector; fa : in Array_of_Faces; 
                 lifted : in Array_of_Lists; tol : in double_float;
                 nbsucc,nbfail : in out Float_Vectors.Vector;
                 mixsub : out Mixed_Subdivision ) is

    res,res_last : Mixed_Subdivision;
    accu : Face_Array(fa'range);
    ma : Matrix(1..n,1..n+1);
    ipvt : Vector(1..n);
    ineqrows : natural;

    procedure Compute_Mixed_Cells
                 ( k,row : in natural; mat : in Matrix; ipvt : in Vector;
                   rowineq : in natural; ineq : in Matrix;
                   mic : in out Face_Array; continue : out boolean );

    -- DESCRIPTION :
    --   Backtrack recursive procedure to enumerate face-face combinations.

    -- ON ENTRY :
    --   k         index for current point set;
    --   row       number of rows already in the matrix;
    --   mat       matrix which determines the inner normal;
    --   ipvt      contains the pivoting information;
    --   rowineq   number of inequality constraints already in ineq;
    --   ineq      matrix for the inequality constraints on the
    --             inner normal v, of type <.,v> >= 0;
    --   mic       contains the current selected faces, up to k-1.

    -- ON RETURN :
    --   mic       updated selected faces.
    --   continue  indicates whether to continue the creation or not.

    procedure Process_Inequalities
                 ( k,rowmat1,rowmat2 : in natural;
                   mat : in Matrix; ipvt : in Vector;
                   rowineq : in out natural; ineq : in out Matrix;
                   mic : in out Face_Array; cont : out boolean ) is

    -- DESCRIPTION :
    --   Updates inequalities and checks feasibility before proceeding.

      fl : boolean := false;
      tmp : List;
      pt,shi : Link_to_Vector;

    begin
      Update_Inequalities(k,rowmat1,rowmat2,n,mat,ipvt,tol,rowineq,ineq,
                          lifted,mic);
      if Check_Feasibility(rowmat2,rowineq,n,ineq,tol)
       then nbfail(k) := nbfail(k) + 1.0;
            cont := true;
       else nbsucc(k) := nbsucc(k) + 1.0;
            Compute_Mixed_Cells(k+1,rowmat2,mat,ipvt,rowineq,ineq,mic,cont);
      end if;
    end Process_Inequalities;

    procedure Compute_Mixed_Cells
                 ( k,row : in natural; mat : in Matrix; ipvt : in Vector;
                   rowineq : in natural; ineq : in Matrix;
                   mic : in out Face_Array; continue : out boolean ) is

      old : Mixed_Subdivision := res_last;
      cont : boolean := true;
      tmpfa : Faces;

    begin
      if k > mic'last
       then Check_and_Update(mic,lifted,mat,ipvt,tol,res,res_last);
            if old /= res_last
             then Process(Head_Of(res_last),continue);
             else continue := true;
            end if;
       else tmpfa := fa(k);
            while not Is_Null(tmpfa) loop  -- enumerate faces of kth polytope
              mic(k) := Head_Of(tmpfa);
              declare                                 -- update the matrices
                fl : boolean;
                newipvt : Vector(ipvt'range) := ipvt;
                newmat : Matrix(mat'range(1),mat'range(2)) := mat;
                newineq : Matrix(ineq'range(1),ineq'range(2)) := ineq;
                newrow : natural := row;
                newrowineq : natural := rowineq;
              begin
                Create_Equalities
                   (n,mic(k),mat,ineq,tol,ipvt,row,rowineq,newmat,newineq,
                    newipvt,newrow,newrowineq,fl);
                if fl
                 then nbfail(k) := nbfail(k) + 1.0;
                 else Process_Inequalities(k,row+1,newrow,newmat,newipvt,
                                           newrowineq,newineq,mic,cont);
                end if;
              end;
              tmpfa := Tail_Of(tmpfa);
              exit when not cont;
            end loop;
            continue := cont;
      end if;
    end Compute_Mixed_Cells;

  begin
    Initialize(n,ma,ipvt);
    ineqrows := Number_of_Inequalities(mix,lifted);
    declare
      ineq : matrix(1..ineqrows,1..n+1);
      cont : boolean;
    begin
      ineq(1,1) := 0.0;
      Compute_Mixed_Cells(accu'first,0,ma,ipvt,0,ineq,accu,cont);
    end;
    mixsub := res;
  end Gen1_Create;

  procedure Create
              ( n : in natural; mix : in Vector; fa : in Array_of_Faces; 
                lifted : in Array_of_Lists; tol : in double_float;
                nbsucc,nbfail : in out Float_Vectors.Vector;
		mixsub : out Mixed_Subdivision ) is

    res,res_last : Mixed_Subdivision;
    accu : Face_Array(fa'range);
    ma : Matrix(1..n,1..n+1);
    ipvt : Vector(1..n);
    ineqrows : natural;

    procedure Compute_Mixed_Cells
                 ( k,row : in natural; mat : in Matrix; ipvt : in Vector;
                   rowineq : in natural; ineq : in Matrix;
                   mic : in out Face_Array );

    -- DESCRIPTION :
    --   Backtrack recursive procedure to enumerate face-face combinations.

    -- ON ENTRY :
    --   k         index for current point set;
    --   row       number of rows already in the matrix;
    --   mat       matrix which determines the inner normal;
    --   ipvt      contains the pivoting information;
    --   rowineq   number of inequality constraints already in ineq;
    --   ineq      matrix for the inequality constraints on the
    --             inner normal v, of type <.,v> >= 0;
    --   mic       contains the current selected faces, up to k-1.

    -- ON RETURN :
    --   mic       updated selected faces.

    procedure Process_Inequalities
                 ( k,rowmat1,rowmat2 : in natural;
                   mat : in matrix; ipvt : in vector;
                   rowineq : in out natural; ineq : in out matrix;
                   mic : in out Face_Array) is

    -- DESCRIPTION :
    --   Updates inequalities and checks feasibility before proceeding.

      tmp : List;
      pt,shi : Link_to_Vector;

    begin
      Update_Inequalities(k,rowmat1,rowmat2,n,mat,ipvt,tol,rowineq,ineq,
                          lifted,mic);
      if Check_Feasibility(rowmat2,rowineq,n,ineq,tol)
       then nbfail(k) := nbfail(k) + 1.0;
       else nbsucc(k) := nbsucc(k) + 1.0;
            Compute_Mixed_Cells(k+1,rowmat2,mat,ipvt,rowineq,ineq,mic);
      end if;
    end Process_Inequalities;

    procedure Compute_Mixed_Cells 
                 ( k,row : in natural; mat : in matrix;
                   ipvt : in vector; rowineq : in natural; ineq : in matrix;
                   mic : in out Face_Array ) is

      tmpfa : Faces;

    begin
      if k > mic'last
       then Check_and_Update(mic,lifted,mat,ipvt,tol,res,res_last);
       else tmpfa := fa(k);
            while not Is_Null(tmpfa) loop  -- enumerate faces of kth polytope
              mic(k) := Head_Of(tmpfa);
              declare                                      -- update matrices
                fl : boolean;
                newipvt : Vector(ipvt'range) := ipvt;
                newmat : Matrix(mat'range(1),mat'range(2)) := mat;
                newineq : Matrix(ineq'range(1),ineq'range(2)) := ineq;
                newrow : natural := row;
                newrowineq : natural := rowineq;
              begin
                Create_Equalities
                    (n,mic(k),mat,ineq,tol,ipvt,row,rowineq,newmat,newineq,
                     newipvt,newrow,newrowineq,fl);
                if fl 
                 then nbfail(k) := nbfail(k) + 1.0;
                 else Process_Inequalities
                        (k,row+1,newrow,newmat,newipvt,newrowineq,newineq,mic);
                end if;
              end;
              tmpfa := Tail_Of(tmpfa);
            end loop;
      end if;
    end Compute_Mixed_Cells;

  begin
    Initialize(n,ma,ipvt);
    ineqrows := Number_of_Inequalities(mix,lifted);
    declare
      ineq : Matrix(1..ineqrows,1..n+1);
    begin
      ineq(1,1) := 0.0;
      Compute_Mixed_Cells(accu'first,0,ma,ipvt,0,ineq,accu);
    end;
    mixsub := res;
  end Create;

end Float_Pruning_Methods;
