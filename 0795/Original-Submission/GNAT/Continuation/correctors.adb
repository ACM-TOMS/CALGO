with Process_io;                     use Process_io;
with Integer_Vectors;
with Complex_Linear_System_Solvers;  use Complex_Linear_System_Solvers;
with Floating_Equalities;            use Floating_Equalities;

package body Correctors is

-- AUXILIARIES FOR IMPLEMENTING MULTIPLE CORRECTORS :

  procedure Equals ( s : in Solu_Info_Array; x : in Vector; i : in natural;
                     d : in double_float; j : in out natural ) is
    eq : boolean;

  begin
    while j < i loop
      eq := true;
      for k in x'range loop
        eq := Is_Equal(s(j).sol.v(k),x(k),d);
        exit when not eq;
      end loop;
      exit when eq;
      j := j + 1;
    end loop;
  end Equals;

  generic
    with procedure Corrector ( s : in out Solu_Info; c : in Corr_Pars );
  procedure Multiple_Silent_Corrector 
              ( s : in out Solu_Info_Array;
                pivot : in out natural; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean );

  -- DESCRIPTION :
  --   Allows to build any multiple silent corrector,
  --   depending on the corrector supplied as generic parameter.

  generic
    with procedure Corrector ( file : in file_type;
                               s : in out Solu_Info; c : in Corr_Pars );
  procedure Multiple_Reporting_Corrector 
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out natural; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean );

  -- DESCRIPTION :
  --   Allows to build any multiple reporting corrector,
  --   depending on the corrector supplied as generic parameter.

  procedure Multiple_Silent_Corrector 
              ( s : in out Solu_Info_Array;
                pivot : in out natural; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    i,j : natural;
    ffail : boolean := false;

  begin
    i := pivot;
    loop
      Corrector(s(i),c);
      ffail := (((s(i).resa > c.epsaf) and then (s(i).cora > c.epsax))
        or else ((s(i).resr > c.epsrf) and then (s(i).corr > c.epsrx)));
      if not ffail
       then j := 1;
            Equals(s,s(i).sol.v,i,dist_sols,j);
            if j /= i
             then ffail := true;
            end if;
      end if;
      exit when ffail;
      i := i + 1;
      if i > s'last
       then i := s'first;
      end if;
      exit when (i = pivot);
    end loop;
    if ffail
     then pivot := i;
    end if;
    fail := ffail;
  end Multiple_Silent_Corrector;

  procedure Multiple_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out natural; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    i,j : natural;
    ffail : boolean := false;

  begin
    i := pivot;
    loop
      Write_path(file,i);
      Corrector(file,s(i),c);
      sWrite(file,s(i).sol.all);
      ffail := (((s(i).resa > c.epsaf) and then (s(i).cora > c.epsax))
        or else ((s(i).resr > c.epsrf) and then (s(i).corr > c.epsrx)));
      if not ffail
       then j := 1;
            Equals(s,s(i).sol.v,i,dist_sols,j);
            if j /= i
             then ffail := true;
            end if;
      end if;
      exit when ffail;
      i := i + 1;
      if i > s'last
       then i := s'first;
      end if;
      exit when (i = pivot);
    end loop;
    if ffail
     then pivot := i;
    end if;
    fail := ffail;
  end Multiple_Reporting_Corrector;

-- TARGET PROCEDURES :

  procedure Affine_Single_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer;
    y : vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Integer_Vectors.Vector(y'range);
    nit : natural := 0;
    normv : double_float;

  begin
    s.resa := Norm(y);  s.cora := 1.0;               -- INITIALIZATION
    normv := Norm(s.sol.v);
    if normv + 1.0 /= 1.0
     then s.corr := s.cora / normv; s.resr := s.resa / normv;
    end if;
    while nit < c.maxit loop  -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      lufac(j,y'last,ipvt,info);
      exit when info /= 0;              -- STOP WHEN SINGULAR JACOBIAN
      Min_Vector(y);
      lusolve(j,y'last,ipvt,y); s.cora := Norm(y);
      Plus_Vector(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);  s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0
       then s.corr := s.cora / normv; s.resr := s.resa / normv;
      end if;
      nit := nit + 1;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                   -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when numeric_error => return;
  end Affine_Single_Loose_Normal_Silent_Corrector;

  procedure Affine_Single_Loose_Normal_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer;
    y : vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Integer_Vectors.Vector(y'range);
    nit : natural := 0;
    normv : double_float;

  begin
    s.resa := Norm(y);  s.cora := 1.0;               -- INITIALIZATION
    normv := Norm(s.sol.v);
    if normv + 1.0 /= 1.0
     then s.corr := s.cora / normv;
          s.resr := s.resa / normv;
    end if;
    while nit < c.maxit loop  -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      lufac(j,y'last,ipvt,info);
      exit when info /= 0;              -- STOP WHEN SINGULAR JACOBIAN
      Min_Vector(y);
      lusolve(j,y'last,ipvt,y); s.cora := Norm(y);
      Plus_Vector(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);  s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0
       then s.corr := s.cora / normv; s.resr := s.resa / normv;
      end if;
      cWrite(file,s.cora,s.corr,s.resa,s.resr);    -- WRITE STATISTICS 
      nit := nit + 1;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                   -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when numeric_error => return;
  end Affine_Single_Loose_Normal_Reporting_Corrector;

  procedure Affine_Single_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    y : vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Integer_Vectors.Vector(y'range);
    nit : natural := 0;
    normv : double_float;

  begin
    s.resa := Norm(y);  s.cora := 1.0;                 -- INITIALIZATION
    normv := Norm(s.sol.v);
    if normv + 1.0 /= 1.0
     then s.corr := s.cora / normv; s.resr := s.resa / normv;
    end if;
    while nit < c.maxit loop    -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      lufco(j,y'last,ipvt,s.rcond);
      exit when s.rcond + 1.0 = s.rcond;  -- STOP WHEN SINGULAR JACOBIAN
      Min_Vector(y);
      lusolve(j,y'last,ipvt,y); s.cora := Norm(y);
      Plus_Vector(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);  s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0
       then s.corr := s.cora / normv; s.resr := s.resa / normv;
      end if;
      nit := nit + 1;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                    -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when numeric_error => return;
  end Affine_Single_Loose_Conditioned_Silent_Corrector;

  procedure Affine_Single_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    y : vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Integer_Vectors.Vector(y'range);
    nit : natural := 0;
    normv : double_float;

  begin
    s.resa := Norm(y);  s.cora := 1.0;                 -- INITIALIZATION
    normv := Norm(s.sol.v);
    if normv + 1.0 /= 1.0
     then s.corr := s.cora / normv; s.resr := s.resa / normv;
    end if;
    while nit < c.maxit loop    -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      lufco(j,y'last,ipvt,s.rcond);
      exit when s.rcond + 1.0 = s.rcond;  -- STOP WHEN SINGULAR JACOBIAN
      Min_Vector(y);
      lusolve(j,y'last,ipvt,y); s.cora := Norm(y);
      Plus_Vector(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);  s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0
       then s.corr := s.cora / normv; s.resr := s.resa / normv;
      end if;
      cWrite(file,s.cora,s.corr,s.resa,s.resr);      -- WRITE STATISTICS
      cWrite(file,s.rcond,s.sol.m);
      nit := nit + 1;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                     -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when numeric_error => return;
  end Affine_Single_Loose_Conditioned_Reporting_Corrector;

  procedure Affine_Single_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer;
    y : vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Integer_Vectors.Vector(y'range);
    nit : natural := 0;
    normv,ncora,nresa,ncorr,nresr : double_float;
    stop : boolean;

  begin
    s.resa := Norm(y);  s.cora := 1.0;              -- INITIALIZATION
    normv := Norm(s.sol.v);
    if normv + 1.0 /= 1.0
     then s.corr := s.cora / normv; s.resr := s.resa / normv;
    end if;
    while nit < c.maxit loop -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      lufac(j,y'last,ipvt,info);
      exit when info /= 0;             -- STOP WHEN SINGULAR JACOBIAN
      Min_Vector(y);
      lusolve(j,y'last,ipvt,y); ncora := Norm(y);
      Plus_Vector(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);  nresa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0
       then ncorr := ncora / normv; nresr := nresa / normv;
      end if;
      nit := nit + 1;
      stop := (((ncora > s.cora) and then (nresa > s.resa))
               or else ((ncorr > s.corr) and then (ncorr > s.corr)));
                                                 -- STOP WHEN DIVERGENCE
      s.cora := ncora; s.resa := nresa; s.corr := ncorr; s.resr := nresr;
      exit when stop;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                  -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when numeric_error => return;
  end Affine_Single_Severe_Normal_Silent_Corrector;

  procedure Affine_Single_Severe_Normal_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer;
    y : vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Integer_Vectors.Vector(y'range);
    nit : natural := 0;
    normv,ncora,nresa,ncorr,nresr : double_float;
    stop : boolean;

  begin
    s.resa := Norm(y);  s.cora := 1.0;              -- INITIALIZATION
    normv := Norm(s.sol.v);
    if normv + 1.0 /= 1.0
     then s.corr := s.cora / normv; s.resr := s.resa / normv;
    end if;
    while nit < c.maxit loop -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      lufac(j,y'last,ipvt,info);
      exit when info /= 0;             -- STOP WHEN SINGULAR JACOBIAN
      Min_Vector(y);
      lusolve(j,y'last,ipvt,y); ncora := Norm(y);
      Plus_Vector(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);  nresa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0
       then ncorr := ncora / normv; nresr := nresa / normv;
      end if;
      cWrite(file,ncora,ncorr,nresa,nresr);       -- WRITE STATISTICS
      nit := nit + 1;
      stop := (((ncora > s.cora) and then (nresa > s.resa))
              or else ((ncorr > s.corr) and then (ncorr > s.corr)));
                                                 -- STOP WHEN DIVERGENCE
      s.cora := ncora; s.resa := nresa; s.corr := ncorr; s.resr := nresr;
      exit when stop;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                 -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when numeric_error => return;
  end Affine_Single_Severe_Normal_Reporting_Corrector;

  procedure Affine_Single_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    y : vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Integer_Vectors.Vector(y'range);
    nit : natural := 0;
    normv,ncora,nresa,ncorr,nresr : double_float;
    stop : boolean;

  begin
    s.resa := Norm(y);  s.cora := 1.0;                  -- INITIALIZATION
    normv := Norm(s.sol.v);
    if normv + 1.0 /= 1.0
     then s.corr := s.cora / normv; s.resr := s.resa / normv;
    end if;
    while nit < c.maxit loop     -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      lufco(j,y'last,ipvt,s.rcond);
      exit when s.rcond + 1.0 = s.rcond;   -- STOP WHEN SINGULAR JACOBIAN
      Min_Vector(y);
      lusolve(j,y'last,ipvt,y); ncora := Norm(y);
      Plus_Vector(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);  nresa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0
       then ncorr := ncora / normv; nresr := nresa / normv;
      end if;
      nit := nit + 1;
      stop := (((ncora > s.cora) and then (nresa > s.resa))
            or else ((ncorr > s.corr) and then (ncorr > s.corr)));
                                                    -- STOP WHEN DIVERGENCE
      s.cora := ncora; s.resa := nresa; s.corr := ncorr; s.resr := nresr;
      exit when stop;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                      -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when numeric_error => return;
  end Affine_Single_Severe_Conditioned_Silent_Corrector;

  procedure Affine_Single_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    y : vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Integer_Vectors.Vector(y'range);
    nit : natural := 0;
    normv,ncora,nresa,ncorr,nresr : double_float;
    stop : boolean;

  begin
    s.resa := Norm(y);  s.cora := 1.0;                  -- INITIALIZATION
    normv := Norm(s.sol.v);
    if normv + 1.0 /= 1.0
     then s.corr := s.cora / normv; s.resr := s.resa / normv;
    end if;
    while nit < c.maxit loop     -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      lufco(j,y'last,ipvt,s.rcond);
      exit when s.rcond + 1.0 = s.rcond;   -- STOP WHEN SINGULAR JACOBIAN
      Min_Vector(y);
      lusolve(j,y'last,ipvt,y); ncora := Norm(y);
      Plus_Vector(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);  nresa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0
       then ncorr := ncora / normv; nresr := nresa / normv;
      end if;
      cWrite(file,ncora,ncorr,nresa,nresr);           -- WRITE STATISTICS
      cWrite(file,s.rcond,s.sol.m);
      nit := nit + 1;
      stop := (((ncora > s.cora) and then (nresa > s.resa))
            or else ((ncorr > s.corr) and then (ncorr > s.corr)));
                                                     -- STOP WHEN DIVERGENCE
      s.cora := ncora; s.resa := nresa; s.corr := ncorr; s.resr := nresr;
      exit when stop;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                      -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when numeric_error => return;
  end Affine_Single_Severe_Conditioned_Reporting_Corrector;

  procedure Affine_Multiple_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out natural; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Loose_Normal_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Loose_Normal_Silent_Corrector;

  procedure Affine_Multiple_Loose_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out natural; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Loose_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Loose_Normal_Reporting_Corrector;

  procedure Affine_Multiple_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out natural; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Loose_Conditioned_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Loose_Conditioned_Silent_Corrector;

  procedure Affine_Multiple_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out natural; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Loose_Conditioned_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Loose_Conditioned_Reporting_Corrector;

  procedure Affine_Multiple_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out natural; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Severe_Normal_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Severe_Normal_Silent_Corrector;

  procedure Affine_Multiple_Severe_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out natural; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Severe_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Severe_Normal_Reporting_Corrector;

  procedure Affine_Multiple_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out natural; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Severe_Conditioned_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Severe_Conditioned_Silent_Corrector;

  procedure Affine_Multiple_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out natural; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Affine_Single_Severe_Conditioned_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Affine_Multiple_Severe_Conditioned_Reporting_Corrector;

  procedure Projective_Single_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer;
    y : vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Integer_Vectors.Vector(y'range);
    nit : natural := 0;
    normv : double_float;

  begin
    s.resa := Norm(y);  s.cora := 1.0;                -- INITIALIZATION
    normv := Norm(s.sol.v);
    if normv + 1.0 /= 1.0
     then s.corr := s.cora / normv; s.resr := s.resa / normv;
    end if;
    while nit < c.maxit loop   -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        j(j'last(1),jj) := s.sol.v(jj);        -- CORRECT PERPENDICULAR
      end loop;
      lufac(j,y'last,ipvt,info);
      exit when info /= 0;               -- STOP WHEN SINGULAR JACOBIAN
      Min_Vector(y);
      y(y'last) := CMPLX(0.0);               -- IGNORE SCALING EQUATION
      lusolve(j,y'last,ipvt,y); s.cora := Norm(y);
      Plus_Vector(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);
      y(y'last) := CMPLX(0.0);  s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0
       then s.corr := s.cora / normv; s.resr := s.resa / normv;
            for jj in s.sol.v'range loop          -- SCALE THE SOLUTION
              s.sol.v(jj) := s.sol.v(jj)/CMPLX(normv);
            end loop;
      end if;
      nit := nit + 1;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
         and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                    -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when numeric_error => return;
  end Projective_Single_Loose_Normal_Silent_Corrector;

  procedure Projective_Single_Loose_Normal_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer;
    y : vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Integer_Vectors.Vector(y'range);
    nit : natural := 0;
    normv : double_float;

  begin
    s.resa := Norm(y);  s.cora := 1.0;                -- INITIALIZATION
    normv := Norm(s.sol.v);
    if normv + 1.0 /= 1.0
     then s.corr := s.cora / normv; s.resr := s.resa / normv;
    end if;
    while nit < c.maxit loop   -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        j(j'last(1),jj) := s.sol.v(jj);        -- CORRECT PERPENDICULAR
      end loop;
      lufac(j,y'last,ipvt,info);
      exit when info /= 0;               -- STOP WHEN SINGULAR JACOBIAN
      Min_Vector(y);
      y(y'last) := CMPLX(0.0);               -- IGNORE SCALING EQUATION
      lusolve(j,y'last,ipvt,y); s.cora := Norm(y);
      Plus_Vector(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);
      y(y'last) := CMPLX(0.0);  s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0
       then s.corr := s.cora / normv; s.resr := s.resa / normv;
            for jj in s.sol.v'range loop          -- SCALE THE SOLUTION
              s.sol.v(jj) := s.sol.v(jj)/CMPLX(normv);
            end loop;
      end if;
      cWrite(file,s.cora,s.corr,s.resa,s.resr);     -- WRITE STATISTICS
      nit := nit + 1;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                    -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when numeric_error => return;
  end Projective_Single_Loose_Normal_Reporting_Corrector;

  procedure Projective_Single_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    y : vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Integer_Vectors.Vector(y'range);
    nit : natural := 0;
    normv : double_float;

  begin
    s.resa := Norm(y);  s.cora := 1.0;                -- INITIALIZATION
    normv := Norm(s.sol.v);
    if normv + 1.0 /= 1.0
     then s.corr := s.cora / normv; s.resr := s.resa / normv;
    end if;
    while nit < c.maxit loop   -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        j(j'last(1),jj) := s.sol.v(jj);        -- CORRECT PERPENDICULAR
      end loop;
      lufco(j,y'last,ipvt,s.rcond);
      exit when s.rcond + 1.0 = s.rcond; -- STOP WHEN SINGULAR JACOBIAN
      Min_Vector(y);
      y(y'last) := CMPLX(0.0);               -- IGNORE SCALING EQUATION
      lusolve(j,y'last,ipvt,y); s.cora := Norm(y);
      Plus_Vector(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);
      y(y'last) := CMPLX(0.0);  s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0
       then s.corr := s.cora / normv; s.resr := s.resa / normv;
            for jj in s.sol.v'range loop          -- SCALE THE SOLUTION
              s.sol.v(jj) := s.sol.v(jj)/CMPLX(normv);
            end loop;
      end if;
      nit := nit + 1;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                    -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when numeric_error => return;
  end Projective_Single_Loose_Conditioned_Silent_Corrector;

  procedure Projective_Single_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    y : vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Integer_Vectors.Vector(y'range);
    nit : natural := 0;
    normv : double_float;

  begin
    s.resa := Norm(y);  s.cora := 1.0;                -- INITIALIZATION
    normv := Norm(s.sol.v);
    if normv + 1.0 /= 1.0
     then s.corr := s.cora / normv; s.resr := s.resa / normv;
    end if;
    while nit < c.maxit loop   -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        j(j'last(1),jj) := s.sol.v(jj);        -- CORRECT PERPENDICULAR
      end loop;
      lufco(j,y'last,ipvt,s.rcond);
      exit when s.rcond + 1.0 = s.rcond; -- STOP WHEN SINGULAR JACOBIAN
      Min_Vector(y);
      y(y'last) := CMPLX(0.0);               -- IGNORE SCALING EQUATION
      lusolve(j,y'last,ipvt,y); s.cora := Norm(y);
      Plus_Vector(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);
      y(y'last) := CMPLX(0.0);  s.resa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0
       then s.corr := s.cora / normv; s.resr := s.resa / normv;
            for jj in s.sol.v'range loop          -- SCALE THE SOLUTION
              s.sol.v(jj) := s.sol.v(jj)/CMPLX(normv);
            end loop;
      end if;
      cWrite(file,s.cora,s.corr,s.resa,s.resr);     -- WRITE STATISTICS
      cWrite(file,s.rcond,s.sol.m);
      nit := nit + 1;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                    -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when numeric_error => return;
  end Projective_Single_Loose_Conditioned_Reporting_Corrector;

  procedure Projective_Single_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer;
    y : vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Integer_Vectors.Vector(y'range);
    nit : natural := 0;
    normv,ncora,nresa,ncorr,nresr : double_float;
    stop : boolean;

  begin
    s.resa := Norm(y);  s.cora := 1.0;               -- INITIALIZATION
    normv := Norm(s.sol.v);
    if normv + 1.0 /= 1.0
     then s.corr := s.cora / normv; s.resr := s.resa / normv;
    end if;
    while nit < c.maxit loop  -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        j(j'last(1),jj) := s.sol.v(jj);       -- CORRECT PERPENDICULAR
      end loop;
      lufac(j,y'last,ipvt,info);
      exit when info /= 0;              -- STOP WHEN SINGULAR JACOBIAN
      Min_Vector(y);
      y(y'last) := CMPLX(0.0);              -- IGNORE SCALING EQUATION
      lusolve(j,y'last,ipvt,y); ncora := Norm(y);
      Plus_Vector(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);
      y(y'last) := CMPLX(0.0);  nresa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0
       then ncorr := ncora / normv; nresr := nresa / normv;
            for jj in s.sol.v'range loop         -- SCALE THE SOLUTION
              s.sol.v(jj) := s.sol.v(jj)/CMPLX(normv);
            end loop;
      end if;
      nit := nit + 1;
      stop := (((ncora > s.cora) and then (nresa > s.resa))
             or else ((ncorr > s.corr) and then (ncorr > s.corr)));
                                                    -- STOP WHEN DIVERGENCE
      s.cora := ncora; s.resa := nresa; s.corr := ncorr; s.resr := nresr;
      exit when stop;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                   -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when numeric_error => return;
  end Projective_Single_Severe_Normal_Silent_Corrector;

  procedure Projective_Single_Severe_Normal_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    info : integer;
    y : vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Integer_Vectors.Vector(y'range);
    nit : natural := 0;
    normv,ncora,nresa,ncorr,nresr : double_float;
    stop : boolean;

  begin
    s.resa := Norm(y);  s.cora := 1.0;               -- INITIALIZATION 
    normv := Norm(s.sol.v);
    if normv + 1.0 /= 1.0
     then s.corr := s.cora / normv; s.resr := s.resa / normv;
    end if;
    while nit < c.maxit loop  -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        j(j'last(1),jj) := s.sol.v(jj);       -- CORRECT PERPENDICULAR
      end loop;
      lufac(j,y'last,ipvt,info);
      exit when info /= 0;              -- STOP WHEN SINGULAR JACOBIAN
      Min_Vector(y);
      y(y'last) := CMPLX(0.0);              -- IGNORE SCALING EQUATION
      lusolve(j,y'last,ipvt,y); ncora := Norm(y);
      Plus_Vector(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);
      y(y'last) := CMPLX(0.0);  nresa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0
       then ncorr := ncora / normv; nresr := nresa / normv;
            for jj in s.sol.v'range loop         -- SCALE THE SOLUTION
              s.sol.v(jj) := s.sol.v(jj)/CMPLX(normv);
            end loop;
      end if;
      cWrite(file,ncora,ncorr,nresa,nresr);        -- WRITE STATISTICS
      nit := nit + 1;
      stop := (((ncora > s.cora) and then (nresa > s.resa))
            or else ((ncorr > s.corr) and then (ncorr > s.corr)));
                                                   -- STOP WHEN DIVERGENCE
      s.cora := ncora; s.resa := nresa; s.corr := ncorr; s.resr := nresr;
      exit when stop;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
         and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                   -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when numeric_error => return;
  end Projective_Single_Severe_Normal_Reporting_Corrector;

  procedure Projective_Single_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info; c : in Corr_Pars ) is

    y : vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Integer_Vectors.Vector(y'range);
    nit : natural := 0;
    normv,ncora,nresa,ncorr,nresr : double_float;
    stop : boolean;

  begin
    s.resa := Norm(y);  s.cora := 1.0;               -- INITIALIZATION
    normv := Norm(s.sol.v);
    if normv + 1.0 /= 1.0
     then s.corr := s.cora / normv; s.resr := s.resa / normv;
    end if;
    while nit < c.maxit loop    -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        j(j'last(1),jj) := s.sol.v(jj);         -- CORRECT PERPENDICULAR
      end loop;
      lufco(j,y'last,ipvt,s.rcond);
      exit when s.rcond + 1.0 = s.rcond;  -- STOP WHEN SINGULAR JACOBIAN
      Min_Vector(y);
      y(y'last) := CMPLX(0.0);                -- IGNORE SCALING EQUATION
      lusolve(j,y'last,ipvt,y); ncora := Norm(y);
      Plus_Vector(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);
      y(y'last) := CMPLX(0.0);  nresa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0
       then ncorr := ncora / normv; nresr := nresa / normv;
            for jj in s.sol.v'range loop           -- SCALE THE SOLUTION
              s.sol.v(jj) := s.sol.v(jj)/CMPLX(normv);
            end loop;
      end if;
      nit := nit + 1;
      stop := (((ncora > s.cora) and then (nresa > s.resa))
             or else ((ncorr > s.corr) and then (ncorr > s.corr)));
                                                    -- STOP WHEN DIVERGENCE
      s.cora := ncora; s.resa := nresa; s.corr := ncorr; s.resr := nresr;
      exit when stop;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                     -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when numeric_error => return;
  end Projective_Single_Severe_Conditioned_Silent_Corrector;

  procedure Projective_Single_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type;
                s : in out Solu_Info; c : in Corr_Pars ) is

    y : vector(s.sol.v'range) := H(s.sol.v,s.sol.t);
    j : Matrix(y'range,y'range);
    ipvt : Integer_Vectors.Vector(y'range);
    nit : natural := 0;
    normv,ncora,nresa,ncorr,nresr : double_float;
    stop : boolean;

  begin
    s.resa := Norm(y);  s.cora := 1.0;               -- INITIALIZATION
    normv := Norm(s.sol.v);
    if normv + 1.0 /= 1.0
     then s.corr := s.cora / normv; s.resr := s.resa / normv;
    end if;
    while nit < c.maxit loop    -- STOP WHEN MAXIMUM #ITERATIONS REACHED
      j := dH(s.sol.v,s.sol.t);
      for jj in y'range loop
        j(j'last(1),jj) := s.sol.v(jj);         -- CORRECT PERPENDICULAR
      end loop;
      lufco(j,y'last,ipvt,s.rcond);
      exit when s.rcond + 1.0 = s.rcond;  -- STOP WHEN SINGULAR JACOBIAN
      Min_Vector(y);
      y(y'last) := CMPLX(0.0);                -- IGNORE SCALING EQUATION
      lusolve(j,y'last,ipvt,y); ncora := Norm(y);
      Plus_Vector(s.sol.v,y);   s.length_path := s.length_path + s.cora;
      y := H(s.sol.v,s.sol.t);
      y(y'last) := CMPLX(0.0);  nresa := Norm(y);
      normv := Norm(s.sol.v);
      if normv + 1.0 /= 1.0
       then ncorr := ncora / normv; nresr := nresa / normv;
            for jj in s.sol.v'range loop           -- SCALE THE SOLUTION
              s.sol.v(jj) := s.sol.v(jj)/CMPLX(normv);
            end loop;
      end if;
      cWrite(file,ncora,ncorr,nresa,nresr);          -- WRITE STATISTICS
      cWrite(file,s.rcond,s.sol.m);
      nit := nit + 1;
      stop := (((ncora > s.cora) and then (nresa > s.resa))
              or else ((ncorr > s.corr) and then (ncorr > s.corr)));
                                                   -- STOP WHEN DIVERGENCE
      s.cora := ncora; s.resa := nresa; s.corr := ncorr; s.resr := nresr;
      exit when stop;
      exit when (((s.cora <= c.epsax) or else (s.resa <= c.epsaf))
        and then ((s.corr <= c.epsrx) or else (s.resr <= c.epsrf)));
    end loop;                     -- STOP WHEN DESIRED PRECISION REACHED
    s.niter := s.niter + nit; s.nsyst := s.nsyst + nit;
  exception
    when numeric_error => return;
  end Projective_Single_Severe_Conditioned_Reporting_Corrector;

  procedure Projective_Multiple_Loose_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out natural; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Loose_Normal_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Loose_Normal_Silent_Corrector;

  procedure Projective_Multiple_Loose_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out natural; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Loose_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Loose_Normal_Reporting_Corrector;

  procedure Projective_Multiple_Loose_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out natural; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Loose_Conditioned_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Loose_Conditioned_Silent_Corrector;

  procedure Projective_Multiple_Loose_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out natural; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Loose_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Loose_Conditioned_Reporting_Corrector;

  procedure Projective_Multiple_Severe_Normal_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out natural; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Severe_Normal_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Severe_Normal_Silent_Corrector;

  procedure Projective_Multiple_Severe_Normal_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out natural; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Loose_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Severe_Normal_Reporting_Corrector;

  procedure Projective_Multiple_Severe_Conditioned_Silent_Corrector
              ( s : in out Solu_Info_Array;
                pivot : in out natural; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Severe_Conditioned_Silent_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Silent_Corrector(Single_Corrector);

  begin
    Corrector(s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Severe_Conditioned_Silent_Corrector;

  procedure Projective_Multiple_Severe_Conditioned_Reporting_Corrector
              ( file : in file_type; s : in out Solu_Info_Array;
                pivot : in out natural; dist_sols : in double_float;
                c : in Corr_Pars; fail : out boolean ) is

    procedure Single_Corrector is
      new Projective_Single_Severe_Conditioned_Reporting_Corrector(Norm,H,dH);
    procedure Corrector is new Multiple_Reporting_Corrector(Single_Corrector);

  begin
    Corrector(file,s,pivot,dist_sols,c,fail);
  end Projective_Multiple_Severe_Conditioned_Reporting_Corrector;

end Correctors;
