with integer_io;                       use integer_io;
with Float_Vectors_io;                 use Float_Vectors_io;
with Complex_Numbers_io,Process_io;    use Complex_Numbers_io,Process_io;

with Mathematical_Functions;           use Mathematical_Functions;
with Floating_Equalities;              use Floating_Equalities;
with Predictors,Correctors;            use Predictors,Correctors;
with Dispatch_Predictors;              use Dispatch_Predictors;
with Continuation_Parameters;          use Continuation_Parameters;
with Directions_of_Solution_Paths;     use Directions_of_Solution_Paths;
with Float_Vectors_of_Vectors;

package body Path_Trackers is

  use Floating_Point_Numbers.double_float_io;

-- AUXILIAIRIES :

  function At_End ( t,target : double_complex; distance,tol : double_float )
                  return boolean is

  -- DESCRIPTION :
  --   Decides whether at end of continuation stage.

  begin
    if distance < tol
     then return Is_Equal(t,target,tol);
     else if modulus(t-target) <= distance
           then return true;
           else return false;
          end if;
    end if;
  end At_End;

  function Stop ( p : Pred_Pars; t,target : double_complex;
                  step : double_float ) return boolean is

  -- DESCRIPTION :
  --   Returns true if either the step size is smaller than p.minstep, or
  --   or alternatively, in case of geometric predictor, if the distance to
  --   the target has become smaller than p.minstep.

  begin
    if step <= p.minstep
     then return true;
     else if (p.predictor_type = 2) or (p.predictor_type = 5)
           then return (modulus(t-target) <= p.minstep);
           else return false;
          end if;
    end if;
  end Stop;

  function Is_Multiple ( a,b,tol : double_float ) return integer is

  -- DESCRIPTION :
  --   Returns a/b if a is a multiple of b, returns 0 in the other case.

    fq : double_float;
    iq : integer;
  begin
    if abs(b) < tol
     then return 0;
     else fq := a/b;
          iq := integer(fq);
          if abs(fq - double_float(iq)) <= tol
           then return iq;
           else return 0;
          end if;
    end if;
  end Is_Multiple;

  procedure Linear_Step_Control 
              ( success : in boolean; p : in Pred_Pars;
                step : in out double_float;
                nsuccess : in out natural; trial : in natural ) is

  -- DESCRIPTION :
  --   Control of step size for linear path following.
  --   With geometric prediction, the ratio (=step) will be enlarged
  --   when not success.  In order not to break the sequence, the ratio
  --   is not reduced when success.

  begin
    if (p.predictor_type = 2) or (p.predictor_type = 5)
     then if success
           then nsuccess := nsuccess + 1;
           else nsuccess := 0;
                if p.expfac < 1.0
                 then step := step + p.expfac*(1.0 - step);
                 else step := step + (1.0 - step)/p.expfac;
                end if;
                if step >= 1.0
                 then step := p.minstep/2.0;  -- that ends the game
                end if;
          end if;
     else if success
           then nsuccess := nsuccess + 1;
                if nsuccess > p.success_steps
                 then step := step*p.expfac;
                      if step > p.maxstep
                       then step := p.maxstep;
                      end if;
                end if;
           else nsuccess := 0;
                if trial mod 3 = 0
                 then step := step*p.redfac;
                end if;
          end if;
    end if;
  end Linear_Step_Control;

  procedure Circular_Step_Control
             ( success : in boolean; p : in Pred_Pars; twopi : in double_float;
               step : in out double_float; nsuccess : in out natural ) is

  -- DESCRIPTION :
  --   Control of step size for circular path following, note that the
  --   step size should be some multiple of pi.

    maxstep : constant double_float := p.maxstep*twopi;

  begin
    if success
     then nsuccess := nsuccess + 1;
          if nsuccess > p.success_steps
           then step := step*2.0;          -- expansion factor = 2
                if step > maxstep
                 then step := maxstep;
                end if;
          end if;
     else nsuccess := 0;
          step := step*0.5;                -- reduction factor = 1/2
    end if;
  end Circular_Step_Control;

  procedure Set_Corrector_Parameters
              ( c : in out Corr_Pars; eps : double_float ) is

  -- DESCRIPTION :
  --   All eps* parameters in c are set to eps.

  begin
    c.epsrx := eps; c.epsax := eps; c.epsrf := eps; c.epsaf := eps;
  end Set_Corrector_Parameters;

  function End_Game_Corrector_Parameters
             ( current : Corr_Pars; distance,tol : double_float ) 
             return Corr_Pars is

  -- DESCRIPTION :
  --   Returns corrector parameters for the end game of the first or the 
  --   second continuation stage, depending on the distance from the target.

    res : Corr_Pars := current;

  begin
    if distance < tol                 -- correct further to detect clustering
     then res := current;
          Set_Corrector_Parameters(res,tol);
     else res := Create_End_Game;     -- or to move to end game more smoothly
    end if;
    return res;
  end End_Game_Corrector_Parameters;

-- MANAGEMENT OF DATA DURING PATH FOLLOWING :

  procedure Linear_Single_Initialize 
                      ( s : in Solu_Info; p : in Pred_Pars;
                        old_t,prev_t : out double_complex;
                        old_solution,prev_solution : out Vector ) is

  -- DESCRIPTION :
  --   Initialization for linear path following of one path.

  -- ON ENTRY :
  --   s                solution at beginning of path;
  --   p                predictor parameters.

  -- ON RETURN :
  --   old_t            back up value for continuation parameter;
  --   prev_t           previous value of continuation parameter;
  --   old_solution     back up value for solution vector;
  --   prev_solution    previous value of solution vector;

  begin
    old_t := s.sol.t; old_solution := s.sol.v;
    if p.predictor_type <= 2                    -- for all secant predictors
     then prev_t := s.sol.t;
          prev_solution := s.sol.v;
    end if;
  end Linear_Single_Initialize;

  procedure Linear_Single_Management
                 ( s : in out Solu_Info; p : in Pred_Pars; c : in Corr_Pars;
                   old_t,prev_t : in out double_complex;
                   old_solution,prev_solution,old_v,prev_v,vv : in out Vector;
                   step : in out double_float; nsuccess,trial : in out natural;
                   success : in out boolean ) is

  -- DESCRIPTION :
  --   Management of data after correction during linear path following.

  -- PARAMETERS :
  --   s                current solution;
  --   p                predictor parameters;
  --   c                corrector parameters;
  --   old_t            back up value for continuation parameter;
  --   prev_t           previous value of continuation parameter;
  --   old_solution     back up value for solution vector;
  --   prev_solution    previous value of solution vector;
  --   old_v            back up value for tangent direction;
  --   prev_v           previous value for tangent direction;
  --   vv               current tangent direction;
  --   step             current step size;
  --   nsuccess         number of consecutive successes;
  --   trial            number of trials after failure;
  --   success          successful correction step.

  begin
    success := (((s.resa <= c.epsaf) or else (s.cora <= c.epsax))
        or else ((s.resr <= c.epsrf) or else (s.corr <= c.epsrx)));
    if ((p.predictor_type <= 2) and then success)         -- secant predictors
     then prev_t := old_t;
          prev_solution := old_solution;
    end if;
    if ((p.predictor_type = 6) and then success)          -- Hermite predictor
     then prev_t := old_t;
          prev_solution := old_solution;
          prev_v := old_v;
    end if;
    if ((p.predictor_type = 1) or (p.predictor_type = 3)) -- complex predictors
     then if success
           then trial := 0;
           else trial := trial + 1;
          end if;
    end if;
    s.nstep := s.nstep + 1;
    if not success
     then s.nfail := s.nfail + 1;
     end if;
    Linear_Step_Control(success,p,step,nsuccess,trial);
    if step < p.minstep
     then return;
    end if;
    if not success
     then s.sol.t := old_t;
          s.sol.v := old_solution;
     else old_t := s.sol.t;
          old_solution := s.sol.v;
          if p.predictor_type = 6
           then old_v := vv;
          end if;
    end if;
  end Linear_Single_Management;

  procedure Linear_Multiple_Initialize 
                 ( s : in Solu_Info_Array; p : in Pred_Pars;
                   t,old_t,prev_t : out double_complex;
                   sa,old_sa,prev_sa : in out Solution_Array ) is

  -- DECRIPTION :
  --   Initialization for linear path following for more than one path.

  begin
    t := s(s'first).sol.t;
    old_t := s(s'first).sol.t;
    Copy(s,sa); Copy(sa,old_sa);
    case p.predictor_type is
      when 0 | 1 | 2 => prev_t := s(s'first).sol.t; Copy(sa,prev_sa);
      when others => null;
    end case;
  end Linear_Multiple_Initialize;

  procedure Linear_Multiple_Management
                  ( s : in out Solu_Info_array;
                    sa,old_sa,prev_sa : in out Solution_Array;
                    t,old_t,prev_t : in out double_complex; p : in Pred_Pars; 
                    step : in out double_float; pivot : in natural; 
                    nsuccess,trial : in out natural; success : in boolean ) is

  -- DESCRIPTION :
  --   Management of data after correction during linear path following.

  -- PARAMETERS :
  --   s            current solutions with information statistics;
  --   sa           current solutions;
  --   old_sa       back up value for solutions;
  --   prev_sa      previous solutions;
  --   t            current value of continuation parameter;
  --   old_t        back up value for continuation parameter;
  --   prev_t       previous value of continuation parameter;
  --   p            predictor parameters;
  --   step         current step size;
  --   pivot        solution where correction is started;
  --   nsuccess     number of consecutive successes;
  --   trial        number of trials after failure;
  --   success      successful correction step.

  begin
    if ((p.predictor_type <= 2) and then success)      -- for secant predictors
     then prev_t := old_t; Copy(old_sa,prev_sa);
    end if;
    if ((p.predictor_type = 1) or (p.predictor_type = 3)) -- complex predictors
     then if success
           then trial := 0;
           else trial := trial + 1;
          end if;
    end if;
    for k in s'range loop
      s(k).nstep := s(k).nstep + 1;
    end loop;
    if not success
     then s(pivot).nfail := s(pivot).nfail + 1;
    end if;
    Linear_Step_Control(success,p,step,nsuccess,trial);
    if step < p.minstep then return; end if;
    if success
     then Copy(sa,old_sa); old_t := t;
     else Copy(old_sa,sa); t := old_t;
    end if;
  end Linear_Multiple_Management;

  procedure Circular_Management
                   ( s : in out Solu_Info; p : in Pred_Pars; c : in Corr_Pars;
                     old_t,prev_t : in out double_complex;
                     start_solution : in Vector;
                     old_solution,prev_solution,w_sum,w_all_sum : in out Vector;
                     twopi,epslop,tol : in double_float;
                     theta,old_theta,step : in out double_float;
                     nsuccess,n_sum,n_all_sum,w_c : in out natural;
                     max_wc : in natural; stop,success : in out boolean ) is

  -- DESCRIPTION :
  --   Management of circular path following.

  -- PARAMETERS :
  --   s                current solution;
  --   p                predictor parameters;
  --   c                corrector parameters;
  --   old_t            back up value for continuation parameter;
  --   prev_t           previous value of continuation parameter;
  --   start_solution   solution vector at start of continuation;
  --   old_solution     back up value for solution vector;
  --   prev_solution    previous value of solution vector;
  --   w_sum            sum of cycle;
  --   w_all_sum        sum of all cycles;
  --   twopi            two times PI;
  --   epslop           tolerance to decide whether two vectors are equal;
  --   theta            current value of theta;
  --   old_theta        back up value for theta;
  --   step             current step size;
  --   nsuccess         number of consecutive successes;
  --   n_sum            number of cycles;
  --   n_all_sum        total number of cycles;
  --   w_c              winding number;
  --   max_wc           upper bound on winding number;
  --   stop             true when winding number has been computed;
  --   success          successful correction step.

    tmp : integer;

  begin
    success := (((s.resa <= c.epsaf) or else (s.cora <= c.epsax))
        or else ((s.resr <= c.epsrf) or else (s.corr <= c.epsrx)));
    if p.predictor_type = 0 and then success
     then prev_t := old_t;
          prev_solution := old_solution;
    end if;
    s.nstep := s.nstep + 1;
    if not success then s.nfail := s.nfail + 1; end if;
    if success
     then old_theta := theta; old_t := s.sol.t;
          old_solution := s.sol.v;
          w_all_sum := w_all_sum + s.sol.v;
          n_all_sum := n_all_sum + 1;
          tmp := Is_Multiple(theta,twopi*p.maxstep,tol);
          if tmp /= 0
           then w_sum := w_sum + s.sol.v; n_sum := n_sum + 1;
          end if;
          w_c := Is_Multiple(theta,twopi,tol);
          if w_c /= 0
           then if Is_Equal(s.sol.v,start_solution,epslop)
                 then w_sum := w_sum - s.sol.v * CMPLX(0.5);
                      w_all_sum := w_all_sum - s.sol.v * CMPLX(0.5);
                      stop := true;
                 else stop := (w_c >= max_wc);
                end if;
          end if;
     else theta := old_theta;
          s.sol.t := old_t;
          s.sol.v := old_solution;
    end if;
    if not stop
     then Circular_Step_Control(success,p,twopi,step,nsuccess);
    end if;
  end Circular_Management;

-- UPDATE OF PATH DIRECTION :

  procedure Update_Direction
                ( proj : in boolean;
                  freqcnt,defer,r,m,estm,cntm : in out natural;
                  thresm : in natural; er : in out integer;
                  t,target : in double_complex; x : in Vector; 
                  dt,s,logs : in out Float_Vectors.Vector;
                  logx,wvl0,wvl1,wvl2 : in out Float_Vectors_of_Vectors.Vector;
                  v,errv : in out Float_Vectors.Vector;
                  err : in out double_float ) is

  -- DESCRIPTION :
  --   Computes an approximation of the direction of the path.

  -- ON ENTRY :
  --   file       to write intermediate data on, may be omitted;
  --   proj       whether solution vector lies in projective space;
  --   freqcnt    counts the frequency of calls;
  --   defer      only if freqcnt = defer, calculations are done;
  --   r          current order in extrapolation formula;
  --   m          current value for multiplicity;
  --   estm       current estimated for multiplicity;
  --   cntm       number of consecutive times estm has been guessed;
  --   thresm     threshold for augmenting m to estm;
  --   er         order of extrapolator on the errors;
  --   t          current value of continuation parameter;
  --   target     target value of continuation parameter;
  --   x          current solution vector;
  --   dt         vector with distances to target;
  --   s          s-values w.r.t. the current value m;
  --   logs       logarithms of the s-values;
  --   logx       logarithms of the solution vectors;
  --   wvl0       consecutive estimates for previous direction;
  --   wvl1       consecutive estimates for current direction;
  --   wvl2       used as work space for wvl0 and wvl1;
  --   v          current approximate direction of the path;
  --   errv       vector of errors used to estimate m;
  --   err        norm of errv.

  -- ON RETURN :
  --   All in-out variables are updated, provided freqcnt = defer.

  begin
    if freqcnt >= defer
     then
       freqcnt := 0;
       if proj
        then Projective_Update_Direction
               (r,m,estm,cntm,thresm,er,t,target,x,dt,s,logs,logx,v,errv,err);
        else Affine_Update_Direction
               (r,m,estm,cntm,thresm,er,t,target,x,
                dt,s,logs,logx,wvl0,wvl1,wvl2,v,errv,err);
       end if;
      -- defer := defer + 1;  -- that is asking for troubles !
     else
       freqcnt := freqcnt + 1;
    end if;
  end Update_Direction;

  procedure Update_Direction
                ( file : in file_type; proj : in boolean;
                  freqcnt,defer,r,m,estm,cntm : in out natural;
                  thresm : in natural; er : in out integer;
                  t,target : in double_complex; x : in Vector; 
                  dt,s,logs : in out Float_Vectors.Vector;
                  logx,wvl0,wvl1,wvl2 : in out Float_Vectors_of_Vectors.Vector;
                  v,errv : in out Float_Vectors.Vector;
                  err : in out double_float ) is

  -- DESCRIPTION :
  --   Computes an approximation of the direction of the path and produces
  --   intermediate output to the file.

  begin
    if freqcnt >= defer
     then
       freqcnt := 0;
       if proj
        then Projective_Update_Direction
               (file,r,m,estm,cntm,thresm,er,t,target,x,
                     dt,s,logs,logx,v,errv,err);
        else Affine_Update_Direction
               (file,r,m,estm,cntm,thresm,er,t,target,x,
                     dt,s,logs,logx,wvl0,wvl1,wvl2,v,errv,err);
       end if;
      -- defer := defer + 1;  -- asking for troubles !
       put(file,"direction : "); put(file,v); new_line(file);
       put(file,"difference to old direction : "); put(file,err);
       new_line(file);
       put(file,"++ current m : "); put(file,m,1); put(file," ++ "); 
       put(file,cntm,1); put(file," times estimated m : "); put(file,estm,1);
       put(file," ++ threshold : "); put(file,thresm,1); put_line(file," ++");
       new_line(file);
     else
       freqcnt := freqcnt + 1;
    end if;
  end Update_Direction;

-- LINEAR PATH FOLLOWING FOR ONE PATH :

  procedure Linear_Single_Normal_Silent_Continue
              ( s : in out Solu_Info; target : in double_complex;
                tol : in double_float; proj : in boolean;
                p : in Pred_Pars; c : in Corr_Pars ) is
 
    old_t,prev_t : double_complex;
    old_solution,prev_solution,old_v,prev_v,vv : Vector(s.sol.v'range);
    step : double_float := p.maxstep;
    nsuccess,trial : natural := 0;
    success : boolean := true;

    procedure Predictor is new Single_Predictor(Norm,dH,dH);

    procedure Corrector ( s : in out Solu_Info; c : in Corr_Pars ) is

      procedure Affine_Corrector is
        new Affine_Single_Severe_Normal_Silent_Corrector(Norm,H,dH);
      procedure Projective_Corrector is
        new Projective_Single_Severe_Normal_Silent_Corrector(Norm,H,dH);

    begin
      if proj
       then Projective_Corrector(s,c);
       else Affine_Corrector(s,c);
      end if;
    end Corrector;

  begin
    Linear_Single_Initialize(s,p,old_t,prev_t,old_solution,prev_solution);
    while not (At_End(s.sol.t,target,p.dist_target,tol) and success) 
                                           and (s.niter <= c.maxtot) loop

      Predictor(s,p,prev_solution,prev_v,vv,prev_t,target,step,tol,trial);
      Corrector(s,c);
      Linear_Single_Management(s,p,c,old_t,prev_t,old_solution,prev_solution,
                               old_v,prev_v,vv,step,nsuccess,trial,success);
      if Stop(p,s.sol.t,target,step) then return; end if;
    end loop;
    declare
      cp : Corr_Pars := End_Game_Corrector_Parameters(c,p.dist_target,tol);
    begin
      Corrector(s,cp);
    end;
  end Linear_Single_Normal_Silent_Continue;

  procedure Linear_Single_Normal_Reporting_Continue
              ( file : in file_type; s : in out Solu_Info;
                target : in double_complex; tol : in double_float;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is
 
    old_t,prev_t : double_complex;
    old_solution,prev_solution,old_v,prev_v,vv : Vector(s.sol.v'range);
    step : double_float := p.maxstep;
    nsuccess,trial : natural := 0;
    success : boolean := true;

    procedure Predictor is new Single_Predictor(Norm,dH,dH);

    procedure Corrector ( s : in out Solu_Info; c : in Corr_Pars ) is

      procedure Affine_Corrector is
        new Affine_Single_Severe_Normal_Reporting_Corrector(Norm,H,dH);
      procedure Projective_Corrector is
        new Projective_Single_Severe_Normal_Reporting_Corrector(Norm,H,dH);

    begin
      if proj
       then Projective_Corrector(file,s,c);
       else Affine_Corrector(file,s,c);
      end if;
    end Corrector;

  begin
    Linear_Single_Initialize(s,p,old_t,prev_t,old_solution,prev_solution);
    sWrite(file,s.sol.all);
    while not (At_End(s.sol.t,target,p.dist_target,tol) and success) 
                                           and (s.niter <= c.maxtot) loop

      Predictor(s,p,prev_solution,prev_v,vv,prev_t,target,step,tol,trial);
      if p.predictor_type = 6
       then put_line(file,"previous and current direction : "); 
            for i in prev_v'range loop
              put(file,prev_v(i),3,3,3);  put(file,"   ");
              put(file,vv(i),3,3,3);      new_line(file);
            end loop;
      end if;
      pWrite(file,step,s.sol.t,s.sol.all);
      Corrector(s,c);
      sWrite(file,s.sol.all);
      Linear_Single_Management(s,p,c,old_t,prev_t,old_solution,prev_solution,
                               old_v,prev_v,vv,step,nsuccess,trial,success);
      if Stop(p,s.sol.t,target,step) then return; end if;
    end loop;
    declare
      cp : Corr_Pars := End_Game_Corrector_Parameters(c,p.dist_target,tol);
    begin
      Corrector(s,cp);
      sWrite(file,s.sol.all);
    end;
  end Linear_Single_Normal_Reporting_Continue;

  procedure Linear_Single_Conditioned_Silent_Continue
              ( s : in out Solu_Info; target : in double_complex;
                tol : in double_float; proj : in boolean;
                rtoric : in natural; 
                v : in out Float_Vectors.Link_to_Vector;
                errorv : in out double_float;
                p : in Pred_Pars; c : in Corr_Pars ) is

    old_t,prev_t : double_complex;
    old_solution,prev_solution,old_v,prev_v,vv : Vector(s.sol.v'range);
    step : double_float := p.maxstep;
    nsuccess,trial : natural := 0;
    success : boolean := true;

    dls,stp : double_float := 0.0;
    tolsing : constant double_float
            := Continuation_Parameters.tol_endg_inverse_condition;
    diff : Float_Vectors.Vector(s.sol.v'range) := (s.sol.v'range => 0.0);

    r : natural := 0;
    er : integer := 0;
    dt,ds,logs : Float_Vectors.Vector(0..rtoric) := (0..rtoric => 0.0);
    logx : Float_Vectors_of_Vectors.Vector(0..rtoric);
    wvl0,wvl1,wvl2 : Float_Vectors_of_Vectors.Vector(1..rtoric);
    errv : Float_Vectors.Vector(0..rtoric) := (0..rtoric => 0.0);

    m : natural := 1;
    thresm : natural := p.success_steps;
    estm : natural := m;
    fcnt,cntm : natural := 0;
    defer : natural := thresm;

    procedure Predictor is new Single_Predictor(Norm,dH,dH);

    procedure Corrector ( s : in out Solu_Info; c : in Corr_Pars ) is

      procedure Affine_Corrector is
        new Affine_Single_Severe_Conditioned_Silent_Corrector(Norm,H,dH);
      procedure Projective_Corrector is
        new Projective_Single_Severe_Conditioned_Silent_Corrector(Norm,H,dH);

    begin
      if proj
       then Projective_Corrector(s,c);
       else Affine_Corrector(s,c);
      end if;
    end Corrector;

  begin
    Linear_Single_Initialize(s,p,old_t,prev_t,old_solution,prev_solution);
    if rtoric > 0
     then s.sol.m := m; 
    end if;
    while not (At_End(s.sol.t,target,p.dist_target,tol) and success)
                                           and (s.niter <= c.maxtot) loop

      if (rtoric > 0) 
       then
         if success and then s.rcond > tolsing
                    and then (errorv < 100.0) -- avoid divergence
          then Update_Direction
                 (proj,fcnt,defer,r,s.sol.m,estm,cntm,thresm,er,s.sol.t,target,
                  s.sol.v,dt,ds,logs,logx,wvl0,wvl1,wvl2,v.all,errv,errorv);
          else er := -2;
         end if;
      end if;

      Predictor(s,p,prev_solution,prev_v,vv,prev_t,target,step,tol,trial);
      Corrector(s,c);
      Linear_Single_Management(s,p,c,old_t,prev_t,old_solution,prev_solution,
                               old_v,prev_v,vv,step,nsuccess,trial,success);
      if Stop(p,s.sol.t,target,step) then return; end if;
    end loop;
    declare
      cp : Corr_Pars := End_Game_Corrector_Parameters(c,p.dist_target,tol);
    begin
      Corrector(s,cp);
    end;
  end Linear_Single_Conditioned_Silent_Continue;

  procedure Linear_Single_Conditioned_Reporting_Continue
              ( file : in file_type; s : in out Solu_Info;
                target : in double_complex; tol : in double_float;
                proj : in boolean; rtoric : in natural;
                v : in out Float_Vectors.Link_to_Vector;
                errorv : in out double_float;
                p : in Pred_Pars; c : in Corr_Pars ) is

    old_t,prev_t : double_complex;
    step : double_float := p.maxstep;
    nsuccess,trial : natural := 0;
    old_solution,prev_solution,old_v,prev_v,vv : Vector(s.sol.v'range);
    success : boolean := true;

    dls,stp : double_float := 0.0;
    tolsing : constant double_float
            := Continuation_Parameters.tol_endg_inverse_condition;
    diff : Float_Vectors.Vector(s.sol.v'range) := (s.sol.v'range => 0.0);

    r : natural := 0;
    er : integer := 0;
    dt,ds,logs : Float_Vectors.Vector(0..rtoric) := (0..rtoric => 0.0);
    logx : Float_Vectors_of_Vectors.Vector(0..rtoric);
    wvl0,wvl1,wvl2 : Float_Vectors_of_Vectors.Vector(1..rtoric);
    errv : Float_Vectors.Vector(0..rtoric) := (0..rtoric => 0.0);

    m : natural := 1;
    thresm : natural := p.success_steps;
    estm : natural := m;
    fcnt,cntm : natural := 0;
    defer : natural := thresm;

    procedure Predictor is new Single_Predictor(Norm,dH,dH);

    procedure Corrector ( s : in out Solu_Info; c : in Corr_Pars ) is

      procedure Affine_Corrector is
        new Affine_Single_Severe_Conditioned_Reporting_Corrector(Norm,H,dH);
      procedure Projective_Corrector is
        new Projective_Single_Severe_Conditioned_Reporting_Corrector(Norm,H,dH);

    begin
      if proj
       then Projective_Corrector(file,s,c);
       else Affine_Corrector(file,s,c);
      end if;
    end Corrector;

  begin
    Linear_Single_Initialize(s,p,old_t,prev_t,old_solution,prev_solution);
    sWrite(file,s.sol.all);          -- writing the start solution
    if rtoric > 0
     then s.sol.m := m;
    end if;
    while not (At_End(s.sol.t,target,p.dist_target,tol) and success)
                                           and (s.niter <= c.maxtot) loop

      if (rtoric > 0) 
       then
         if success and then s.rcond > tolsing
                    and then (errorv < 100.0) -- avoid divergence
          then Update_Direction(file,
                  proj,fcnt,defer,r,s.sol.m,estm,cntm,thresm,er,s.sol.t,target,
                  s.sol.v,dt,ds,logs,logx,wvl0,wvl1,wvl2,v.all,errv,errorv);
          else er := -2;
         end if;
      end if;

      Predictor(s,p,prev_solution,prev_v,vv,prev_t,target,step,tol,trial);
      pWrite(file,step,s.sol.t,s.sol.all);
      Corrector(s,c);
      sWrite(file,s.sol.all);
      Linear_Single_Management(s,p,c,old_t,prev_t,old_solution,prev_solution,
                               old_v,prev_v,vv,step,nsuccess,trial,success);
      if Stop(p,s.sol.t,target,step) then return; end if;
    end loop;
    declare
      cp : Corr_Pars := End_Game_Corrector_Parameters(c,p.dist_target,tol);
    begin
      Corrector(s,cp);
      sWrite(file,s.sol.all);
    end;
  end Linear_Single_Conditioned_Reporting_Continue;

-- LINEAR PATH FOLLOWING FOR A NUMBER OF PATHS :

  procedure Linear_Multiple_Normal_Silent_Continue
              ( s : in out Solu_Info_Array;
                target : in double_complex; tol,dist_sols : in double_float;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is

    t,old_t,prev_t : double_complex;
    sa,old_sa,prev_sa : Solution_Array(s'range);
    step : double_float := p.maxstep;
    cnt,nsuccess,trial : natural := 0;
    pivot : natural := s'first;
    success,fail : boolean := true;

    procedure Predictor is new Multiple_Predictor(Norm,dH,dH);

    procedure Affine_Corrector is
      new Affine_Multiple_Loose_Normal_Silent_Corrector(Norm,H,dH);
    procedure Projective_Corrector is
      new Projective_Multiple_Loose_Normal_Silent_Corrector(Norm,H,dH);

    procedure Correct ( s : in out Solu_Info_Array;
                        pivot : in out natural; dist_sols : in double_float;
                        c : in Corr_Pars; fail : out boolean ) is
    begin
      Copy(sa,s);
      if proj
       then Projective_Corrector(s,pivot,dist_sols,c,fail);
       else Affine_Corrector(s,pivot,dist_sols,c,fail);
      end if;
    end Correct;

  begin
    Linear_Multiple_Initialize(s,p,t,old_t,prev_t,sa,old_sa,prev_sa);
    while not (At_End(t,target,p.dist_target,tol) and success) 
                            and (s(s'first).niter <= c.maxtot) loop

      Predictor(s,p,sa,prev_sa,t,prev_t,target,step,tol,dist_sols,trial);
      Correct(s,pivot,dist_sols,c,fail); Copy(s,sa);
      success := not fail;
      Linear_Multiple_Management(s,sa,old_sa,prev_sa,t,old_t,prev_t,p,step,
                                 pivot,nsuccess,trial,success);
      if step < p.minstep then return; end if;
    end loop;
    declare
      cp : Corr_Pars := End_Game_Corrector_Parameters(c,p.dist_target,tol);
    begin
      Correct(s,pivot,dist_sols,cp,fail);
    end;
  end Linear_Multiple_Normal_Silent_Continue;

  procedure Linear_Multiple_Normal_Reporting_Continue
              ( file : in file_type; s : in out Solu_Info_Array;
                target : in double_complex; tol,dist_sols : in double_float;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is

    t,old_t,prev_t : double_complex;
    sa,old_sa,prev_sa : Solution_Array(s'range);
    step : double_float := p.maxstep;
    pivot : natural := s'first;
    cnt,nsuccess,trial : natural := 0;
    success,fail : boolean := true;

    procedure Predictor is new Multiple_Predictor(Norm,dH,dH);

    procedure Affine_Corrector is
      new Affine_Multiple_Loose_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Projective_Corrector is
      new Projective_Multiple_Loose_Normal_Reporting_Corrector(Norm,H,dH);

    procedure Correct ( file : in file_type; s : in out Solu_Info_Array;
                        pivot : in out natural; dist_sols : in double_float;
                        c : in Corr_Pars; fail : out boolean ) is
    begin
      Copy(sa,s);
      if proj
       then Projective_Corrector(file,s,pivot,dist_sols,c,fail);
       else Affine_Corrector(file,s,pivot,dist_sols,c,fail);
      end if;
    end Correct;

  begin
    Linear_Multiple_Initialize(s,p,t,old_t,prev_t,sa,old_sa,prev_sa);
    for k in s'range loop                       -- write the start solutions
      sWrite(file,sa(k).all);
    end loop;
    while not (At_End(t,target,p.dist_target,tol) and success) 
                            and (s(s'first).niter <= c.maxtot) loop

      Predictor(s,p,sa,prev_sa,t,prev_t,target,step,tol,dist_sols,trial);
      pWrite(file,step,t);
      Correct(file,s,pivot,dist_sols,c,fail); Copy(s,sa);
      success := not fail;
      Linear_Multiple_Management(s,sa,old_sa,prev_sa,t,old_t,prev_t,p,step,
                                 pivot,nsuccess,trial,success);
      if step < p.minstep then return; end if;
    end loop;
    declare
      cp : Corr_Pars := End_Game_Corrector_Parameters(c,p.dist_target,tol);
    begin
      Correct(file,s,pivot,dist_sols,cp,fail);
    end;
  end Linear_Multiple_Normal_Reporting_Continue;

  procedure Linear_Multiple_Conditioned_Silent_Continue
              ( s : in out Solu_Info_Array;
                target : in double_complex; tol,dist_sols : in double_float;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is

    t,old_t,prev_t : double_complex;
    sa,old_sa,prev_sa : Solution_Array(s'range);
    step : double_float := p.maxstep;
    pivot : natural := s'first;
    cnt,nsuccess,trial : natural := 0;
    success,fail : boolean := true;

    procedure Predictor is new Multiple_Predictor(Norm,dH,dH);

    procedure Affine_Corrector is
      new Affine_Multiple_Loose_Conditioned_Silent_Corrector(Norm,H,dH);
    procedure Projective_Corrector is
      new Projective_Multiple_Loose_Conditioned_Silent_Corrector(Norm,H,dH);

    procedure Correct ( s : in out Solu_Info_Array;
                        pivot : in out natural; dist_sols : in double_float;
                        c : in Corr_Pars; fail : out boolean ) is
    begin
      Copy(sa,s);
      if proj
       then Projective_Corrector(s,pivot,dist_sols,c,fail);
       else Affine_Corrector(s,pivot,dist_sols,c,fail);
      end if;
    end Correct;

  begin
    Linear_Multiple_Initialize(s,p,t,old_t,prev_t,sa,old_sa,prev_sa);
    while not (At_End(t,target,p.dist_target,tol) and success) 
                            and (s(s'first).niter <= c.maxtot) loop

      Predictor(s,p,sa,prev_sa,t,prev_t,target,step,tol,dist_sols,trial);
      Correct(s,pivot,dist_sols,c,fail); Copy(s,sa);
      success := not fail;
      Linear_Multiple_Management(s,sa,old_sa,prev_sa,t,old_t,prev_t,p,step,
                                 pivot,nsuccess,trial,success);
      if step < p.minstep then return; end if;
    end loop;
    declare
      cp : Corr_Pars := End_Game_Corrector_Parameters(c,p.dist_target,tol);
    begin
      Correct(s,pivot,dist_sols,cp,fail);
    end;
  end Linear_Multiple_Conditioned_Silent_Continue;

  procedure Linear_Multiple_Conditioned_Reporting_Continue
              ( file : in file_type; s : in out Solu_Info_Array;
                target : in double_complex; tol,dist_sols : in double_float;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is

    t,old_t,prev_t : double_complex;
    sa,old_sa,prev_sa : Solution_Array(s'range);
    step : double_float := p.maxstep;
    pivot : natural := s'first;
    cnt,nsuccess,trial : natural := 0;
    success,fail : boolean := true;

    procedure Predictor is new Multiple_Predictor(Norm,dH,dH);

    procedure Affine_Corrector is 
      new Affine_Multiple_Loose_Conditioned_Reporting_Corrector(Norm,H,dH);
    procedure Projective_Corrector is 
      new Projective_Multiple_Loose_Conditioned_Reporting_Corrector(Norm,H,dH);

    procedure Correct ( file : in file_type; s : in out Solu_Info_Array;
                        pivot : in out natural; dist_sols : in double_float;
                        c : in Corr_Pars; fail : out boolean ) is
    begin
      Copy(sa,s);
      if proj
       then Projective_Corrector(file,s,pivot,dist_sols,c,fail);
       else Affine_Corrector(file,s,pivot,dist_sols,c,fail);
      end if;
    end Correct;

  begin
    Linear_Multiple_Initialize(s,p,t,old_t,prev_t,sa,old_sa,prev_sa);
    for k in s'range loop                        -- write start solutions
      sWrite(file,sa(k).all);
    end loop;
    while not (At_End(t,target,p.dist_target,tol) and success) 
                            and (s(s'first).niter <= c.maxtot) loop

      Predictor(s,p,sa,prev_sa,t,prev_t,target,step,tol,dist_sols,trial);
      pWrite(file,step,t);
      Correct(file,s,pivot,dist_sols,c,fail); Copy(s,sa);
      success := not fail;
      Linear_Multiple_Management(s,sa,old_sa,prev_sa,t,old_t,prev_t,p,step,
                                 pivot,nsuccess,trial,success);
      if step < p.minstep then return; end if;
    end loop;
    declare
      cp : Corr_Pars := End_Game_Corrector_Parameters(c,p.dist_target,tol);
    begin
      Correct(file,s,pivot,dist_sols,cp,fail);
    end;
  end Linear_Multiple_Conditioned_Reporting_Continue;

-- CIRCULAR PATH FOLLOWING FOR ONE PATH :

  procedure Circular_Single_Normal_Reporting_Continue
              ( file : in file_type; s : in out Solu_Info;
                target : in double_complex; tol,epslop : in double_float;
                wc : out natural; max_wc : in natural; sum,all_sum : out Vector;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is
 
    old_t,prev_t : double_complex;
    t0_min_target : double_complex := s.sol.t - target;
    theta,old_theta : double_float := 0.0;
    twopi : constant double_float := 2.0*PI;
    step : double_float := twopi*p.maxstep;
    old_solution,prev_solution,start_solution : Vector(s.sol.v'range);
    w_c,nsuccess : natural := 0;
    success : boolean := true;
    stop : boolean := false;
    w_sum,w_all_sum : Vector(s.sol.v'range) := CMPLX(0.5)*s.sol.v;
    n_sum,n_all_sum : natural := 0;

    procedure T_C_Predictor is new Tangent_Circular_Predictor(Norm,dH,dH);

    procedure Affine_Corrector is
      new Affine_Single_Severe_Normal_Reporting_Corrector(Norm,H,dH);
    procedure Projective_Corrector is
      new Projective_Single_Severe_Normal_Reporting_Corrector(Norm,H,dH);

  begin
    old_t := s.sol.t; old_solution := s.sol.v;        -- INITIALIZATION
    start_solution := s.sol.v;
    if p.predictor_type = 0 
     then prev_t := s.sol.t; prev_solution := s.sol.v;
    end if;
    sWrite(file,s.sol.all);                -- write the start solution
    while (s.niter <= c.maxtot) loop 
      case p.predictor_type is
        when 0 => Secant_Circular_Predictor(s.sol.v,prev_solution,s.sol.t,theta,
                                    prev_t,t0_min_target,target,step,tol);
        when 2 => T_C_Predictor(s.sol.v,s.sol.t,theta, t0_min_target,target,
                                step,tol);
                  s.nsyst := s.nsyst + 1;
        when others => null; -- these choices make no sense !!!
      end case;
      pWrite(file,step,s.sol.t,s.sol.all);
      if proj
       then Projective_Corrector(file,s,c);
       else Affine_Corrector(file,s,c);
      end if;
      sWrite(file,s.sol.all);
      Circular_Management
           (s,p,c,old_t,prev_t,start_solution,old_solution,prev_solution,
            w_sum,w_all_sum,twopi,epslop,tol,theta,old_theta,
            step,nsuccess,n_sum,n_all_sum,w_c,max_wc,stop,success);
      exit when stop;
      if step < p.minstep then return; end if;
    end loop;
    wc := w_c;
    if n_all_sum /= 0
     then all_sum := w_all_sum*CMPLX(1.0/double_float(n_all_sum));
    end if;
    if n_sum /= 0
     then sum := w_sum*CMPLX(1.0/double_float(n_sum));
     elsif n_all_sum /= 0
         then all_sum := w_all_sum*CMPLX(1.0/double_float(n_all_sum));
    end if;
  end Circular_Single_Normal_Reporting_Continue;

  procedure Circular_Single_Conditioned_Reporting_Continue
              ( file : in file_type; s : in out Solu_Info;
                target : in double_complex; tol,epslop : in double_float;
                wc : out natural; max_wc : in natural; sum,all_sum : out Vector;
                proj : in boolean; p : in Pred_Pars; c : in Corr_Pars ) is
 
    old_t,prev_t : double_complex;
    theta,old_theta : double_float := 0.0;
    twopi : constant double_float := 2.0*PI;
    step : double_float := twopi*p.maxstep;
    t0_min_target : double_complex := s.sol.t - target;
    old_solution,prev_solution,start_solution : Vector(s.sol.v'range);
    w_c,nsuccess : natural := 0;
    success : boolean := true;
    stop : boolean := false;
    w_sum,w_all_sum : Vector(s.sol.v'range) := CMPLX(0.5)*s.sol.v;
    n_sum,n_all_sum : natural := 0;

    procedure T_C_Predictor is new Tangent_Circular_Predictor(Norm,dH,dH);

    procedure Affine_Corrector is
      new Affine_Single_Severe_Conditioned_Reporting_Corrector(Norm,H,dH);
    procedure Projective_Corrector is
      new Projective_Single_Severe_Conditioned_Reporting_Corrector(Norm,H,dH);

  begin
    old_t := s.sol.t; old_solution := s.sol.v;            -- INITIALIZATION
    start_solution := s.sol.v;
    if p.predictor_type = 0
     then prev_t := s.sol.t; prev_solution := old_solution;
    end if;
    sWrite(file,s.sol.all);              -- writing the start solution
    while s.niter <= c.maxtot loop
      case p.predictor_type is
        when 0 => Secant_Circular_Predictor(s.sol.v,prev_solution,s.sol.t,theta,
                                    prev_t,t0_min_target,target,step,tol);
        when 2 => T_C_Predictor(s.sol.v,s.sol.t,theta,t0_min_target,
                                target,step,tol);
                  s.nsyst := s.nsyst + 1;
        when others => null; -- these choices make no sense !!!
      end case;
      pWrite(file,step,s.sol.t,s.sol.all);
      if proj
       then Projective_Corrector(file,s,c);
       else Affine_Corrector(file,s,c);
      end if;
      sWrite(file,s.sol.all);
      Circular_Management
          (s,p,c,old_t,prev_t,start_solution,old_solution,prev_solution,
           w_sum,w_all_sum,twopi,epslop,tol,theta,old_theta,
           step,nsuccess,n_sum,n_all_sum,w_c,max_wc,stop,success);
      exit when stop;
      if step < p.minstep then return; end if;
    end loop;
    wc := w_c;
    if n_all_sum /= 0
     then all_sum := w_all_sum*CMPLX(1.0/double_float(n_all_sum));
    end if;
    if n_sum /= 0
     then sum := w_sum*CMPLX(1.0/double_float(n_sum));
     elsif n_all_sum /= 0
         then sum := w_all_sum*CMPLX(1.0/double_float(n_all_sum));
    end if;
  end Circular_Single_Conditioned_Reporting_Continue;

end Path_Trackers;
