with Predictors;                   use Predictors;

package body Dispatch_Predictors is

  procedure Single_Predictor
              ( s : in out Solu_Info; p : in Pred_Pars;
                prev_x,prev_v : in Vector; v : in out Vector;
                prev_t,target : in double_complex;
                step,tol : in double_float; trial : in out natural ) is

    procedure TR_Predictor is new Tangent_Single_Real_Predictor(Norm,dH,dH);
    procedure TC_Predictor is new Tangent_Single_Complex_Predictor(Norm,dH,dH);
    procedure TG_Predictor is new Tangent_Geometric_Predictor(Norm,dH,dH);
    procedure HR_Predictor is new Hermite_Single_Real_Predictor(Norm,dH,dH);

  begin
    case p.predictor_type is
      when 0 => Secant_Single_Real_Predictor
                  (s.sol.v,prev_x,s.sol.t,prev_t,target,step,tol,p.power);
      when 1 => Secant_Single_Complex_Predictor
                  (s.sol.v,prev_x,s.sol.t,prev_t,target,step,tol,
                   p.dist_target,trial);
      when 2 => Secant_Geometric_Predictor
                  (s.sol.v,prev_x,s.sol.t,prev_t,target,step,tol);
      when 3 => TR_Predictor(s.sol.v,s.sol.t,target,step,tol,p.power);
                s.nsyst := s.nsyst + 1;
      when 4 => TC_Predictor
                  (s.sol.v,s.sol.t,target,step,tol,p.dist_target,trial);
                s.nsyst := s.nsyst + 1;
      when 5 => TG_Predictor(s.sol.v,s.sol.t,target,step,tol);
                s.nsyst := s.nsyst + 1;
      when 6 => HR_Predictor
                  (s.sol.v,prev_x,s.sol.t,prev_t,target,v,prev_v,step,tol);
                s.nsyst := s.nsyst + 1;
      when others => null;
    end case;
  end Single_Predictor;

  procedure Multiple_Predictor
              ( s : in out Solu_Info_Array; p : in Pred_Pars;
                sa : in out Solution_Array; prev_sa : in Solution_Array;
                t : in out double_complex; prev_t,target : in double_complex;
                step,tol,dist : in double_float; trial : in natural ) is

    cnt : natural := 0;

    procedure TR_Predictor is new Tangent_Multiple_Real_Predictor(Norm,dH,dH);
    procedure TC_Predictor is
      new Tangent_Multiple_Complex_Predictor(Norm,dH,dH);

  begin
    case p.predictor_type is
      when 0 => Secant_Multiple_Real_Predictor
                  (sa,prev_sa,t,prev_t,target,step,tol,dist,p.power);
      when 1 => Secant_Multiple_Complex_Predictor
                  (sa,prev_sa,t,prev_t,target,step,tol,dist,
                   p.dist_target,trial);
      when 3 => TR_Predictor(sa,t,target,step,tol,dist,cnt,p.power);
                for k in s'range loop
                  s(k).nsyst := s(k).nsyst + 1;
                end loop;
      when 4 => TC_Predictor
                  (sa,t,target,step,tol,dist,p.dist_target,trial,cnt);
                for k in s'range loop
                  s(k).nsyst := s(k).nsyst + 1;
                end loop;
      when others => null;
    end case;
  end Multiple_Predictor;

end Dispatch_Predictors;
