with Floating_Point_Numbers;    use Floating_Point_Numbers;
with Solutions;                 use Solutions;

package Continuation_Data is

-- DESCRIPTION :
--   In order to keep the parameters and variables manageble,
--   they have been grouped into records.

-- DATA STRUCTURES FOR PARAMETERS :

  type Pred_Pars is record    -- contains the parameters for the predictor

    minstep,maxstep : double_float;  -- minimum and maximum step size
    expfac,redfac : double_float;    -- expansion and reduction factor
                                     --  for step length control
    success_steps : natural;     -- number of successful steps before expansion
    predictor_type : natural;    -- type of predictor used
    dist_target : double_float;  -- distance to target
    power : positive;            -- power of t in (polyhedral) homotopy

  end record;

  type Corr_Pars is record    -- contains the parameters for the corrector

    epsrx,epsax,epsrf,epsaf : double_float;  
                              -- desired precisions for x and its residual f(x)
                              -- once relative (r) and once absolute (a)

    maxit,maxtot : natural;   -- maximum number of corrector iterations
                              -- for one step and for the whole path
  end record;

-- DATASTRUCTURES FOR VARIABLES :

  type Solu_Info is record    -- contains information about the solution

    sol : Link_to_Solution;   -- the solution: vector, t and multiplicity

    corr,cora,resr,resa,rcond : double_float; 
                              -- last correction (cor) and residual (res), 
                              -- once relative (r) and once absolute (a)
                              -- and estimate for inverse condition of jacobian

    length_path : double_float;  -- length of the path

    nstep,nfail,niter,nsyst : natural;  -- various counters :
                              -- number of steps, failures, corrector
                              -- iterations and number of linear systems solved
  end record;

  type Solu_Info_Array is array ( integer range <> ) of Solu_Info;

-- CREATERS :

  function Shallow_Create ( s : Link_to_Solution ) return Solu_Info;
  function Deep_Create    ( s : Solution ) return Solu_Info;
  function Shallow_Create ( s : Solution_Array ) return Solu_Info_Array;
  function Deep_Create    ( s : Solution_Array ) return Solu_Info_Array;
  function Shallow_Create ( s : Solution_List )  return Solu_Info_Array;
  function Deep_Create    ( s : Solution_List )  return Solu_Info_Array;

  function Shallow_Create ( s : Solu_Info ) return Link_to_Solution;
  function Deep_Create    ( s : Solu_Info ) return Solution;
  function Shallow_Create ( s : Solu_Info_Array ) return Solution_Array;
  function Deep_Create    ( s : Solu_Info_Array ) return Solution_Array;
  function Shallow_Create ( s : Solu_Info_Array ) return Solution_List;
  function Deep_Create    ( s : Solu_Info_Array ) return Solution_List;

  -- DESCRIPTION :
  --   A shallow create copies the pointer to the solution, while
  --   a deep create allocates memory for a copy of the solution.

-- OPERATIONS ON Solu_Info :

  procedure Copy_Info ( s1 : in Solu_Info; s2 : in out Solu_Info );
  procedure Copy_Solu ( s1 : in Solu_Info; s2 : in out Solu_Info );
  procedure Copy      ( s1 : in Solu_Info; s2 : in out Solu_Info );

  -- DESCRIPTION : 
  --   Copies the information, the solution or everything from s1 to s2.

  procedure Init_Info ( s : in out Solu_Info );

  -- DESCRIPTION :
  --   Initializes the information of the solution.

  procedure Add_Info ( s1 : in out Solu_Info; s2 : in Solu_Info );

  -- DESCRIPTION :
  --   Adds the information in the counters of s2 to s1.

  procedure Update_Info ( s1 : in out Solu_Info; s2 : in Solu_Info );

  -- DESCRIPTION :
  --   Adds the information in the counters of s2 to s1 and copies the
  --   other information from s2 to s1.

-- OPERATIONS ON Solu_Info_Array :

  procedure Copy ( s : in Solu_Info_Array; sa : in out Solution_Array );
  procedure Copy ( sa : in Solution_Array; s : in out Solu_Info_Array );

  -- DESCRIPTION : Copies s into sa or vice versa.

-- DESTRUCTORS :

  procedure Clear ( s : in out Solu_Info );
  procedure Clear ( s : in out Solu_Info_Array );

  -- DESCRIPTION :
  --   This is clear is only needed after a deep create.

end Continuation_Data;
